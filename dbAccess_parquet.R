## dbAccess_parquet.R
## --------------------------------------------------------------------------
## Drop-in replacement for dbAccess.R that reads the Parquet files (produced by
## sqliteToParquet.R) via DuckDB, instead of the SQLite database.
##
##   source("dbAccess_parquet.R")
##   db <- dp_connect()            # DuckDB connection; Parquet stays on disk
##   meta <- dp_metadata(db)
##
## The five helper names and signatures are IDENTICAL to dbAccess.R, so app.R
## needs only to change which file it sources -- nothing else. DuckDB reads the
## Parquet files directly with predicate + column pushdown, so a query touches
## only the row groups and columns it needs; nothing is loaded whole.
##
## SAME DESIGN PRINCIPLE as the SQLite version: each helper returns ONLY the
## storage+selection object the original reactive consumed. All statistics
## (p.adjust, alpha filtering, the explicit n=, min-over-time, edge-sign) stay
## in the reactives, unchanged. `alpha` in dp_net_pvals is the same performance
## pre-filter as the original, not a statistical decision.
## --------------------------------------------------------------------------

library(DBI)
library(duckdb)

## ---- connection -----------------------------------------------------------
## Returns a DuckDB connection with one read_parquet-backed VIEW per table, so
## the query bodies can reference plain table names (proteins, drug_pvalue, ...)
## exactly like the SQLite version. Views are lazy: no data is read until a
## query runs, and then only the needed row groups/columns.
dp_connect <- function(parquet_dir = "data/parquet",
                       in_memory_network = TRUE,
                       memory_limit = "4GB") {
  con <- dbConnect(duckdb::duckdb())           # in-memory engine; data on disk
  ## Ceiling must leave room for the resident network table (~1.5-2 GB) plus
  ## query working memory. Lower it if you set in_memory_network = FALSE.
  dbExecute(con, sprintf("SET memory_limit = '%s';", memory_limit))
  dbExecute(con, "SET threads = 2;")

  ## Small tables stay as lazy VIEWs over the Parquet files (read on demand,
  ## negligible cost). The large network table is the one queried repeatedly on
  ## every links refresh, so we MATERIALIZE it into an in-memory table once at
  ## startup: that keeps it resident instead of re-reading the Parquet file per
  ## query, which is the actual speed win (DuckDB does NOT use indexes for the
  ## `source IN (...) OR target IN (...)` shape -- it vectorizes a scan -- so no
  ## index is created; it would only cost build time and RAM for no benefit).
  ##
  ## Set in_memory_network = FALSE to keep protein_pvalue lazy on disk (lower
  ## RAM, slower links). The ~58M-row network table costs ~1.5-2 GB resident.
  view_tables  <- c("proteins", "treatments", "exp_times", "drug_order", "meta",
                    "drug_lookup", "drug_pvalue", "drug_coef", "protein_coef")
  for (t in view_tables) {
    f <- file.path(parquet_dir, paste0(t, ".parquet"))
    if (!file.exists(f)) stop("Missing Parquet file: ", f)
    dbExecute(con, sprintf(
      "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_parquet('%s');",
      t, normalizePath(f, winslash = "/")))
  }

  f_net <- file.path(parquet_dir, "protein_pvalue.parquet")
  if (!file.exists(f_net)) stop("Missing Parquet file: ", f_net)
  if (isTRUE(in_memory_network)) {
    ## one full read at startup -> resident in-memory table
    dbExecute(con, sprintf(
      "CREATE OR REPLACE TABLE protein_pvalue AS SELECT * FROM read_parquet('%s');",
      normalizePath(f_net, winslash = "/")))
  } else {
    dbExecute(con, sprintf(
      "CREATE OR REPLACE VIEW protein_pvalue AS SELECT * FROM read_parquet('%s');",
      normalizePath(f_net, winslash = "/")))
  }
  con
}

## Close with dp_disconnect(db) (shuts down the DuckDB instance too).
dp_disconnect <- function(con) {
  dbDisconnect(con, shutdown = TRUE)
}

## ---- metadata (load once at startup) --------------------------------------
dp_metadata <- function(con) {
  proteins   <- dbGetQuery(con, "SELECT protein, name FROM proteins ORDER BY protein")
  treatments <- dbGetQuery(con, "SELECT treatment, raw_label, display FROM treatments ORDER BY treatment")
  exp_times  <- dbGetQuery(con, "SELECT idx, hours FROM exp_times ORDER BY idx")$hours
  drug_order <- dbGetQuery(con, "SELECT ord, raw_label FROM drug_order ORDER BY ord")$raw_label
  meta       <- dbGetQuery(con, "SELECT key, value FROM meta")
  mv <- setNames(meta$value, meta$key)
  list(
    prot_names_short = setNames(proteins$name, NULL),
    treatment        = treatments$raw_label,
    treatNames       = setNames(treatments$raw_label, treatments$display),
    expTimes         = exp_times,
    drugOrder        = drug_order,
    nProtein         = as.integer(mv[["nProtein"]]),
    nTreatment       = as.integer(mv[["nTreatment"]]),
    nTimes           = as.integer(mv[["nTimes"]])
  )
}

.in_clause <- function(ids) paste(ids, collapse = ",")

## ===========================================================================
## DRUG EFFECTS -- replaces the load + selPvecs lines of pvec().
## Returns the [nTreatment x nTimes x |selection|] array of RAW drug p-values,
## identical to allPvecs[, , selection].
## ===========================================================================
dp_drug_selPvecs <- function(con, selection, nTreatment, nTimes) {
  raw <- dbGetQuery(con, sprintf(
    "SELECT protein, treatment, time_idx, pvalue
       FROM drug_pvalue WHERE protein IN (%s)", .in_clause(selection)))
  selPvecs <- array(NA_real_, dim = c(nTreatment, nTimes, length(selection)))
  pcol <- match(raw$protein, selection)
  selPvecs[cbind(raw$treatment, raw$time_idx, pcol)] <- raw$pvalue
  selPvecs
}

## ===========================================================================
## PROTEIN NETWORK -- replaces the per-transition selection block of Links_all().
## Returns a list (one data.frame per transition) with columns source, target,
## pvalue: the `Pval_sel` object, RAW p-values pruned at raw p < alpha.
##
## Queries the in-memory `protein_pvalue` table created in dp_connect() (when
## in_memory_network = TRUE). DuckDB vectorizes a scan with the predicate
## applied per row; the table being resident -- not re-read from Parquet per
## call -- is what makes this fast and flat across alpha. The result is
## identical whether the table is in memory or a lazy Parquet view.
## ===========================================================================
dp_net_pvals <- function(con, selection, alpha, nProtein, nNetTimes) {
  inSel <- .in_clause(selection)
  lapply(seq_len(nNetTimes), function(tp) {
    dbGetQuery(con, sprintf(
      "SELECT source, target, pvalue
         FROM protein_pvalue
        WHERE time_idx = %d
          AND (source IN (%s) OR target IN (%s))
          AND pvalue < %.17g",
      tp, inSel, inSel, alpha))
  })
}

## ===========================================================================
## DRUG COEFFICIENTS -- replaces the load + dEff slice in output$DrugEffects.
## Returns the dEff matrix [nProtein x nTimes] for one treatment; absent cells NA.
## ===========================================================================
dp_drug_dEff <- function(con, treatment_idx, nProtein, nTimes) {
  raw <- dbGetQuery(con, sprintf(
    "SELECT protein, time_idx, coef FROM drug_coef WHERE treatment = %d",
    treatment_idx))
  dEff <- matrix(NA_real_, nrow = nProtein, ncol = nTimes)
  if (nrow(raw)) dEff[cbind(raw$protein, raw$time_idx)] <- raw$coef
  dEff
}

## ===========================================================================
## PROTEIN COEFFICIENTS -- replaces the load + bhat lookup in TempGraph().
## Returns the full bhat vector (length nProtein) for one (target, time_idx);
## absent entries 0. time_idx is 2..nTimes, matching expTimes[t+1].
## ===========================================================================
dp_protein_bhat <- function(con, target, time_idx, nProtein) {
  raw <- dbGetQuery(con, sprintf(
    "SELECT source, coef FROM protein_coef WHERE target = %d AND time_idx = %d",
    target, time_idx))
  bhat <- numeric(nProtein)
  if (nrow(raw)) bhat[raw$source] <- raw$coef
  bhat
}

## All protein coefficients for a set of target proteins at one time, in ONE
## query. Returns a data.frame(target, source, coef) -- only the stored
## (positive) edges; callers treat absent (target,source) pairs as 0.
dp_protein_coef_bulk <- function(con, targets, time_idx) {
  if (length(targets) == 0)
    return(data.frame(target = integer(0), source = integer(0), coef = numeric(0)))
  dbGetQuery(con, sprintf(
    "SELECT target, source, coef FROM protein_coef
       WHERE time_idx = %d AND target IN (%s)",
    time_idx, paste(targets, collapse = ",")))
}
