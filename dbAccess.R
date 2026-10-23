library(DBI)
library(RSQLite)

dp_connect <- function(path = "drugprot.sqlite") {
  con <- dbConnect(RSQLite::SQLite(), path, flags = SQLITE_RO)  # read-only
  dbExecute(con, "PRAGMA query_only = ON;")
  dbExecute(con, "PRAGMA cache_size = -100000;")  # ~100 MB page cache per session
  con
}

## ---- metadata (load once at startup, replaces several .RData) -------------
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
    nProtein         = mv[["nProtein"]],
    nTreatment       = mv[["nTreatment"]],
    nTimes           = mv[["nTimes"]]
  )
}

.in_clause <- function(ids) paste(ids, collapse = ",")

## Drug p-values and coefficients exist as direct effects (conditional on the
## lagged proteome) and total effects (drug design alone). At 6 h both coincide.
## match.arg also keeps the value out of the SQL string unless it is valid.
.effect_type <- function(effect_type) match.arg(effect_type, c("direct", "total"))

## ===========================================================================
## DRUG EFFECTS
## ===========================================================================
dp_drug_selPvecs <- function(con, selection, nTreatment, nTimes,
                             effect_type = "direct") {
  raw <- dbGetQuery(con, sprintf(
    "SELECT protein, treatment, time_idx, pvalue
       FROM drug_pvalue WHERE effect_type = '%s' AND protein IN (%s)",
    .effect_type(effect_type), .in_clause(selection)))
  selPvecs <- array(NA_real_, dim = c(nTreatment, nTimes, length(selection)))
  pcol <- match(raw$protein, selection)
  selPvecs[cbind(raw$treatment, raw$time_idx, pcol)] <- raw$pvalue
  selPvecs
}

## ===========================================================================
## PROTEIN NETWORK
## ===========================================================================
dp_net_pvals <- function(con, selection, alpha, nProtein, nNetTimes) {
  inSel <- .in_clause(selection)
  lapply(seq_len(nNetTimes), function(tp) {
    ## (source IN ... OR target IN ...) triggers SQLite's MULTI-INDEX OR over
    ## idx_net_source / idx_net_target -- no scan of the ~58M-row table.
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
## DRUG COEFFICIENTS
## ===========================================================================
## Drug coefficients (mean debiased effect) of all treatments for a protein
## selection, in ONE query. Returns data.frame(treatment, protein, time_idx, coef).
dp_drug_coef_sel <- function(con, selection, effect_type = "direct") {
  dbGetQuery(con, sprintf(
    "SELECT treatment, protein, time_idx, coef FROM drug_coef
      WHERE effect_type = '%s' AND protein IN (%s)",
    .effect_type(effect_type), .in_clause(selection)))
}

dp_drug_dEff <- function(con, treatment_idx, nProtein, nTimes,
                         effect_type = "direct") {
  raw <- dbGetQuery(con, sprintf(
    "SELECT protein, time_idx, coef FROM drug_coef
      WHERE effect_type = '%s' AND treatment = %d",
    .effect_type(effect_type), treatment_idx))
  dEff <- matrix(NA_real_, nrow = nProtein, ncol = nTimes)
  if (nrow(raw)) dEff[cbind(raw$protein, raw$time_idx)] <- raw$coef
  dEff
}

## ===========================================================================
## PROTEIN COEFFICIENTS
## ===========================================================================
dp_protein_bhat <- function(con, target, time_idx, nProtein) {
  raw <- dbGetQuery(con, sprintf(
    "SELECT source, coef FROM protein_coef WHERE target = %d AND time_idx = %d",
    target, time_idx))
  bhat <- numeric(nProtein)           # zeros by default
  if (nrow(raw)) bhat[raw$source] <- raw$coef
  bhat
}
