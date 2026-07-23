# Drug-Prot

An interactive R Shiny application for querying statistical evidence of drug
effects, drug–drug interactions, and directed temporal protein dependencies in
a large-scale perturbation proteomics dataset.

**Live application:** <https://ulme.shinyapps.io/DrugProt/>

Drug-Prot lets you specify a set of proteins of interest and returns, for that
set: corrected p-values for the effect of each of 63 single drugs and 59 drug
combinations on protein expression at 6, 24, and 48 hours post-treatment; and a
directed temporal dependency network linking the queried proteins to any other
protein in the dataset. Because multiple-testing correction is applied only
over the queried set, focused analyses gain substantial power relative to
proteome-wide searches.

All statistical evidence (approximately 60 million p-values across 5,392
proteins, 122 treatments, and three time points) is precomputed, so queries are
near-instantaneous and the underlying proteomic dataset is never needed at
query time.

## Table of contents

- [What the app provides](#what-the-app-provides)
- [Repository contents](#repository-contents)
- [Requirements](#requirements)
- [Data files](#data-files)
- [Installation and running locally](#installation-and-running-locally)
- [Using the application](#using-the-application)
- [Data model](#data-model)
- [Storage backends](#storage-backends)
- [Deployment](#deployment)
- [Reproducing the evidence](#reproducing-the-evidence)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## What the app provides

**Drug effects.** For the selected protein set, the minimum corrected p-value
over the set and chosen time points for each single drug and each
experimentally available drug pair, shown as an interactive 122 × 122 heatmap
(diagonal: single-drug evidence; off-diagonal: pairwise interaction evidence).
Drug pairs without experimental data are coded with a sentinel value of 2 and
rendered in black. Companion tables list significant single drugs and
significant drug pairs in order of increasing p-value.

**Protein network.** A directed temporal dependency subgraph linking the
queried proteins to other proteins in the dataset. Two views are provided: a
temporal graph, in which each node is a (protein, time point) pair and edges run
from t = 6 to t = 24 and from t = 24 to t = 48; and a summary graph collapsing
the temporal graph across time. Both are rendered as zoomable, draggable,
force-directed networks. Edge direction reflects the temporal ordering of the
measurements and is not protected against unmeasured confounding.

**Downloads.** Corrected drug-effect p-values and protein–protein edge p-values
for the queried set as CSV, and the summary and temporal graphs as standalone
interactive HTML.

## Repository contents

| Path | Description |
|---|---|
| `app.R` | The Shiny application: UI definition and server logic. |
| `dbAccess_parquet.R` | Data access layer backed by Parquet files read through DuckDB. This is the backend the deployed app uses. |
| `dbAccess.R` | Alternative data access layer backed by a single SQLite database. Exposes the same five helper functions with identical signatures, so it is a drop-in replacement. |
| `examples/exampleSet.txt` | Example protein-set upload file (nuclear envelope / lamina-associated proteins), one HGNC symbol per line. |
| `www/drugprot-logo.svg` | Application logo. |
| `DrugProt.Rproj` | RStudio project file. |
| `LICENSE` | BSD 3-Clause License. |

Note that the data files themselves are **not** in this repository — see
[Data files](#data-files) below.

## Requirements

- R (≥ 4.1 recommended)
- R packages:
  - `shiny`, `shinydashboard`, `shinyWidgets`
  - `plotly`, `networkD3`, `grid`
  - `readxl`
  - `DBI`, `duckdb` (Parquet backend)
  - `RSQLite` (only if using the SQLite backend)

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyWidgets",
  "plotly", "networkD3", "readxl",
  "DBI", "duckdb"
))
# only for the SQLite backend:
# install.packages("RSQLite")
```

## Data files

The precomputed statistical evidence is distributed separately from the code
because of its size. Download the archive from:

<https://doi.org/10.5281/zenodo.21508013>

The archive contains raw, uncorrected p-values and effect estimates for every
protein, drug, and time point, together with a README describing the columns.

Place the Parquet files in a directory named `parquet/` at the root of this
repository (alongside `app.R`), so that the layout is:

```
DrugProt/
├── app.R
├── dbAccess_parquet.R
├── parquet/
│   ├── proteins.parquet
│   ├── treatments.parquet
│   ├── exp_times.parquet
│   ├── drug_order.parquet
│   ├── meta.parquet
│   ├── drug_lookup.parquet
│   ├── drug_pvalue.parquet
│   ├── drug_coef.parquet
│   ├── protein_coef.parquet
│   └── protein_pvalue.parquet
└── ...
```

`protein_pvalue.parquet` is by far the largest file (~58 million rows); the
others are small.

## Installation and running locally

```bash
git clone https://github.com/markusul/DrugProt.git
cd DrugProt
# download and unpack the Parquet data into ./parquet/ (see above)
```

Then, from R in the repository root:

```r
shiny::runApp()
```

The application opens with the **About** tab, which describes what Drug-Prot
does and how to interpret the results. Go to **Settings** to select proteins.

### Memory configuration

The connection is opened in `app.R` with:

```r
db <- dp_connect("parquet", in_memory_network = FALSE, memory_limit = "3GB")
```

- `in_memory_network = FALSE` keeps the ~58M-row network table as a lazy
  Parquet view on disk. This minimises resident memory and is the setting used
  for the hosted deployment.
- `in_memory_network = TRUE` materialises that table into memory once at
  startup (~1.5–2 GB resident), which makes repeated network queries
  noticeably faster. Use this if you have the RAM available locally, and raise
  `memory_limit` accordingly.

Results are identical either way.

## Using the application

**Select proteins** (Settings tab), in one of three ways:

1. Type or select from the searchable picker of all 5,392 measured proteins.
2. Click an example button — *LMNA & RAB7A* (the two most IC50-predictive
   proteins from the paper) or *26 IC50 proteins* (the full predictive set).
3. Upload a plain-text file with one protein identifier per line. See
   `examples/exampleSet.txt` for the expected format.

Proteins are named by HGNC gene symbol (e.g. `RAB7A`, `LMNA`, `MAP2K1`).
Proteins quantified together as a group are joined with a slash
(e.g. `CALM1/CALM2/CALM3`). Uploaded names must match those in the dataset.

**Configure the inference** (Settings tab):

- *Significance level* α — the threshold below which a corrected p-value is
  called significant (default 0.05).
- *Correction method for drug effects* — default Holm, controlling the
  family-wise error rate.
- *Correction method for protein effects* — default Benjamini–Hochberg,
  controlling the false discovery rate, which suits the much larger number of
  protein–protein tests.

All correction methods available in R's `p.adjust` are offered. Tests are
two-sided against the null of no effect.

**View results** in the *Drug Effects* and *Protein Network* tabs. Each tab
includes a collapsible "How to read..." box explaining the display.

## Data model

The backend exposes ten tables (as Parquet files or SQLite tables). Proteins,
treatments, and time points are referenced throughout by integer index; the
lookup tables map those indices to names.

| Table | Contents |
|---|---|
| `proteins` | `protein` (index), `name` (HGNC symbol) |
| `treatments` | `treatment` (index), `raw_label`, `display` |
| `exp_times` | `idx`, `hours` (6, 24, 48) |
| `drug_order` | `ord`, `raw_label` — display ordering for the heatmap |
| `drug_lookup` | `id`, `name` — resolves drug IDs appearing in interaction labels |
| `meta` | `key`, `value` — `nProtein`, `nTreatment`, `nTimes` |
| `drug_pvalue` | `protein`, `treatment`, `time_idx`, `pvalue` — raw p-values for drug effects |
| `drug_coef` | `protein`, `treatment`, `time_idx`, `coef` — estimated drug effects |
| `protein_pvalue` | `source`, `target`, `time_idx`, `pvalue` — raw p-values for temporal dependencies |
| `protein_coef` | `source`, `target`, `time_idx`, `coef` — estimated dependency coefficients |

Both data access layers return **raw, uncorrected** p-values. All statistical
decisions — multiple-testing correction, α thresholding, minimum-over-time
aggregation, edge signing — happen in the Shiny reactives, not in the data
layer. The `alpha` argument to `dp_net_pvals()` is a performance pre-filter to
avoid returning tens of millions of irrelevant rows, not a statistical choice.

Two storage conventions are worth noting. Sentinel p-values marking drug pairs
without experimental data are dropped at build time rather than stored, so
absence from `protein_pvalue` means "no experiment". Only nonzero entries of
`protein_coef` are stored, so `dp_protein_bhat()` returns zero for any
(target, source) pair not present in the table.

### Data access API

Both `dbAccess_parquet.R` and `dbAccess.R` provide:

```r
dp_connect(...)                                          # open a connection
dp_metadata(con)                                         # names, labels, counts
dp_drug_selPvecs(con, selection, nTreatment, nTimes)     # drug p-values for a protein set
dp_net_pvals(con, selection, alpha, nProtein, nNetTimes) # network p-values for a protein set
dp_drug_dEff(con, treatment_idx, nProtein, nTimes)       # drug coefficients for one treatment
dp_protein_bhat(con, target, time_idx, nProtein)         # dependency coefficients for one target
```

The Parquet backend additionally provides `dp_disconnect(con)` and
`dp_protein_coef_bulk(con, targets, time_idx)`.

## Storage backends

The application ships with two interchangeable backends:

- **Parquet + DuckDB** (`dbAccess_parquet.R`) — the default and the one used in
  deployment. DuckDB reads the Parquet files directly with predicate and column
  pushdown, so a query touches only the row groups and columns it needs.
- **SQLite** (`dbAccess.R`) — a single-file alternative. Network queries rely on
  SQLite's multi-index OR optimisation over `idx_net_source` / `idx_net_target`
  to avoid scanning the full ~58M-row table.

To switch, change which file `app.R` sources near the top:

```r
# source("dbAccess.R")
# db   <- dp_connect()
# meta <- dp_metadata(db)

source("dbAccess_parquet.R")
db <- dp_connect("parquet", in_memory_network = FALSE, memory_limit = "3GB")
meta <- dp_metadata(db)
```

The function names and signatures are identical, so nothing else changes.

## Deployment

The hosted instance runs on shinyapps.io. `rsconnect/` contains the deployment
configuration. To deploy your own instance:

```r
library(rsconnect)
deployApp()
```

The Parquet directory must be included in the bundle. Note that the app's
memory footprint is dominated by the network table; `in_memory_network = FALSE`
is recommended for hosted deployment.

## Reproducing the evidence

This repository contains the **query application** only. The code that fits the
statistical models and assembles the queryable database is in a separate
repository:

<https://github.com/markusul/DrugProt-Paper>

The full pipeline runs as follows:

```
Perturbation proteomics data (Sun et al., 2025)
  │
  ├─ model fitting: de-sparsified Lasso via the hdi package,
  │  producing p-values and coefficients as .RData under results/
  │
  ├─ R/buildDatabase.R      → data/drugprot.sqlite
  │     Collects drug-effect p-values, protein-network p-values, and
  │     drug/protein coefficients into long, indexed SQLite tables.
  │     Sentinel p-values (== 2, "no experiment") are dropped rather
  │     than stored; only nonzero protein coefficients are kept.
  │
  ├─ R/sqliteToParquet.R    → data/parquet/
  │     Writes one compressed Parquet file per SQLite table.
  │
  └─ copy data/parquet/ into this repository as ./parquet/
        (the app calls dp_connect("parquet", ...); see Data files above)
```

The underlying perturbation proteomics dataset is described in Sun et al.
(2025), <https://doi.org/10.1101/2025.02.07.637070>.

## Citation

If you use Drug-Prot, please cite:

> Ulmer, M., Sun, R., Qian, L., Aebersold, R., Guo, T., and Bühlmann, P.
> (2026). Drug-Prot: A query system for statistical inference of drug effects
> and interactions in dynamic proteomic networks.

Please also cite the underlying dataset:

> Sun, R., Qian, L., Li, Y., et al. (2025). A perturbation proteomics-based
> foundation model for virtual cell construction. bioRxiv.
> <https://doi.org/10.1101/2025.02.07.637070>

### Methods references

- Zhang, C.-H. and Zhang, S. S. (2014). Confidence intervals for low
  dimensional parameters in high dimensional linear models. *J. R. Stat. Soc.
  B* 76(1):217–242. (de-sparsified Lasso)
- van de Geer, S., Bühlmann, P., Ritov, Y., and Dezeure, R. (2014). On
  asymptotically optimal confidence regions and tests for high-dimensional
  models. *Ann. Statist.* 42(3):1166–1202.
- Dezeure, R., Bühlmann, P., Meier, L., and Meinshausen, N. (2015).
  High-dimensional inference: confidence intervals, p-values and R-software
  hdi. *Statist. Sci.* 30(4):533–558.
- Bühlmann, P. (2013). Statistical significance in high-dimensional linear
  models. *Bernoulli* 19(4):1212–1242. (group p-values)
- Holm, S. (1979). A simple sequentially rejective multiple test procedure.
  *Scand. J. Statist.* 6(2):65–70. (FWER control for drug effects)
- Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate.
  *J. R. Stat. Soc. B* 57(1):289–300. (FDR control for protein effects)

## License

BSD 3-Clause License. See [`LICENSE`](LICENSE).

## Contact

Markus Ulmer — <markus.ulmer@stat.math.ethz.ch>
Seminar for Statistics, ETH Zürich
