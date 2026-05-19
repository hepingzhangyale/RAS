# RAS: Regional Association Score for GWAS

The **RAS** package implements the Regional Association Score method for genome-wide association studies. It converts per-SNP effect sizes into a genomic −log₁₀(*p*) time series and applies changepoint detection to locate peaks that mark significant association regions. The method supports both continuous and binary traits.

If you use this package in your research, please cite:

> Y. Jiang & H. Zhang, Empowering genome-wide association studies via a visualizable test based on the regional association score, *Proc. Natl. Acad. Sci. U.S.A.* 122(9) e2419721122 (2025). https://doi.org/10.1073/pnas.2419721122

---

## Installation

### Method 1 — From GitHub via `remotes` (recommended)

`remotes` is a lightweight alternative to `devtools` with no extra dependencies.

```r
# install remotes if you don't have it
install.packages("remotes")

remotes::install_github("hepingzhangyale/RAS")
```

### Method 2 — From GitHub via `devtools`

```r
# install devtools if you don't have it
install.packages("devtools")

devtools::install_github("hepingzhangyale/RAS")
```

### Method 3 — Clone and install from source

Use this when you need a specific branch, want to inspect the code before installing,
or are working on a machine without direct GitHub access.

```bash
# 1. clone the repository
git clone https://github.com/hepingzhangyale/RAS.git

# 2. install from local source (run inside R)
```

```r
devtools::install("path/to/RAS")   # replace with your local clone path
# e.g. devtools::install("C:/Users/you/RAS")
```

> **Windows note:** source installation requires
> [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to be installed.
> Download the version matching your R installation and make sure it is on the PATH.

### Method 4 — Windows binary zip (no compiler required)

Download the pre-built `.zip` from the
[Releases](https://github.com/hepingzhangyale/RAS/releases) page, then install locally:

```r
install.packages(
  "RAS_0.1.8.zip",          # path to the downloaded zip
  repos = NULL,
  type  = "win.binary"
)
```

---

## Dependencies

RAS imports two packages that are installed automatically:

| Package | Role |
|---------|------|
| `segmented` | Segmented regression and Davies test for changepoint detection |
| `parallel` | CPU core detection used by `ras_memory()` diagnostics |

If automatic installation fails, install them manually first:

```r
install.packages(c("segmented", "parallel"))
```

---

## Quick start

```r
library(RAS)

result <- ras(
  geno, phenotype, covariates,
  covariate_cols = c("age", "sex", paste0("pc", 1:10)),
  is_continuous  = TRUE,
  chrom          = 1,
  save_dir       = "results/"
)

print(result)              # detected changepoint positions
plot(result)               # full-chromosome scan profile
plot(result, zoom = TRUE)  # zoomed view around each changepoint
```

## Step-by-step (advanced)

```r
# Step 1: compute averaged -log10(p) profile
scan <- ras_scan(
  geno, phenotype, covariates,
  covariate_cols = c("age", "sex", paste0("pc", 1:10)),
  is_continuous  = TRUE,
  chrom = 1, save_dir = "results/"
)

# Step 2: first-pass changepoint detection
detected <- ras_detect(
  scan$x, scan$y,
  window_size                    = 3000,
  slope.p.values.threshold.left  = 1e-10,
  slope.p.values.threshold.right = 1e-20
)

# Step 3: second-pass validation
final <- ras_validate(
  detected, x = scan$x, y = scan$y,
  this.skip         = 10,
  p.value.threshold = 1e-10
)

# Step 4: plot
result <- structure(
  list(scan = scan, detection = final, chrom = 1, save_dir = "results/"),
  class = "ras"
)
plot(result)
```

## Documentation

```r
?RAS          # package overview and full pipeline description
?ras          # main one-call entry point
?ras_scan     # Stage 1: scan
?ras_detect   # Stage 2: first-pass changepoint detection
?ras_validate # Stage 3: second-pass validation
?plot.ras     # plotting
?ras_memory   # memory and CPU diagnostics
```

---

## Citation

If you use RAS in your research, please cite:

Y. Jiang & H. Zhang, Empowering genome-wide association studies via a visualizable test based on the regional association score, *Proc. Natl. Acad. Sci. U.S.A.* 122(9) e2419721122 (2025). https://doi.org/10.1073/pnas.2419721122

## Authors

- Jiahe Jin &lt;jiahe.jin@yale.edu&gt; (maintainer)
- Yiran Jiang &lt;yiran.jiang@uky.edu&gt;
- Heping Zhang &lt;heping.zhang@yale.edu&gt;
