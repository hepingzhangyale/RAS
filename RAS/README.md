# RAS: Regional Association Score for GWAS

The RAS package implements the Regional Association Score method for genome-wide association studies. It converts per-SNP effect sizes into a genomic −log₁₀(p) time series and applies changepoint detection to locate peaks that mark significant association regions.

## Installation

```r
# install from GitHub
devtools::install_github("your-github-username/RAS")

# dependency
install.packages("segmented")
```

## Quick start

```r
library(RAS)

# one call does everything
result <- ras(geno, phenotype, covariates,
              covariate_cols = c("age", "sex", paste0("pc", 1:10)),
              is_continuous  = TRUE,
              chrom = 1, save_dir = "results/")

print(result)             # detected changepoint positions
plot(result)              # full-chromosome scan profile
plot(result, zoom = TRUE) # zoomed view around each changepoint
```

## Step-by-step (advanced)

```r
scan     <- ras_scan(geno, phenotype, covariates, ...)
detected <- ras_detect(scan$x, scan$y)
final    <- ras_validate(detected, x = scan$x, y = scan$y, this.skip = 10)
```

## Documentation

```r
?RAS   # package overview and full pipeline details
?ras   # main function
```

## Authors

- Jiahe Jin <jiahe.jin@yale.edu> (maintainer)
- Yiran Jiang <yiran.jiang@uky.edu>
- Heping Zhang <heping.zhang@yale.edu>
