# RAS example — backfat thickness in Duroc pigs

An end-to-end, reproducible RAS **how-to** on a real Duroc pig GWAS dataset
(352 pigs, 36,120 SNPs). It starts from raw PLINK files, previews the data, and
shows how to run RAS on chromosome 12.

> ⚠️ **Results withheld pending publication.** This public notebook shows the
> data handling, the PLINK pre-processing, and how to *call* RAS. The detected
> region, result figures, and biological interpretation are omitted because they
> are part of a manuscript in preparation. Run the notebook yourself on the
> Figshare data to reproduce the full output.

## ▶ Open the notebook: [`RAS_pig_demo.ipynb`](RAS_pig_demo.ipynb)

The notebook displays the data, explains each RAS stage, and runs the analysis
**two ways** — the one-call `ras()` and the step-by-step
`ras_scan → ras_detect → ras_validate`. The result-producing cells are left
unexecuted (their code is shown, outputs hidden); run them yourself to reproduce
the full analysis.

A plain-script version of the same analysis is in
[`run_ras_pig.R`](run_ras_pig.R).

## Setup (to run it yourself)

- **R packages:** `install.packages("readxl")` and
  `remotes::install_github("hepingzhangyale/RAS")`
- **PLINK 1.9** (used for the marginal-GWAS step):
  <https://www.cog-genomics.org/plink/>
- **Data:** download the four files from Figshare
  ([10.6084/m9.figshare.4263317](https://doi.org/10.6084/m9.figshare.4263317))
  into `data/` — see [`data/README.md`](data/README.md).

## References

- **Method:** Y. Jiang & H. Zhang, *Empowering genome-wide association studies via
  a visualizable test based on the regional association score,* PNAS 122(9)
  e2419721122 (2025). https://doi.org/10.1073/pnas.2419721122
- **Dataset (source):** P. G. Eusebi et al., *A genome-wide association analysis
  for carcass traits in a commercial Duroc pig population,* Animal Genetics
  48:466–469 (2017). https://doi.org/10.1111/age.12545
- **SSC12 / ACACA locus (`ALGA0066299`):** S. Kim et al., *QTL fine mapping for
  intramuscular fat and fatty acid composition … on SSC12 …,* Czech J. Anim. Sci.
  64:180–188 (2019). https://doi.org/10.17221/50/2018-CJAS
