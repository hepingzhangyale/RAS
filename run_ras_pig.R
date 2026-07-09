# =============================================================================
# RAS end-to-end example: backfat thickness in Duroc pigs
# -----------------------------------------------------------------------------
# Data : PLINK binary GWAS genotypes (.bed/.bim/.fam) + phenotype spreadsheet
#        for 352 Duroc pigs, 36,120 SNPs (Illumina PorcineSNP60), 18 autosomes.
#        Source: Figshare 4263317 (Duroc pig genotype and phenotype dataset).
# Goal : locate the genomic region associated with backfat thickness (BFT34R)
#        and show that RAS pinpoints a single regional peak on chromosome 12.
#
# Requirements: install.packages(c("readxl"))  and the RAS package.
#   remotes::install_github("hepingzhangyale/RAS")
# Run:  Rscript run_ras_pig.R      (from this folder)
# =============================================================================

suppressMessages({
  library(readxl)
  library(RAS)
})

set.seed(1)                                   # RAS uses random 50/50 splits
data_dir    <- "data"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. A tiny base-R reader for PLINK 1.0 .bed files (SNP-major)
#    Returns an n x m dosage matrix (count of allele A1): 0, 1, 2, or NA.
# -----------------------------------------------------------------------------
read_plink_bed <- function(bed, n, m) {
  bytes_per_snp <- ceiling(n / 4)
  con <- file(bed, "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  stopifnot(magic[1] == as.raw(0x6c),         # PLINK magic bytes
            magic[2] == as.raw(0x1b),
            magic[3] == as.raw(0x01))          # 0x01 = SNP-major
  raw_all <- readBin(con, "raw", bytes_per_snp * m)
  dim(raw_all) <- c(bytes_per_snp, m)
  # 2-bit codes -> dosage: 00 = 2, 01 = NA, 10 = 1, 11 = 0
  lut  <- c(2L, NA_integer_, 1L, 0L)
  geno <- matrix(NA_integer_, n, m)
  for (j in seq_len(m)) {
    b <- as.integer(raw_all[, j])
    g <- integer(bytes_per_snp * 4)
    g[seq(1, by = 4, length.out = bytes_per_snp)] <- lut[(b %%  4) + 1]
    g[seq(2, by = 4, length.out = bytes_per_snp)] <- lut[((b %/% 4)  %% 4) + 1]
    g[seq(3, by = 4, length.out = bytes_per_snp)] <- lut[((b %/% 16) %% 4) + 1]
    g[seq(4, by = 4, length.out = bytes_per_snp)] <- lut[((b %/% 64) %% 4) + 1]
    geno[, j] <- g[1:n]
  }
  geno
}

# -----------------------------------------------------------------------------
# 2. Load genotypes and SNP map
# -----------------------------------------------------------------------------
fam <- read.table(file.path(data_dir, "GWAS_genotypes.fam"),
                  header = FALSE, fill = TRUE)[, 1:6]
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
bim <- read.table(file.path(data_dir, "GWAS_genotypes.bim"), header = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
n <- nrow(fam); m <- nrow(bim)
cat(sprintf("Loaded map: %d samples x %d SNPs, %d autosomes\n",
            n, m, length(unique(bim$CHR))))

geno <- read_plink_bed(file.path(data_dir, "GWAS_genotypes.bed"), n, m)

# mean-impute the few missing calls (0.13%)
col_mean <- colMeans(geno, na.rm = TRUE)
na_idx   <- which(is.na(geno), arr.ind = TRUE)
if (nrow(na_idx)) geno[na_idx] <- col_mean[na_idx[, 2]]

# -----------------------------------------------------------------------------
# 3. Load phenotypes and align to the genotype sample order
# -----------------------------------------------------------------------------
ph <- as.data.frame(read_excel(
  file.path(data_dir, "gwas_duroc_phenotypes.xlsx"), sheet = 1))[, 1:9]
colnames(ph)[1:2] <- c("FID", "IID")
ph <- ph[match(paste(fam$FID, fam$IID), paste(ph$FID, ph$IID)), ]

trait <- "BFT34R"                              # backfat thickness at 3/4 rib
y     <- suppressWarnings(as.numeric(ph[[trait]]))
cat(sprintf("Trait %s: %d non-missing, range [%.0f, %.0f]\n",
            trait, sum(!is.na(y)), min(y, na.rm = TRUE), max(y, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 4. Genotype principal components as covariates (population structure)
#    All animals are male, so sex is not informative; the top genotype PCs are
#    the standard covariates to control for structure and relatedness.
# -----------------------------------------------------------------------------
af    <- colMeans(geno) / 2
maf   <- pmin(af, 1 - af)
pcs   <- prcomp(scale(geno[, maf > 0.05]), center = FALSE, scale. = FALSE)
covdf <- data.frame(pc1 = pcs$x[, 1], pc2 = pcs$x[, 2], pc3 = pcs$x[, 3])
cat(sprintf("PC covariates: PC1-3 explain %.1f%% of genotypic variance\n",
            100 * sum(pcs$sdev[1:3]^2) / sum(pcs$sdev^2)))

# -----------------------------------------------------------------------------
# 5. (context) marginal single-SNP GWAS with PLINK on chromosome 12
#    Runs PLINK 1.9 (--linear, additive quantitative model). Shows the signal
#    as a *cluster* of SNPs -- the pattern RAS summarises. Requires PLINK:
#    https://www.cog-genomics.org/plink/  (set `plink` to the executable).
# -----------------------------------------------------------------------------
sel    <- which(bim$CHR == 12)
geno12 <- geno[, sel]
bim12  <- bim[sel, ]
cat(sprintf("\nChromosome 12: %d SNPs, %.1f-%.1f Mb\n",
            length(sel), min(bim12$BP) / 1e6, max(bim12$BP) / 1e6))

plink <- "plink"        # PLINK 1.9 executable (on PATH, or give a full path)

# PLINK phenotype file: FID IID BFT34R  (missing coded as -9)
pheno_tab <- data.frame(FID = fam$FID, IID = fam$IID,
                        BFT34R = ifelse(is.na(y), -9, y))
write.table(pheno_tab, file.path(results_dir, "pheno_BFT34R.txt"),
            row.names = FALSE, quote = FALSE)

system2(plink, c("--bfile", file.path(data_dir, "GWAS_genotypes"),
                 "--chr", "12",
                 "--pheno", file.path(results_dir, "pheno_BFT34R.txt"),
                 "--pheno-name", "BFT34R",
                 "--linear", "--allow-no-sex",
                 "--out", file.path(results_dir, "chr12_plink")),
        stdout = FALSE, stderr = FALSE)

gw <- read.table(file.path(results_dir, "chr12_plink.assoc.linear"), header = TRUE)
gw <- gw[gw$TEST == "ADD" & is.finite(gw$P), ]      # additive test rows

bonf <- -log10(0.05 / m)                              # genome-wide threshold
png(file.path(results_dir, "chr12_plink_gwas.png"),
    width = 12 * 150, height = 5 * 150, res = 150)
plot(gw$BP / 1e6, -log10(gw$P), pch = 19, cex = 0.5, col = "#2c7bb6",
     ylim = c(0, bonf * 1.05),
     xlab = "Chr 12 position (Mb)", ylab = expression(-log[10](p)),
     main = "Marginal single-SNP GWAS (PLINK --linear)  -  backfat thickness (BFT34R)")
abline(h = bonf, col = "red", lty = 2)
text(2, bonf, "Bonferroni 0.05", col = "red", pos = 3, cex = 0.9)
dev.off()
top <- gw[which.max(-log10(gw$P)), ]
cat(sprintf("Top marginal SNP (PLINK): %s at %.2f Mb  (-log10p = %.2f)\n",
            top$SNP, top$BP / 1e6, -log10(top$P)))

# -----------------------------------------------------------------------------
# 6. RAS: regional association scan + changepoint detection on chromosome 12
# -----------------------------------------------------------------------------
res <- ras(
  geno           = geno12,
  phenotype      = y,
  covariates     = covdf,
  covariate_cols = c("pc1", "pc2", "pc3"),
  is_continuous  = TRUE,
  num_rep        = 10,          # average over 10 train/holdout splits
  skip1          = 5,           # profile-grid stride (SNP units)
  skip2          = 10,
  chrom          = 12,
  save_dir       = results_dir,
  run_plots      = FALSE
)
print(res)

# -----------------------------------------------------------------------------
# 7. Report and plot
# -----------------------------------------------------------------------------
tau <- res$detection$tau_hats
cat("\nDetected regional changepoint(s):\n")
for (t in tau) {
  if (is.na(t) || t < 1 || t > nrow(bim12)) next
  cat(sprintf("  SNP index %d  ->  %s at %.2f Mb\n",
              t, bim12$SNP[t], bim12$BP[t] / 1e6))
}

plot(res, device = "png")                 # -> results/chr-12-cp-plot.png
plot(res, zoom = TRUE, device = "png")    # -> results/chr-12-zoom.png
cat("\nPlots written to '", results_dir, "/'.\n", sep = "")
