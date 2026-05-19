#' Compute Per-SNP GWAS Effect Size Weights
#'
#' Runs a per-SNP regression on the training split to obtain effect-size
#' estimates used as polygenic score (PGS) weights in the forward scan.
#'
#' @param geno Matrix of genotype dosages (\eqn{n} samples \eqn{\times}
#'   \eqn{N} variants). Only the rows indexed by \code{this.sample} are used.
#' @param phenotype1 Numeric vector of length \eqn{|\texttt{this.sample}|}.
#'   Training-split phenotype values.  For continuous traits these should be
#'   covariate-residualised values; for binary traits, the raw 0/1 indicator.
#' @param this.sample Integer vector. Row indices of the training split in
#'   \code{geno}.
#' @param this.df Data frame of covariates with rows aligned to \code{geno}.
#'   Required for binary traits; used to subset and build the regression data
#'   frame.  Ignored for continuous traits.
#' @param is_continuous Logical. \code{TRUE} for quantitative traits,
#'   \code{FALSE} for binary (case/control) traits.
#' @param covariate_formula Character. The right-hand side of the regression
#'   formula including the SNP term \code{this.x}.  For example,
#'   \code{"this.x + age + sex + pc1"}.  If \code{NULL}, a default formula
#'   with sex, age, age-squared, age-sex interaction, and the first ten
#'   ancestry PCs is used.
#'
#' @details
#' For \strong{continuous traits} (\code{is_continuous = TRUE}) the function
#' fits \code{lm(phenotype1 ~ this.x)} for each SNP, excluding samples with
#' missing dosage via \code{na.omit}.  The covariates are assumed to have
#' already been residualised from \code{phenotype1} before this function is
#' called (as done in \code{\link{run_ras_scan}}).
#'
#' For \strong{binary traits} (\code{is_continuous = FALSE}) the function
#' builds a data frame combining the training-split rows of \code{this.df}
#' with the dosage column (\code{this.x}) and phenotype (\code{phenotype1}),
#' removes incomplete cases with \code{na.omit}, and fits
#' \code{lm(phenotype1 ~ covariate_formula)}.  Note that \code{lm} rather
#' than \code{glm(family = binomial)} is used at the GWAS weight stage
#' to obtain a linear probability approximation to the effect size; this is
#' consistent with standard GWAS practice for weight derivation.
#'
#' A progress message is printed every 10,000 SNPs.
#'
#' @return
#' A numeric matrix with \eqn{N} rows and 4 columns:
#' \describe{
#'   \item{\code{Estimate}}{Per-SNP effect-size estimate (slope of the dosage
#'     term).  The first column is used as PGS weights by
#'     \code{\link{compute_pgs_matrix}}.}
#'   \item{\code{Std. Error}}{Standard error of the estimate.}
#'   \item{\code{t value}}{t-statistic.}
#'   \item{\code{Pr(>|t|)}}{Two-sided p-value.}
#' }
#'
#' @seealso
#' \code{\link{compute_pgs_matrix}} for the next step, which multiplies
#' these weights by the hold-out genotype matrix.
#' \code{\link{ras_scan}} for the recommended high-level entry point that
#' calls this function automatically.
#'
#' @examples
#' set.seed(1)
#' n_samp <- 60; n_snp <- 10
#' geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE,
#'                       prob = c(0.6, 0.3, 0.1)),
#'                nrow = n_samp, ncol = n_snp)
#' pheno_train <- rnorm(n_samp)
#'
#' coef_mat <- compute_gwas_weights(
#'   geno              = geno,
#'   phenotype1        = pheno_train,
#'   this.sample       = seq_len(n_samp),
#'   this.df           = data.frame(row.names = seq_len(n_samp)),
#'   is_continuous     = TRUE
#' )
#' dim(coef_mat)        # n_snp x 4
#' head(coef_mat, 3)
#' @export
compute_gwas_weights <- function(geno, phenotype1, this.sample, this.df,
                                 is_continuous,
                                 covariate_formula = NULL) {
  if (is.null(covariate_formula)) {
    covariate_formula <- "this.x + sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
  }

  N <- ncol(geno)
  coef.mat <- matrix(NA, nrow = N, ncol = 4)

  cat("Starting GWAS ... \n")

  for (i in 1:N) {
    if (i %% 10000 == 0) { print(i) }

    this.x <- geno[this.sample, i]

    if (is_continuous) {
      fit <- lm(phenotype1[is.na(this.x) == FALSE] ~ na.omit(this.x))
    } else {
      this.df.sub          <- this.df[this.sample, ]
      this.df.sub$this.x   <- this.x
      this.df.sub$phenotype1 <- phenotype1
      clean_df             <- na.omit(this.df.sub)
      fit_formula          <- as.formula(paste("phenotype1 ~", covariate_formula))
      fit                  <- lm(fit_formula, data = clean_df)
    }

    coef.mat[i, ] <- summary(fit)$coefficients[2, ]
  }

  colnames(coef.mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  return(coef.mat)
}


#' Build the Per-Individual PGS Contribution Matrix
#'
#' Multiplies each hold-out individual's genotype dosages by the per-SNP PGS
#' weights to produce a matrix of weighted dosage contributions used in the
#' forward scan.
#'
#' @param geno Matrix of genotype dosages (\eqn{n} samples \eqn{\times}
#'   \eqn{N} variants).  \code{NA} values are replaced with \code{0}
#'   in-place before multiplication.
#' @param this.leftout Integer vector. Row indices of the hold-out split in
#'   \code{geno}.
#' @param pgs.weights Numeric vector of length \eqn{N}. Per-SNP effect-size
#'   weights, typically the \code{Estimate} column from
#'   \code{\link{compute_gwas_weights}}.
#'
#' @details
#' Missing genotype values (\code{NA}) in \code{geno} are imputed to zero
#' before computing the product.  This is a simple mean-imputation equivalent
#' under a centred dosage scale and is sufficient for the forward scan, where
#' the primary goal is to aggregate regional signals rather than obtain
#' individual-level accuracy.
#'
#' The function modifies \code{geno} in the calling frame due to R's
#' copy-on-modify semantics: a copy of \code{geno} is made when the
#' \code{NA}-replacement assignment is executed.  Callers working with very
#' large genotype matrices should be aware that this temporarily doubles the
#' memory footprint of \code{geno}; use \code{\link{ras_memory}} to
#' check feasibility before running.
#'
#' @return
#' A numeric matrix of dimensions
#' \eqn{|\texttt{this.leftout}| \times N}, where entry
#' \eqn{[i, j]} is the weighted dosage contribution of SNP \eqn{j} for
#' hold-out individual \eqn{i} (i.e., \code{geno[this.leftout[i], j] *
#' pgs.weights[j]}).
#'
#' @seealso
#' \code{\link{compute_gwas_weights}} for the step that produces
#' \code{pgs.weights}.
#' \code{\link{screen_forward_max_region}} for the step that consumes this
#' matrix.
#' \code{\link{ras_memory}} for pre-flight memory estimation.
#'
#' @examples
#' set.seed(2)
#' geno    <- matrix(sample(0:2, 50 * 20, replace = TRUE), nrow = 50, ncol = 20)
#' weights <- rnorm(20)
#' leftout <- 1:15
#'
#' pgs_mat <- compute_pgs_matrix(geno, leftout, weights)
#' dim(pgs_mat)   # 15 x 20
#'
#' ## Entry [i, j] equals geno[leftout[i], j] * weights[j]
#' stopifnot(pgs_mat[1, 1] == geno[leftout[1], 1] * weights[1])
#' @export
compute_pgs_matrix <- function(geno, this.leftout, pgs.weights) {
  geno[is.na(geno)] <- 0
  pgs.mat <- matrix(data = NA, nrow = length(this.leftout), ncol = ncol(geno))
  for (i in 1:length(this.leftout)) {
    pgs.mat[i, ] <- geno[this.leftout[i], ] * pgs.weights
  }
  return(pgs.mat)
}
