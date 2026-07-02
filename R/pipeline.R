#' Run the RAS Scan (Repetition Screening Loop)
#'
#' Computes the averaged \eqn{-\log_{10}(p)}-value profile over
#' \strong{num_rep} independent 50/50 train/hold-out splits.
#'
#' @param geno Matrix of genotype dosages (\eqn{n} samples \eqn{\times}
#'   \eqn{N} variants). Rows must be aligned with \strong{phenotype} and
#'   \strong{covariates}.
#' @param phenotype Numeric vector of length \eqn{n}. Raw phenotype values.
#'   For continuous traits the function residualises covariates from the
#'   training phenotype internally each repetition; pass the original
#'   un-residualised values here.
#' @param covariates Data frame with \eqn{n} rows aligned with \strong{geno}.
#'   Must contain all columns named in \strong{covariate_cols}; additional
#'   columns (e.g., an ID column) are silently ignored.
#' @param covariate_cols Character vector of column names in \strong{covariates}
#'   to include in the regression models.
#' @param is_continuous Logical. \code{TRUE} for quantitative traits,
#'   \code{FALSE} for binary (case/control) traits.
#' @param num_rep Integer. Number of independent 50/50 splits to average over.
#'   Default \code{5}.
#' @param skip1 Integer. Step size for the primary SNP position grid
#'   (\code{seq(1, N, by = skip1)}).  Default \code{10}.
#' @param skip2 Integer. Sub-step used to grow the scan window at each grid
#'   position.  Default \code{20}.
#' @param chrom Integer. Chromosome number used in saved file names.
#'   Default \code{1}.
#' @param save_dir Character. Directory for intermediate \code{.rds} output
#'   files; created recursively if it does not exist.  Default
#'   \code{"./result"}.
#' @param min_window_size Integer. Minimum scan window half-size passed to
#'   \code{\link{screen_forward_max_region}}.  Default \code{5}.
#' @param max_window_size Integer. Maximum scan window half-size passed to
#'   \code{\link{screen_forward_max_region}}.  Default \code{100}.
#' @param scan_test Character. Per-window test for binary traits, passed to
#'   \code{\link{screen_forward_max_region}}.  \code{"glm"} (default) fits a
#'   full logistic model per window (Wald p-value); \code{"score"} fits the
#'   null logistic model once and uses a Rao score test per window
#'   (asymptotically equivalent, much faster at scale).  Has no effect for
#'   continuous traits.
#'
#' @details
#' For each of the \strong{num_rep} repetitions the function performs three
#' steps:
#'
#' \enumerate{
#'   \item \strong{GWAS weights.}  A fresh 50/50 random split of all
#'     \eqn{n} samples is drawn.  \code{\link{compute_gwas_weights}} fits a
#'     per-SNP regression on the training half and returns an
#'     \eqn{N \times 4} coefficient matrix; the first column (effect-size
#'     estimates) is used as PGS weights.
#'   \item \strong{PGS contribution matrix.}
#'     \code{\link{compute_pgs_matrix}} multiplies each hold-out
#'     individual's dosage vector by the per-SNP weights, producing an
#'     \eqn{n_{\text{holdout}} \times N} matrix of weighted contributions.
#'   \item \strong{Forward scan.}
#'     \code{\link{screen_forward_max_region}} slides an expanding window
#'     across the genome.  At each position it accumulates weighted dosages,
#'     regresses them against the hold-out phenotype, and records the minimum
#'     \eqn{p}-value over window sizes in
#'     \code{[min_window_size, max_window_size]}.  The result is a vector of
#'     \eqn{-\log_{10}(p)} values on the grid \code{seq(1, N, by = skip1)}.
#' }
#'
#' The \eqn{-\log_{10}(p)} vectors from all repetitions are summed and
#' divided by \strong{num_rep} to form the final profile.
#'
#' For \strong{continuous traits} (\code{is_continuous = TRUE}), covariates
#' are residualised from the training phenotype using only training
#' individuals before the GWAS step, so no hold-out information leaks into
#' the effect-size estimates.  The residualisation is repeated from the
#' original \code{phenotype} vector each repetition.  For \strong{binary
#' traits} (\code{is_continuous = FALSE}) the raw phenotype is passed
#' directly to the GWAS step, and the scan step fits a logistic model via
#' \code{glm(family = binomial())}.
#'
#' After each repetition the large intermediate matrices (\code{coef.mat},
#' \code{pgs.mat}) are removed from the R session and
#' \code{\link{release_memory}} is called to return free heap pages to the
#' OS.  On Linux/glibc this calls \code{malloc_trim(0)} and can
#' substantially reduce RSS between repetitions.
#'
#' @return
#' Invisibly returns a named list with two elements:
#' \describe{
#'   \item{\code{x}}{Integer vector of length \eqn{\lceil N /
#'     \texttt{skip1} \rceil}. The SNP position index grid
#'     \code{seq(1, N, by = skip1)}.}
#'   \item{\code{y}}{Numeric vector, same length as \code{x}. Averaged
#'     \eqn{-\log_{10}(p)}-value profile across all \strong{num_rep}
#'     repetitions.}
#' }
#' The following \code{.rds} files are written to \code{save_dir}:
#' \describe{
#'   \item{\code{chr-<chrom>_coef_mat-<rep>.rds}}{The \eqn{N \times 4}
#'     GWAS coefficient matrix from repetition \code{rep}.}
#'   \item{\code{mean_p_values_chr<chrom>_reps1-<num_rep>.rds}}{The
#'     averaged \eqn{-\log_{10}(p)} vector \code{y}.}
#' }
#'
#' @seealso
#' \code{\link{ras}} for a single-call wrapper that also runs
#' changepoint detection and plotting.
#' \code{\link{compute_gwas_weights}}, \code{\link{compute_pgs_matrix}},
#' \code{\link{screen_forward_max_region}} for the constituent steps.
#' \code{\link{ras_detect}},
#' \code{\link{ras_validate}} for downstream changepoint analysis.
#' \code{\link{ras_memory}} to check memory requirements before
#' loading large genotype data.
#' \code{\link{release_memory}} for OS-level heap reclamation.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n_samp <- 80; n_snp <- 60
#' geno <- matrix(
#'   sample(0:2, n_samp * n_snp, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
#'   nrow = n_samp, ncol = n_snp)
#' pheno <- rnorm(n_samp)
#' cov_df <- data.frame(age = rnorm(n_samp), sex = rbinom(n_samp, 1, 0.5))
#' cov_df$age_squared <- cov_df$age^2
#' cov_df$age_sex     <- cov_df$age * cov_df$sex
#' for (i in 1:10) cov_df[[paste0("pc", i)]] <- rnorm(n_samp)
#'
#' scan <- run_ras_scan(
#'   geno            = geno,
#'   phenotype       = pheno,
#'   covariates      = cov_df,
#'   covariate_cols  = c("age", "sex", "age_squared", "age_sex",
#'                       paste0("pc", 1:10)),
#'   is_continuous   = TRUE,
#'   num_rep         = 2,
#'   skip1           = 5,
#'   skip2           = 5,
#'   min_window_size = 2,
#'   max_window_size = 10,
#'   chrom           = 1,
#'   save_dir        = tempdir()
#' )
#' plot(scan$x, scan$y, type = "l",
#'      xlab = "SNP index", ylab = expression(-log[10](p)),
#'      main = "RAS scan profile")
#' }
#' @export
ras_scan <- function(geno, phenotype, covariates, covariate_cols,
                         is_continuous,
                         num_rep        = 5,
                         skip1          = 10,
                         skip2          = 20,
                         chrom          = 1,
                         save_dir       = "./result",
                         min_window_size = 5,
                         max_window_size = 100,
                         scan_test       = c("glm", "score")) {

  scan_test <- match.arg(scan_test)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  n <- nrow(geno)
  N <- ncol(geno)

  # Restrict covariate df to the required columns only (e.g. strips any ID col)
  cov.df <- covariates[, covariate_cols, drop = FALSE]

  # Two separate formula strings — the GWAS model includes the SNP term
  # ("this.x"), while the scan model uses the PGS term ("this.pgs") instead.
  gwas_formula <- paste("this.x +", paste(covariate_cols, collapse = " + "))
  scan_formula  <- paste(covariate_cols, collapse = " + ")

  this.seq      <- seq(1, N, by = skip1)
  full.p.values <- rep(0, length(this.seq))

  for (this.rep in seq_len(num_rep)) {
    cat("==============================\n")
    cat(" Repetition", this.rep, "/", num_rep, "\n")

    # 50/50 random split
    this.sample  <- sort(sample(n, n %/% 2))
    this.leftout <- setdiff(seq_len(n), this.sample)

    # ── Phenotype preparation ────────────────────────────────────────────────
    # Continuous: regress covariates out of the TRAINING split to get residuals.
    # Residualisation is done fresh each repetition from the original phenotype,
    # using only training individuals (no holdout data leakage).
    if (is_continuous) {
      train.df            <- cov.df[this.sample, , drop = FALSE]
      train.df$phenotype1 <- phenotype[this.sample]
      lm0        <- lm(
        as.formula(paste("phenotype1 ~", paste(covariate_cols, collapse = " + "))),
        data = train.df
      )
      phenotype1 <- as.numeric(lm0$residuals)
    } else {
      phenotype1 <- phenotype[this.sample]
    }
    phenotype2 <- phenotype[this.leftout]

    # ── Step a: per-SNP GWAS on training split ───────────────────────────────
    coef.mat <- compute_gwas_weights(
      geno              = geno,
      phenotype1        = phenotype1,
      this.sample       = this.sample,
      this.df           = cov.df,
      is_continuous     = is_continuous,
      covariate_formula = gwas_formula
    )
    pgs.weights <- coef.mat[, 1]
    saveRDS(coef.mat, file.path(save_dir,
      paste0("chr-", chrom, "_coef_mat-", this.rep, ".rds")))

    # ── Step b: PGS contribution matrix for holdout individuals ──────────────
    pgs.mat <- compute_pgs_matrix(
      geno         = geno,
      this.leftout = this.leftout,
      pgs.weights  = pgs.weights
    )

    # ── Step c: forward scan on holdout covariates + holdout phenotype ───────
    scan.df            <- cov.df[this.leftout, , drop = FALSE]
    scan.df$phenotype2 <- phenotype2

    return.p.values <- screen_forward_max_region(
      geno              = geno,
      pgs.mat           = pgs.mat,
      this.df           = scan.df,
      num_signals       = -1,
      start.point       = 1,
      save.directory    = save_dir,
      this.chrome       = chrom,
      min_window_size   = min_window_size,
      max_window_size   = max_window_size,
      isSimulation      = FALSE,
      this.repetition   = this.rep,
      screening_round   = 1,
      isPlot            = FALSE,
      skip1             = skip1,
      skip2             = skip2,
      is_continuous     = is_continuous,
      covariate_formula = scan_formula,
      scan_test         = scan_test
    )

    full.p.values <- full.p.values + return.p.values

    rm(coef.mat, pgs.weights, pgs.mat, return.p.values, scan.df)
    release_memory(verbose = TRUE)
  }

  mean.p.values <- full.p.values / num_rep
  saveRDS(mean.p.values,
    file.path(save_dir,
      paste0("mean_p_values_chr", chrom, "_reps1-", num_rep, ".rds")))
  cat("Scanning complete!\n")

  invisible(list(x = this.seq, y = mean.p.values))
}


#' RAS: Regional Association Score Analysis
#'
#' Detect significant association regions in GWAS data.
#' Returns an object of class \code{"ras"}.
#'
#' @section Simple (recommended):
#' \preformatted{
#' ## one call does everything
#' result <- ras(geno, phenotype, covariates,
#'               covariate_cols = c("age", "sex", paste0("pc", 1:10)),
#'               is_continuous  = TRUE,
#'               chrom = 1, save_dir = "results/")
#'
#' print(result)             # detected changepoint positions
#' plot(result)              # full-chromosome scan profile
#' plot(result, zoom = TRUE) # zoomed view around each changepoint
#' }
#'
#' For step-by-step control see \code{\link{ras_scan}},
#' \code{\link{ras_detect}}, \code{\link{ras_validate}}.
#'
#' @inheritParams ras_scan
#'
#' @param cp_p_threshold Numeric. Davies test p-value threshold for
#'   first-pass candidate changepoints.  Default \code{0.01}.
#' @param cp_window_size Integer. Sliding window width for first-pass
#'   detection.  Default \code{3000}.
#' @param cp_min_length Integer. Minimum segment length before and after a
#'   candidate changepoint.  Default \code{10}.
#' @param cp_slope_check_window Integer. Half-width of the local window used
#'   to verify slope direction around each candidate.  Default \code{30}.
#' @param cp_slope_left Numeric. One-tailed p-value threshold for the
#'   left-side slope test (tests that the slope to the left of the candidate
#'   is significantly positive).  Default \code{1e-10}.
#' @param cp_slope_right Numeric. One-tailed p-value threshold for the
#'   right-side slope test (tests that the slope to the right is
#'   significantly negative).  Default \code{1e-20}.
#' @param second_window_size Integer. Half-window size for second-pass local
#'   Davies tests.  Default \code{50}.
#' @param second_p_threshold Numeric. Davies test p-value threshold for
#'   second-pass validation.  Default \code{1e-10}.
#' @param run_plots Logical. If \code{TRUE} (default), saves diagnostic plots
#'   via \code{\link{plot_ras_scan}} and \code{\link{plot_ras_zoom_regions}}.
#' @param plot_device Character. Output device for plots: \code{"pdf"}
#'   (default), \code{"png"}, or \code{"screen"}.
#' @param plot_p_threshold Numeric. Significance reference line drawn on plots
#'   (on the \eqn{-\log_{10}} scale).  Default \code{8}.
#' @param plot_y_cap Numeric or \code{NULL}. If provided, the y-axis is
#'   capped at this value and positions above the cap are annotated with
#'   arrows.  Default \code{NULL} (no cap).
#' @param min_signal Numeric. Minimum \eqn{-\log_{10}(p)} scan value required
#'   to retain a validated changepoint (passed to \code{\link{ras_validate}});
#'   also used as the low-signal colour boundary in all diagnostic plots.
#'   Default \code{2.5}.
#'
#' @details
#' For a stage-by-stage description of the pipeline see \code{\link{RAS}}.
#'
#' @return
#' Invisibly returns a named list with two elements:
#' \describe{
#'   \item{\code{scan}}{A list with elements:
#'     \describe{
#'       \item{\code{x}}{Integer vector. SNP position index grid
#'         \code{seq(1, N, by = skip1)}.}
#'       \item{\code{y}}{Numeric vector (same length as \code{x}). Averaged
#'         \eqn{-\log_{10}(p)}-value profile across all repetitions.}
#'     }
#'   }
#'   \item{\code{detection}}{A list returned by \code{\link{ras_validate}}
#'     with elements:
#'     \describe{
#'       \item{\code{tau_hats}}{Integer vector. Validated changepoint positions
#'         (re-mapped to genomic coordinates via \strong{this.start} and
#'         \strong{this.skip}).}
#'       \item{\code{all.changepoints}}{Integer vector. All candidate positions
#'         examined during the first pass, re-mapped to genomic coordinates.}
#'       \item{\code{all.p.values}}{Numeric vector. \eqn{-\log_{10}(p)} Davies
#'         test values for every candidate in \code{all.changepoints}.}
#'       \item{\code{left.slopes}}{Numeric vector. Estimated left-side slopes
#'         at each validated changepoint.}
#'       \item{\code{right.slopes}}{Numeric vector. Estimated right-side slopes
#'         at each validated changepoint.}
#'     }
#'   }
#' }
#' Plot files (if \code{run_plots = TRUE}) are written to \code{save_dir}
#' with names \code{chr-<chrom>-cp-plot.<ext>},
#' \code{chr-<chrom>-cp-p-values-plot.<ext>}, and
#' \code{chr-<chrom>-zoom.<ext>}.
#'
#' @seealso
#' \code{\link{ras_scan}} for the scan-only step (useful when tuning
#' changepoint parameters separately).
#' \code{\link{ras_detect}},
#' \code{\link{ras_validate}} for the changepoint detection steps.
#' \code{\link{plot.ras}} for the plotting step.
#' \code{\link{ras_memory}} to check memory requirements before
#' loading large genotype data.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n_samp <- 80; n_snp <- 60
#' geno <- matrix(
#'   sample(0:2, n_samp * n_snp, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
#'   nrow = n_samp, ncol = n_snp)
#' pheno <- rnorm(n_samp)
#' cov_df <- data.frame(age = rnorm(n_samp), sex = rbinom(n_samp, 1, 0.5))
#' cov_df$age_squared <- cov_df$age^2
#' cov_df$age_sex     <- cov_df$age * cov_df$sex
#' for (i in 1:10) cov_df[[paste0("pc", i)]] <- rnorm(n_samp)
#'
#' result <- ras(
#'   geno              = geno,
#'   phenotype         = pheno,
#'   covariates        = cov_df,
#'   covariate_cols    = c("age", "sex", "age_squared", "age_sex",
#'                         paste0("pc", 1:10)),
#'   is_continuous     = TRUE,
#'   num_rep           = 2,
#'   skip1             = 5,
#'   skip2             = 5,
#'   min_window_size   = 2,
#'   max_window_size   = 10,
#'   chrom             = 1,
#'   save_dir          = tempdir(),
#'   run_plots         = TRUE,
#'   plot_device       = "png"
#' )
#'
#' cat("Detected changepoints:", result$detection$tau_hats, "\n")
#' plot(result$scan$x, result$scan$y, type = "l",
#'      xlab = "SNP index", ylab = expression(-log[10](p)),
#'      main = "RAS scan profile")
#' }
#' @export
ras <- function(geno, phenotype, covariates, covariate_cols,
                             is_continuous,
                             num_rep               = 5,
                             skip1                 = 10,
                             skip2                 = 20,
                             chrom                 = 1,
                             save_dir              = "./result",
                             min_window_size       = 5,
                             max_window_size       = 100,
                             scan_test             = c("glm", "score"),
                             cp_p_threshold        = 0.01,
                             cp_window_size        = 3000,
                             cp_min_length         = 10,
                             cp_slope_check_window = 30,
                             cp_slope_left         = 1e-10,
                             cp_slope_right        = 1e-20,
                             second_window_size    = 50,
                             second_p_threshold    = 1e-10,
                             min_signal            = 2.5,
                             run_plots             = TRUE,
                             plot_device           = "pdf",
                             plot_p_threshold      = 8,
                             plot_y_cap            = NULL) {

  scan_test <- match.arg(scan_test)

  # ── Stage 1: scan ─────────────────────────────────────────────────────────
  scan <- ras_scan(
    geno            = geno,
    phenotype       = phenotype,
    covariates      = covariates,
    covariate_cols  = covariate_cols,
    is_continuous   = is_continuous,
    num_rep         = num_rep,
    skip1           = skip1,
    skip2           = skip2,
    chrom           = chrom,
    save_dir        = save_dir,
    min_window_size = min_window_size,
    max_window_size = max_window_size,
    scan_test       = scan_test
  )

  x <- scan$x
  y <- scan$y

  # ── Stage 2: first-pass changepoint detection ──────────────────────────────
  cat("==============================\n")
  cat(" Changepoint Detection (Pass 1)\n")
  cp_result <- ras_detect(
    x                              = x,
    y                              = y,
    p.values.threshold             = cp_p_threshold,
    min.length                     = cp_min_length,
    window_size                    = cp_window_size,
    slope_check_window_size        = cp_slope_check_window,
    slope.p.values.threshold.left  = cp_slope_left,
    slope.p.values.threshold.right = cp_slope_right
  )

  # ── Stage 3: second-pass validation ───────────────────────────────────────
  cat("==============================\n")
  cat(" Changepoint Validation (Pass 2)\n")
  detection <- ras_validate(
    this.result        = cp_result,
    x                  = x,
    y                  = y,
    this.start         = 1,
    this.skip          = skip1,
    second_window_size = second_window_size,
    p.value.threshold  = second_p_threshold,
    min_signal         = min_signal
  )

  n_detected <- length(detection$tau_hats)
  cat("==============================\n")
  if (n_detected == 0) {
    cat(" No changepoints detected.\n")
  } else {
    cat(sprintf(" Detected %d changepoint(s) at position(s): %s\n",
                n_detected, paste(detection$tau_hats, collapse = ", ")))
  }

  result <- structure(
    list(scan = scan, detection = detection, chrom = chrom, save_dir = save_dir),
    class = "ras"
  )

  # ── Stage 4: plots ────────────────────────────────────────────────────────
  if (run_plots) {
    cat(" Generating plots...\n")
    plot.ras(result, zoom = FALSE, device = plot_device,
             p.threshold = plot_p_threshold, y_cap = plot_y_cap,
             min_signal  = min_signal)
    plot.ras(result, zoom = TRUE,  device = plot_device,
             p.threshold = plot_p_threshold, min_signal = min_signal)
    cat(sprintf(" Plots saved to: %s\n", save_dir))
  }

  cat("==============================\n")
  cat(" Pipeline complete.\n")

  invisible(result)
}
