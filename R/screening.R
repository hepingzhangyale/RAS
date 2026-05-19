#' Forward Scan to Compute the RAS Profile
#'
#' Slides an expanding window across the genome, accumulating weighted
#' genotype contributions and testing association with the hold-out phenotype
#' at each position.
#'
#' @param geno Matrix of genotype dosages (\eqn{n} samples \eqn{\times}
#'   \eqn{N} variants).
#' @param pgs.mat Matrix of pre-computed per-variant PGS contributions
#'   (\eqn{n_{\text{holdout}} \times N}), as returned by
#'   \code{\link{compute_pgs_matrix}}.
#' @param this.df Data frame containing the hold-out phenotype in a column
#'   named \code{phenotype2} and any covariate columns named in
#'   \code{covariate_formula}.  Rows must correspond to hold-out individuals.
#' @param num_signals Integer. Number of true signals; used only for
#'   simulation plots (\code{isSimulation = TRUE}) to draw reference lines.
#'   Pass \code{-1} for real data.
#' @param start.point Integer. Starting SNP index for the scan.
#'   Default \code{1}.
#' @param save.directory Character. Directory for optional PDF plots.
#'   Default \code{"simulation"}.
#' @param this.chrome Integer. Chromosome number for non-simulation plot
#'   filenames.  Default \code{1}.
#' @param min_window_size Integer. Minimum scan window half-size.
#'   Default \code{5}.
#' @param max_window_size Integer. Maximum scan window half-size.
#'   Default \code{100}.
#' @param isSimulation Logical. If \code{TRUE}, use simulation-mode plot
#'   filenames and draw true-signal reference lines.  Default \code{TRUE}.
#' @param this.repetition Integer. Repetition index used in plot filenames.
#'   Default \code{1}.
#' @param screening_round Integer. Number of forward screening rounds.
#'   Default \code{1}.
#' @param isPlot Logical. If \code{TRUE}, save a PDF of the scan profile.
#'   Default \code{TRUE}.
#' @param skip1 Integer. Step size for the primary SNP position grid.
#'   Default \code{100}.
#' @param skip2 Integer. Step size for growing the window within each grid
#'   position.  Default \code{5}.
#' @param is_continuous Logical. If \code{TRUE}, fits \code{lm}; otherwise
#'   \code{glm(family = binomial())}.
#' @param covariate_formula Character. Right-hand side covariates (without
#'   the PGS term \code{this.pgs}) in the scan regression formula.  If
#'   \code{NULL}, a default formula with sex, age, age-squared, age-sex
#'   interaction, and the first ten ancestry PCs is used.
#' @param signal.starts Integer vector. True signal start positions for
#'   simulation reference lines.  Default \code{NULL}.
#' @param signal.window.size Integer. Width of each true signal window for
#'   simulation reference lines.  Default \code{NULL}.
#'
#' @details
#' At each grid position \eqn{j \in \{1, 1 + \texttt{skip1}, \ldots, N\}}
#' the function builds a polygenic score (PGS) by accumulating columns of
#' \code{pgs.mat} from a growing symmetric window
#' \eqn{[j - w + 1,\; j + w - 1]} for
#' \eqn{w \in \{0, \texttt{min\_window\_size}, \ldots, \texttt{max\_window\_size}\}}.
#' The accumulation is incremental: only the newly added columns are summed
#' at each step, avoiding redundant computation.
#'
#' At window size \eqn{w = 0}, only column \eqn{j} of \code{pgs.mat} is
#' used (the single-SNP PGS).  For each window size the PGS is tested
#' against \code{phenotype2} via \code{lm} (continuous) or
#' \code{glm(family = binomial())} (binary).  The minimum p-value over all
#' window sizes is stored for position \eqn{j}.
#'
#' The final return value is \eqn{-\log_{10}} of these per-position minimum
#' p-values.
#'
#' @return
#' Numeric vector of length \eqn{\lceil N / \texttt{skip1} \rceil}.
#' Each element is the \eqn{-\log_{10}(p)}-value at the corresponding grid
#' position, where the p-value is the minimum over all tested window sizes.
#'
#' @seealso
#' \code{\link{compute_pgs_matrix}} for the step that produces
#' \code{pgs.mat}.
#' \code{\link{ras_scan}} for the recommended high-level entry point that
#' calls this function automatically for each repetition.
#'
#' @examples
#' \donttest{
#' set.seed(3)
#' n_samp <- 80; n_snp <- 50
#' geno <- matrix(sample(0:2, n_samp * n_snp, replace = TRUE,
#'                       prob = c(0.6, 0.3, 0.1)),
#'                nrow = n_samp, ncol = n_snp)
#' weights  <- rnorm(n_snp)
#' leftout  <- 31:80
#' pgs_mat  <- compute_pgs_matrix(geno, leftout, weights)
#'
#' scan_df           <- data.frame(matrix(rnorm(50 * 12), 50, 12))
#' colnames(scan_df) <- c("age", "sex", "age_squared", "age_sex",
#'                        paste0("pc", 1:8))
#' scan_df$phenotype2 <- rnorm(50)
#'
#' p_vec <- screen_forward_max_region(
#'   geno              = geno,
#'   pgs.mat           = pgs_mat,
#'   this.df           = scan_df,
#'   num_signals       = -1,
#'   is_continuous     = TRUE,
#'   covariate_formula = paste(c("age", "sex"), collapse = " + "),
#'   skip1             = 5,
#'   skip2             = 5,
#'   min_window_size   = 2,
#'   max_window_size   = 8,
#'   isPlot            = FALSE
#' )
#' plot(seq(1, n_snp, by = 5), p_vec, type = "l",
#'      xlab = "SNP index", ylab = expression(-log[10](p)))
#' }
#' @export
screen_forward_max_region <- function(geno, pgs.mat, this.df, num_signals,
                              start.point = 1, save.directory = "simulation", this.chrome = 1,
                              min_window_size = 5, max_window_size = 100,
                              isSimulation = TRUE, this.repetition = 1, screening_round = 1,
                              isPlot = TRUE, skip1 = 100, skip2 = 5,
                              is_continuous,
                              covariate_formula = NULL,
                              signal.starts = NULL, signal.window.size = NULL) {
  if (is.null(covariate_formula)) {
    covariate_formula <- "sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
  }

  N <- ncol(geno)
  n <- nrow(geno)

  overall.min.pvalues <- c()
  this.seq <- seq(1, N, by = skip1)

  p.values <- rep(NA, screening_round * (length(this.seq)))

  for (kk in 1:screening_round) {

    big.count <- 0

    for (this.start in this.seq) {
      big.count <- big.count + 1

      sub.seq <- seq(min_window_size, max_window_size, by = skip2)
      sub.seq <- c(0, sub.seq)
      if (sub.seq[length(sub.seq)] != max_window_size) {
        sub.seq <- c(sub.seq, max_window_size)
      }

      sub.p.values <- rep(NA, length(sub.seq))

      count <- 0
      this.left.0 <- this.start
      this.right.0 <- this.start

      for (ws in sub.seq) {
        count <- count + 1

        if (kk == 1 && count == 1) {
          this.pgs <- pgs.mat[, (this.start)]
        } else {
          this.left.1 <- this.start - ws + 1
          this.left.1[this.left.1 < 1] <- 1

          this.right.1 <- this.start + ws - 1
          this.right.1[this.right.1 > N] <- N

          if (this.left.1 != this.left.0) {
            this.left.seq <- this.left.1:(this.left.0 - 1)
            for (j in this.left.seq) {
              this.pgs <- this.pgs + pgs.mat[, j]
            }
          }

          if (this.right.1 != this.right.0) {
            this.right.seq <- (this.right.0 + 1):this.right.1
            for (j in this.right.seq) {
              this.pgs <- this.pgs + pgs.mat[, j]
            }
          }

          this.left.0 <- this.left.1
          this.right.0 <- this.right.1
        }

        this.df$this.pgs <- this.pgs

        full_formula <- as.formula(paste("phenotype2 ~ this.pgs +", covariate_formula))
        if (is_continuous) {
          fit <- lm(full_formula, data = this.df)
        } else {
          fit <- glm(full_formula, family = binomial(), data = this.df)
        }

        if (nrow(summary(fit)$coefficients) > 1) {
          sub.p.values[(kk - 1) * length(this.seq) + count] <- summary(fit)$coefficients[2, 4]
        } else {
          print("warning")
          sub.p.values[(kk - 1) * length(this.seq) + count] <- 1
        }
      }

      p.values[big.count] <- min(sub.p.values)
    }
  }

  if (isPlot) {
    if (isSimulation) {
      pdf(file = paste0(save.directory, "/forward-", this.repetition, "-", kk, ".pdf"))
      plot(this.seq, -log(p.values[((kk - 1) * length(this.seq) + 1):(kk * length(this.seq))], base = 10),
           ylab = "-log(p.values)", type = "l", xlab = "Position",
           main = paste0("Forward, Round ", kk))
      if (num_signals >= 1) {
        abline(v = signal.starts, col = "red", lty = 2)
        abline(v = signal.starts + signal.window.size - 1, col = "red", lty = 2)
      }
    } else {
      pdf(file = paste0(save.directory, "/BPM_chr", this.chrome, "-", this.start,
                        "_forward-", this.repetition, "-", kk, ".pdf"))
      plot(this.seq, -log(p.values[((kk - 1) * length(this.seq) + 1):(kk * length(this.seq))], base = 10),
           ylab = "-log(p.values)", type = "l", xlab = "Position",
           main = paste0("BPM Chrome-", this.chrome, " Starting from ", this.start,
                         "-- Forward, Round ", kk))
    }
    dev.off()
  }

  return(-log(p.values, base = 10))
}
