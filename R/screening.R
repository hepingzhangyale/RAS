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
#' @param is_continuous Logical. If \code{TRUE}, the per-window PGS p-value is
#'   obtained from a linear model via an exact Frisch-Waugh-Lovell
#'   residualisation (numerically identical to \code{lm(phenotype2 ~ this.pgs +
#'   covariates)} but the covariate factorisation is done once instead of per
#'   window).  If \code{FALSE}, a binary trait is used (see \code{scan_test}).
#' @param scan_test Character. Test used for the binary branch
#'   (\code{is_continuous = FALSE}).  \code{"glm"} (default) fits a full
#'   logistic regression per window and reports the Wald p-value of the PGS
#'   term, preserving the original behaviour.  \code{"score"} fits the
#'   covariate-only null logistic model \strong{once} and evaluates a Rao
#'   score test for the PGS term at each window, which is asymptotically
#'   equivalent but avoids per-window iterative fitting (typically 1-2 orders
#'   of magnitude faster at scale).  Ignored when \code{is_continuous = TRUE}
#'   (the continuous branch always uses \code{lm}).
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
                              scan_test = c("glm", "score"),
                              signal.starts = NULL, signal.window.size = NULL) {
  if (is.null(covariate_formula)) {
    covariate_formula <- "sex + age + age_squared + age_sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
  }
  scan_test <- match.arg(scan_test)

  N <- ncol(geno)
  n <- nrow(geno)

  # ── Score-test pre-computation (binary branch only) ───────────────────────
  # The Rao score test for adding the PGS term to a logistic model needs only
  # the covariate-only NULL fit, which is identical across all windows. Fit it
  # ONCE here and cache the quantities each window reuses:
  #   r0 = y - p0      (null residual; score numerator is pgs' r0)
  #   w0 = p0(1 - p0)  (null working weights)
  #   Mmat = (Z'W0 Z)^{-1} Z'W0  so the W0-orthogonalised pgs is
  #          g_perp = g - Z (Mmat g), and the score stat is
  #          T = (g' r0)^2 / sum(w0 * g_perp^2) ~ chi-sq(1).
  score_mode <- (!is_continuous && scan_test == "score")
  if (score_mode) {
    cov_form <- stats::as.formula(paste("~", covariate_formula))
    Zmat <- stats::model.matrix(cov_form, data = this.df)   # incl. intercept
    yvec <- this.df$phenotype2
    keep <- stats::complete.cases(Zmat) & !is.na(yvec)
    Zmat <- Zmat[keep, , drop = FALSE]
    yvec <- as.numeric(yvec[keep])
    fit0 <- stats::glm.fit(Zmat, yvec, family = stats::binomial())
    p0   <- fit0$fitted.values
    w0   <- p0 * (1 - p0)
    r0   <- yvec - p0
    Mmat <- solve(crossprod(Zmat, Zmat * w0), t(Zmat * w0))  # (Z'W0Z)^-1 Z'W0
  }

  # ── FWL pre-computation (continuous branch) ───────────────────────────────
  # By Frisch-Waugh-Lovell, the PGS coefficient (and its SE/t/p) in the full
  # lm(phenotype2 ~ this.pgs + covariates) equals the simple regression of the
  # covariate-residualised phenotype on the covariate-residualised PGS, using
  # the full-model residual df = n - rank(Z) - 1.  The covariate design Z is
  # identical across all windows, so factor it ONCE: cache qr(Z), the
  # residualised phenotype ryC and Syy.  Per window only the PGS is residualised
  # (one qr.resid).  Output is exact-equal to the original per-window lm.
  cont_mode <- is_continuous
  if (cont_mode) {
    cov_formC <- stats::as.formula(paste("~", covariate_formula))
    ZmatC <- stats::model.matrix(cov_formC, data = this.df)  # incl. intercept
    yvecC <- this.df$phenotype2
    keepC <- stats::complete.cases(ZmatC) & !is.na(yvecC)
    ZmatC <- ZmatC[keepC, , drop = FALSE]
    yvecC <- as.numeric(yvecC[keepC])
    qrZC  <- qr(ZmatC)
    ryC   <- qr.resid(qrZC, yvecC)
    SyyC  <- sum(ryC * ryC)
    dfC   <- length(yvecC) - (qrZC$rank + 1L)   # n - rank(full design)
  }

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

        if (score_mode) {
          # Rao score test: orthogonalise this window's pgs against the
          # covariates in the null-fit W0 metric, then T = U^2 / V ~ chi-sq(1).
          g     <- this.pgs[keep]
          gperp <- g - as.numeric(Zmat %*% (Mmat %*% g))
          U     <- sum(g * r0)               # = sum(gperp * r0), since Z'r0 = 0
          V     <- sum(w0 * gperp * gperp)
          if (is.finite(V) && V > 0) {
            Tstat <- (U * U) / V
            sub.p.values[(kk - 1) * length(this.seq) + count] <-
              stats::pchisq(Tstat, df = 1, lower.tail = FALSE)
          } else {
            sub.p.values[(kk - 1) * length(this.seq) + count] <- 1
          }
        } else if (cont_mode) {
          # Exact FWL: residualise this window's pgs against the covariates,
          # then a simple regression on residuals reproduces lm's pgs p-value.
          g    <- this.pgs[keepC]
          rg   <- qr.resid(qrZC, g)
          Srxx <- sum(rg * rg)
          if (is.finite(Srxx) && Srxx > 0 && dfC > 0L) {
            Srxy <- sum(rg * ryC)
            beta <- Srxy / Srxx
            RSS  <- max(SyyC - beta * Srxy, 0)
            SE   <- sqrt((RSS / dfC) / Srxx)
            tval <- beta / SE
            sub.p.values[(kk - 1) * length(this.seq) + count] <-
              2 * stats::pt(-abs(tval), df = dfC)
          } else {
            sub.p.values[(kk - 1) * length(this.seq) + count] <- 1
          }
        } else {
          this.df$this.pgs <- this.pgs

          full_formula <- as.formula(paste("phenotype2 ~ this.pgs +", covariate_formula))
          fit <- glm(full_formula, family = binomial(), data = this.df)

          if (nrow(summary(fit)$coefficients) > 1) {
            sub.p.values[(kk - 1) * length(this.seq) + count] <- summary(fit)$coefficients[2, 4]
          } else {
            print("warning")
            sub.p.values[(kk - 1) * length(this.seq) + count] <- 1
          }
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
