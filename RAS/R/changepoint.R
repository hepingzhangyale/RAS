#' Detect a Single Changepoint via Segmented Regression
#'
#' Fits a segmented linear model to detect one breakpoint in the relationship
#' between \code{x} and \code{y} up to index \code{t}, and tests its
#' significance using the Davies test.
#'
#' @param x Numeric vector. Predictor (e.g., SNP position indices).
#' @param y Numeric vector. Response (e.g., \eqn{-\log_{10}(p)}-values from
#'   a scan). Must be the same length as \code{x}.
#' @param t Integer. Number of observations to use from the start of \code{x}
#'   and \code{y}. Must satisfy \code{t >= 4} for segmented regression to be
#'   identifiable.
#'
#' @details
#' The function first fits an ordinary linear model \code{y ~ x} over the
#' first \code{t} observations, then calls \code{\link[segmented]{segmented}}
#' to estimate a single breakpoint.  If the segmented fit succeeds, the Davies
#' test (\code{\link[segmented]{davies.test}}) is applied to the base linear
#' model to assess whether the slope change is significant.
#'
#' Both the \code{segmented} call and the Davies test are wrapped in
#' \code{try(..., silent = TRUE)}: if either fails (e.g., due to collinear
#' data or a degenerate window), the function returns \code{p.values = 1} and
#' \code{NULL} for the breakpoint and slopes rather than propagating an error.
#'
#' The breakpoint position is returned as an index into the full \code{x}
#' vector (not just the first \code{t} elements), found by locating the
#' element of \code{x[1:t]} closest to the estimated breakpoint coordinate.
#'
#' @return
#' A named list with four elements:
#' \describe{
#'   \item{\code{break.points}}{Integer. Index of the estimated breakpoint in
#'     \code{x[1:t]}, or \code{NULL} if no breakpoint was found.}
#'   \item{\code{p.values}}{Numeric. Davies test p-value for the breakpoint,
#'     or \code{1} if the fit failed or no breakpoint was found.}
#'   \item{\code{slope.left}}{Numeric. Estimated slope to the left of the
#'     breakpoint, or \code{NULL} if not found.}
#'   \item{\code{slope.right}}{Numeric. Estimated slope to the right of the
#'     breakpoint, or \code{NULL} if not found.}
#' }
#'
#' @references
#' Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
#' present only under the alternative. \emph{Biometrika}, \strong{74}(1),
#' 33--43. \doi{10.1093/biomet/74.1.33}
#'
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. \emph{Statistics in Medicine}, \strong{22}(19),
#' 3055--3071. \doi{10.1002/sim.1545}
#'
#' @seealso
#' \code{\link{ras_detect}} which calls this
#' function repeatedly in a sliding-window loop.
#' \code{\link{slope_test}} for the one-tailed slope verification step.
#' \code{\link[segmented]{segmented}}, \code{\link[segmented]{davies.test}}
#' for the underlying segmented-regression routines.
#'
#' @examples
#' set.seed(1)
#' x <- 1:60
#' y <- c(seq(0, 6, length.out = 30),
#'        seq(6, 2, length.out = 30)) + rnorm(60, sd = 0.4)
#' result <- get_break_points(x, y, t = 60)
#' cat("Breakpoint index:", result$break.points, "\n")
#' cat("Davies p-value:  ", result$p.values,    "\n")
#' cat("Left slope:      ", result$slope.left,   "\n")
#' cat("Right slope:     ", result$slope.right,  "\n")
#' @importFrom segmented segmented slope davies.test
#' @export
get_break_points <- function(x, y, t) {
  this.df <- data.frame(y = y[1:t], x = x[1:t])
  fit_lm <- lm(y ~ x, data = this.df)

  fit_segmented <- try({
    suppressWarnings(segmented(fit_lm, seg.Z = ~x, npsi = 1))
  }, silent = TRUE)

  if (inherits(fit_segmented, "try-error")) {
    return(list(break.points = NULL, p.values = 1, slope.left = NULL, slope.right = NULL))
  }

  if (is.null(fit_segmented$psi)) {
    return(list(break.points = NULL, p.values = 1, slope.left = NULL, slope.right = NULL))
  }

  this.slopes <- slope(fit_segmented)$x[, 1]
  tau_hat <- which.min(abs(x[1:t] - fit_segmented$psi[1, 2]))

  # Davies test for existence of a breakpoint
  this.test <- try(davies.test(fit_lm), silent = TRUE)
  if (inherits(this.test, "try-error")) {
    return(list(break.points = NULL, p.values = 1, slope.left = NULL, slope.right = NULL))
  }
  this.p.value <- as.numeric(this.test$p.value)
  if (!is.finite(this.p.value) || is.na(this.p.value)) this.p.value <- 1

  return(list(break.points = tau_hat, p.values = this.p.value,
              slope.left = this.slopes[1], slope.right = this.slopes[2]))
}


#' One-Tailed Slope Test Through the Origin
#'
#' Fits a no-intercept linear model \code{y ~ x - 1} and performs a
#' one-tailed t-test on the slope coefficient.
#'
#' @param x Numeric vector. Predictor values, typically centred at the
#'   candidate changepoint so that the origin corresponds to the changepoint
#'   position.
#' @param y Numeric vector. Response values, same length as \code{x},
#'   typically centred at the changepoint y-value.
#' @param lower.tail Logical. Passed to \code{\link[stats]{pt}}.  Use
#'   \code{FALSE} to test \eqn{H_0\!: \beta \ge 0} against
#'   \eqn{H_1\!: \beta < 0} (left-slope test); use \code{TRUE} for the
#'   opposite direction (right-slope test).
#'
#' @details
#' The model \code{lm(y ~ x - 1)} forces the regression line through the
#' origin, which is appropriate after centring both \code{x} and \code{y}
#' at the candidate changepoint.  Under this parameterisation a positive
#' left slope and a negative right slope correspond to the peak shape
#' expected at a true association region.
#'
#' The function is called by \code{\link{changePoint_detection_window_scanning}}
#' on both sides of each Davies-significant candidate, using asymmetric
#' thresholds (\code{cp_slope_left = 1e-10}, \code{cp_slope_right = 1e-20}
#' by default) to require a steeper descending edge than ascending edge.
#'
#' @return
#' Numeric scalar. The one-tailed p-value for the slope coefficient.
#'
#' @seealso
#' \code{\link{get_break_points}} for the Davies test step that precedes this
#' slope check. \code{\link{ras_detect}} for the
#' full detection workflow.
#'
#' @examples
#' # Left-side slope: x goes from negative to 0, y should be rising (positive slope)
#' x_left <- -10:0
#' y_left <- x_left * 0.8 + rnorm(11, sd = 0.3)
#' slope_test(x_left, y_left, lower.tail = FALSE)  # expect small p (slope > 0)
#'
#' # Right-side slope: x goes from 0 to positive, y should be falling (negative slope)
#' x_right <- 0:10
#' y_right <- x_right * (-0.8) + rnorm(11, sd = 0.3)
#' slope_test(x_right, y_right, lower.tail = TRUE)  # expect small p (slope < 0)
#' @export
slope_test <- function(x, y, lower.tail) {
  model <- lm(y ~ x - 1)
  t_statistic <- coef(summary(model))["x", "t value"]
  p_value <- pt(t_statistic, df = df.residual(model), lower.tail = lower.tail)
  p_value
}


#' First-Pass Changepoint Detection via Sliding Window
#'
#' Scans a \eqn{-\log_{10}(p)}-value sequence using a sliding window to
#' detect positions where the slope changes significantly from positive to
#' negative.
#'
#' @param x Numeric vector. Predictor sequence (e.g., SNP position indices).
#' @param y Numeric vector. Response sequence (e.g., \eqn{-\log_{10}(p)}
#'   values from \code{\link{run_ras_scan}}).  Must be the same length as
#'   \code{x}.
#' @param p.values.threshold Numeric. Davies test p-value threshold for
#'   nominating a candidate changepoint.  Default \code{0.01}.
#' @param min.length Integer. Minimum number of observations required on each
#'   side of a candidate changepoint.  Default \code{10}.
#' @param skip Integer. Step size when iterating the sliding window start
#'   position.  Default \code{1}.
#' @param window_size Integer. Number of observations in each sliding window.
#'   Default \code{3000}.
#' @param slope_check_window_size Integer. Half-width of the local region used
#'   to verify slope direction via \code{\link{slope_test}}.  Default
#'   \code{30}.
#' @param slope.p.values.threshold Numeric. Reserved combined slope threshold
#'   (currently unused in filtering).  Default \code{1e-8}.
#' @param slope.p.values.threshold.left Numeric. One-tailed p-value threshold
#'   for the left-side slope test.  A candidate is accepted only if the slope
#'   to its left is significantly positive (\eqn{p <} this value).  Default
#'   \code{1e-10}.
#' @param slope.p.values.threshold.right Numeric. One-tailed p-value threshold
#'   for the right-side slope test.  A candidate is accepted only if the slope
#'   to its right is significantly negative (\eqn{p <} this value).  Default
#'   \code{1e-20}.
#'
#' @details
#' At each window start position the function calls
#' \code{\link{get_break_points}} on the window to estimate and test a single
#' breakpoint.  A candidate is retained when three conditions are all met:
#' \enumerate{
#'   \item The Davies test p-value is below \code{p.values.threshold}.
#'   \item The left-side slope (estimated by segmented regression) is
#'     positive.
#'   \item The right-side slope is negative.
#' }
#' Candidates passing these three filters are then subjected to one-tailed
#' slope tests (\code{\link{slope_test}}) using centred sub-sequences of
#' half-width \code{slope_check_window_size}.  Only candidates that also pass
#' both slope p-value thresholds are recorded.
#'
#' After the sliding-window loop, each accepted changepoint is refined to the
#' nearest local peak in \code{y} via \code{\link{get_local_maximum}}.
#'
#' The function records all candidates examined (\code{all.changepoints},
#' \code{all.p.values}) in addition to the accepted ones, so that downstream
#' functions can display the full candidate landscape.
#'
#' @return
#' A named list with eight elements:
#' \describe{
#'   \item{\code{tau_hats}}{Integer vector. Accepted changepoint indices
#'     (refined to local peaks).}
#'   \item{\code{p.values}}{Numeric vector. Davies test p-values at accepted
#'     changepoints.}
#'   \item{\code{slope.left}}{Numeric vector. Left-side slopes at accepted
#'     changepoints.}
#'   \item{\code{slope.right}}{Numeric vector. Right-side slopes at accepted
#'     changepoints.}
#'   \item{\code{all.changepoints}}{Integer vector. All candidate positions
#'     examined, including those that failed the slope tests.}
#'   \item{\code{all.p.values}}{Numeric vector. Davies p-values for all
#'     examined candidates; set to \code{1} for candidates that failed the
#'     slope direction or slope significance filters.}
#'   \item{\code{slope.angle}}{Numeric vector. Interior angle (degrees) at
#'     each accepted changepoint, computed from the left and right slope
#'     estimates.}
#'   \item{\code{previous_tau_hats}}{Integer vector. Copy of \code{tau_hats}
#'     before any downstream modification; used by
#'     \code{\link{second_scanning}}.}
#' }
#'
#' @references
#' Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
#' present only under the alternative. \emph{Biometrika}, \strong{74}(1),
#' 33--43. \doi{10.1093/biomet/74.1.33}
#'
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. \emph{Statistics in Medicine}, \strong{22}(19),
#' 3055--3071. \doi{10.1002/sim.1545}
#'
#' @seealso
#' \code{\link{ras_validate}} for the second-pass validation step.
#' \code{\link{get_break_points}} for the per-window segmented regression.
#' \code{\link{slope_test}} for the one-tailed slope verification.
#' \code{\link{get_local_maximum}} for the peak-refinement step.
#' \code{\link{ras}} for the recommended end-to-end entry point.
#'
#' @examples
#' set.seed(42)
#' x <- 1:300
#' y <- c(seq(0, 8, length.out = 150),
#'        seq(8, 1, length.out = 150)) + rnorm(300, sd = 0.5)
#' result <- ras_detect(
#'   x, y,
#'   window_size             = 150,
#'   slope_check_window_size = 20,
#'   slope.p.values.threshold.left  = 1e-3,
#'   slope.p.values.threshold.right = 1e-3
#' )
#' cat("Detected changepoints:", result$tau_hats, "\n")
#' @export
ras_detect <- function(x, y, p.values.threshold = 0.01,
                                                  min.length = 10, skip = 1,
                                                  window_size = 3000,
                                                  slope_check_window_size = 30,
                                                  slope.p.values.threshold = 1e-8,
                                                  slope.p.values.threshold.left = 1e-10,
                                                  slope.p.values.threshold.right = 1e-20) {
  tau_hats <- c()
  p.values <- c()
  slope.left <- c()
  slope.right <- c()
  end.index <- c()
  all.changepoints <- c()
  all.p.values <- c()
  slope.angle <- c()
  start.t <- min.length - window_size

  while (TRUE) {
    this.seq <- seq(start.t, length(x) - min.length + 1, by = skip)

    if (this.seq[length(this.seq)] != length(x) - min.length + 1) {
      this.seq <- c(this.seq, length(x) - min.length + 1)
    }

    for (start.t in this.seq) {
      t <- window_size

      if (start.t < 0) {
        t <- window_size + start.t
        start.t <- 1
      }

      if (start.t > (length(x) - window_size + 1)) {
        t <- length(x) - start.t + 1
      }

      this.result <- get_break_points(x[start.t:length(x)], y[start.t:length(x)], t)

      if (!is.null(this.result$break.points)) {
        this.result$break.points <- start.t - 1 + this.result$break.points
        all.changepoints <- c(all.changepoints, this.result$break.points)
        all.p.values <- c(all.p.values, as.numeric(this.result$p.values))
      }

      if (is.null(this.result$break.points)) next

      v1 <- this.result$slope.left > 0
      v2 <- this.result$slope.right < 0

      if (this.result$p.values <= p.values.threshold) {
        if (v1 && v2) {
          this.tau <- this.result$break.points
          this.lower.index <- max(1, this.tau - slope_check_window_size)
          this.upper.index <- min(length(x), this.tau + slope_check_window_size)

          this.pvalue1 <- slope_test(x[this.lower.index:this.tau] - x[this.tau],
                                     y[this.lower.index:this.tau] - y[this.tau], FALSE)
          this.pvalue2 <- slope_test(x[this.tau:this.upper.index] - x[this.tau],
                                     y[this.tau:this.upper.index] - y[this.tau], TRUE)

          if (this.pvalue1 < slope.p.values.threshold.left &&
              this.pvalue2 < slope.p.values.threshold.right) {
            cat("ChangePoint Detected at ", this.result$break.points,
                " with p-value", this.result$p.values,
                ", Left-slope ", this.result$slope.left,
                "Right-slope ", this.result$slope.right, " \n")
            cat(this.pvalue1, " ", this.pvalue2, "\n")

            tau_hats <- c(tau_hats, this.result$break.points)
            p.values <- c(p.values, this.result$p.values)
            slope.left <- c(slope.left, this.result$slope.left)
            slope.right <- c(slope.right, this.result$slope.right)

            angle_degrees1 <- (atan(this.result$slope.left) * 180) / pi
            angle_degrees2 <- (atan(this.result$slope.right) * 180) / pi
            theta_degrees <- angle_degrees2 - angle_degrees1 + 180
            slope.angle <- c(slope.angle, theta_degrees)

            end.index <- c(end.index, start.t + t - 1)
            start.t <- this.result$break.points
            break
          } else {
            all.p.values[length(all.p.values)] <- 1
          }
        } else {
          all.p.values[length(all.p.values)] <- 1
        }
      }
    }

    if (start.t == this.seq[length(this.seq)]) {
      break
    }
  }

  # Refine each detected changepoint to the nearest local peak in y
  tau_hats <- sapply(tau_hats, get_local_maximum, y = y)

  return(list(tau_hats = tau_hats, p.values = p.values,
              slope.left = slope.left, slope.right = slope.right,
              all.changepoints = all.changepoints, all.p.values = all.p.values,
              slope.angle = slope.angle, previous_tau_hats = tau_hats))
}


#' Second-Pass Changepoint Validation
#'
#' Validates candidate changepoints from
#' \code{\link{changePoint_detection_window_scanning}} by re-running local
#' Davies tests in windows around each candidate.
#'
#' @param this.result List. Output from
#'   \code{\link{changePoint_detection_window_scanning}}.
#' @param x Numeric vector. Predictor sequence used in the original scan.
#' @param y Numeric vector. Response sequence used in the original scan.
#' @param this.start Integer. Genomic start position for index re-mapping to
#'   chromosome coordinates.  Default \code{1}.
#' @param this.skip Integer. Step size used in the original scan
#'   (\code{skip1}), needed for index re-mapping.  Default \code{30}.
#' @param second_window_size Integer. Half-window size for local Davies test
#'   re-validation.  Default \code{50}.
#' @param p.value.threshold Numeric. Davies test p-value threshold for
#'   second-pass acceptance.  A candidate is retained if \emph{either} the
#'   right-side or left-side local Davies test passes.  Default \code{1e-10}.
#'
#' @details
#' For each candidate \eqn{\hat{\tau}} in \code{this.result$tau_hats} the
#' function fits two local linear models and applies the Davies test to each:
#' one on the window
#' \eqn{[\hat{\tau},\; \min(\hat{\tau} + \texttt{second\_window\_size},\, n)]}
#' and one on
#' \eqn{[\max(\hat{\tau} - \texttt{second\_window\_size},\, 1),\; \hat{\tau}]}.
#' If either Davies p-value is below \code{p.value.threshold}, the candidate
#' is accepted; otherwise its \code{all.p.values} entry is set to zero to
#' suppress it in downstream plots.
#'
#' Windows with fewer than four observations are not tested (Davies test
#' requires at least four points) and their p-value is set to \code{1.0}.
#'
#' After filtering, any remaining candidate with \code{y[tau_hat] <= 2.5} is
#' removed, as such positions are below the minimum signal threshold.
#'
#' Accepted changepoint indices and all candidate indices are re-mapped from
#' the scan grid to genomic coordinates via
#' \eqn{\texttt{this.start} + (\text{index} - 1) \times \texttt{this.skip}}.
#'
#' @return
#' A named list with five elements:
#' \describe{
#'   \item{\code{tau_hats}}{Integer vector. Validated changepoint positions,
#'     re-mapped to genomic coordinates.}
#'   \item{\code{all.changepoints}}{Integer vector. All candidate positions
#'     examined, re-mapped to genomic coordinates.}
#'   \item{\code{all.p.values}}{Numeric vector. \eqn{-\log_{10}(p)} Davies
#'     values for each candidate; set to \code{0} for rejected candidates.}
#'   \item{\code{left.slopes}}{Numeric vector. Left-side slopes at accepted
#'     changepoints (carried over from first-pass output).}
#'   \item{\code{right.slopes}}{Numeric vector. Right-side slopes at accepted
#'     changepoints.}
#' }
#'
#' @references
#' Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
#' present only under the alternative. \emph{Biometrika}, \strong{74}(1),
#' 33--43. \doi{10.1093/biomet/74.1.33}
#'
#' @seealso
#' \code{\link{ras_detect}} for the first-pass
#' detection step whose output this function takes as input.
#' \code{\link{plot.ras}} for visualising the validated changepoints.
#' \code{\link{ras}} for the recommended end-to-end entry point.
#'
#' @examples
#' set.seed(42)
#' x <- 1:300
#' y <- c(seq(0, 8, length.out = 150),
#'        seq(8, 1, length.out = 150)) + rnorm(300, sd = 0.5)
#'
#' cp_result <- ras_detect(
#'   x, y,
#'   window_size             = 150,
#'   slope_check_window_size = 20,
#'   slope.p.values.threshold.left  = 1e-3,
#'   slope.p.values.threshold.right = 1e-3
#' )
#'
#' final <- ras_validate(
#'   cp_result, x = x, y = y,
#'   this.skip          = 1,
#'   second_window_size = 30,
#'   p.value.threshold  = 1e-3
#' )
#' cat("Validated changepoints:", final$tau_hats, "\n")
#' @importFrom segmented davies.test
#' @export
ras_validate <- function(this.result, x, y, this.start = 1, this.skip = 30,
                            second_window_size = 50, p.value.threshold = 1e-10) {

  unique_all.changepoints <- unique(this.result$all.changepoints)

  if (length(unique_all.changepoints) == 0) {
    return(list(all.changepoints = c(), tau_hats = c(), all.p.values = c(),
                left.slopes = c(), right.slopes = c()))
  }

  min_all.p.values <- sapply(unique_all.changepoints,
                             function(cp) {
                               vals <- as.numeric(this.result$all.p.values[this.result$all.changepoints == cp])
                               min(vals, na.rm = TRUE)
                             })
  tau_hats <- c()
  left.slopes <- c()
  right.slopes <- c()

  all.changepoints <- unique_all.changepoints
  all.p.values <- -log(pmax(min_all.p.values, .Machine$double.eps), base = 10)

  this.count <- 0

  for (tau_hat in this.result$tau_hats) {
    this.count <- this.count + 1

    this.df <- data.frame(
      y = y[tau_hat:min((tau_hat + second_window_size), length(y))],
      x = x[tau_hat:min((tau_hat + second_window_size), length(x))]
    )
    fit_lm <- lm(y ~ x, data = this.df)

    if (dim(this.df)[1] >= 4) {
      this.test1 <- davies.test(fit_lm)
    } else {
      this.test1 <- list(p.value = 1.0)
    }

    this.df <- data.frame(
      y = y[max((tau_hat - second_window_size), 1):tau_hat],
      x = x[max((tau_hat - second_window_size), 1):tau_hat]
    )
    fit_lm <- lm(y ~ x, data = this.df)

    if (dim(this.df)[1] >= 4) {
      this.test2 <- davies.test(fit_lm)
    } else {
      this.test2 <- list(p.value = 1.0)
    }

    cat(this.test1$p.value, this.test2$p.value, "\n")

    if (this.test1$p.value < p.value.threshold || this.test2$p.value < p.value.threshold) {
      left.slopes <- c(left.slopes, this.result$slope.left[which(tau_hat == this.result$tau_hats)])
      right.slopes <- c(right.slopes, this.result$slope.right[which(tau_hat == this.result$tau_hats)])
      tau_hats <- c(tau_hats, tau_hat)
    } else {
      all.p.values[which(tau_hat == all.changepoints)] <- 0
      all.p.values[which(this.result$previous_tau_hats[this.count] == all.changepoints)] <- 0
    }
  }

  this.remove <- which(y[tau_hats] <= 2.5)
  if (length(this.remove) >= 1) {
    tau_hats <- tau_hats[-this.remove]
  }

  all.changepoints <- this.start + (all.changepoints - 1) * this.skip
  tau_hats <- this.start + (tau_hats - 1) * this.skip

  return(list(all.changepoints = all.changepoints, tau_hats = tau_hats,
              all.p.values = all.p.values, left.slopes = left.slopes,
              right.slopes = right.slopes))
}
