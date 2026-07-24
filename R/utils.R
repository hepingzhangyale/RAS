#' Estimate Memory Requirements and Check System Readiness
#'
#' Detects the current machine's RAM and CPU configuration, computes
#' per-stage peak memory estimates for the RAS pipeline, and issues a
#' go / no-go verdict.
#'
#' @param n_total Integer. Total number of samples (rows of \code{geno}).
#' @param n_train Integer. Number of training-split samples (typically
#'   \code{n_total / 2}).
#' @param n_holdout Integer. Number of hold-out samples (rows of
#'   \code{pgs.mat}; typically \code{n_total - n_train}).
#' @param n_snps Integer. Number of variants (columns of \code{geno}).
#' @param bytes_per_element Numeric. Bytes per matrix element.  Default
#'   \code{8} (R \code{double}).
#' @param abort Logical. If \code{TRUE} and estimated peak memory exceeds
#'   available RAM, calls \code{\link[base]{stop}} so the pipeline cannot
#'   proceed.  Default \code{FALSE}.
#'
#' @details
#' Peak memory is estimated for three pipeline stages:
#' \describe{
#'   \item{Stage 1 — \code{\link{compute_gwas_weights}}}{Holds the full
#'     \code{geno} matrix plus a small \eqn{N \times 4} coefficient matrix.
#'     Peak \eqn{\approx} \code{geno_mb + coefmat_mb}.}
#'   \item{Stage 2 — \code{\link{compute_pgs_matrix}}}{Worst-case holds two
#'     copies of \code{geno} (R copy-on-modify triggered by
#'     \code{geno[is.na(geno)] <- 0}) plus the output \code{pgs.mat}.
#'     Peak \eqn{\approx} \code{2 * geno_mb + pgsmat_mb}.}
#'   \item{Stage 3 — \code{\link{screen_forward_max_region}}}{Holds
#'     \code{geno} and \code{pgs.mat} simultaneously.
#'     Peak \eqn{\approx} \code{geno_mb + pgsmat_mb}.}
#' }
#' The overall estimated peak is the maximum across the three stages.
#'
#' Available RAM is queried via \code{wmic} on Windows and
#' \code{/proc/meminfo} on Linux; the function degrades gracefully (prints a
#' caution message and returns \code{can_proceed = TRUE}) if the query fails.
#'
#' The function prints a formatted report to the console and invisibly returns
#' the numeric estimates for programmatic use.
#'
#' @return
#' Invisibly returns a named list with two elements:
#' \describe{
#'   \item{\code{memory}}{Named numeric vector with elements
#'     \code{geno_mb}, \code{pgsmat_mb}, \code{coefmat_mb},
#'     \code{stage1_peak_mb}, \code{stage2_peak_mb}, \code{stage3_peak_mb},
#'     \code{overall_peak_mb}, \code{available_mb}, \code{total_mb}.}
#'   \item{\code{system}}{Named list with elements \code{cores_physical},
#'     \code{cores_logical}, \code{r_version}, \code{platform}, and
#'     \code{can_proceed} (logical: \code{TRUE} if RAM is sufficient or
#'     cannot be measured).}
#' }
#'
#' @seealso
#' \code{\link{ras_scan}}, \code{\link{ras}} for the
#' functions whose memory use is being estimated.
#' \code{\link{release_memory}} to reclaim heap memory after each repetition.
#'
#' @examples
#' ## Check feasibility for 5,000 samples and 500,000 SNPs
#' ras_memory(
#'   n_total   = 5000,
#'   n_train   = 2500,
#'   n_holdout = 2500,
#'   n_snps    = 500000
#' )
#'
#' ## Abort if memory is insufficient (wrapped in try() so the example runs)
#' try(ras_memory(
#'   n_total   = 100000,
#'   n_train   = 50000,
#'   n_holdout = 50000,
#'   n_snps    = 1000000,
#'   abort     = TRUE
#' ))
#' @export
ras_memory <- function(n_total, n_train, n_holdout, n_snps,
                                bytes_per_element = 8,
                                abort = FALSE) {
  to_mb <- function(bytes) bytes / 1024^2

  # ── Memory estimates ──────────────────────────────────────────────────────
  geno_mb    <- to_mb(n_total   * n_snps * bytes_per_element)
  pgsmat_mb  <- to_mb(n_holdout * n_snps * bytes_per_element)
  coefmat_mb <- to_mb(n_snps    * 4      * bytes_per_element)

  stage1_mb <- geno_mb + coefmat_mb
  stage2_mb <- 2 * geno_mb + pgsmat_mb   # copy-on-modify worst case
  stage3_mb <- geno_mb + pgsmat_mb

  overall_peak_mb <- max(stage1_mb, stage2_mb, stage3_mb)

  # ── System detection ──────────────────────────────────────────────────────
  free_mb  <- .get_free_memory_mb()
  total_mb <- .get_total_memory_mb()

  cores_physical <- tryCatch(
    parallel::detectCores(logical = FALSE),
    error = function(e) NA_integer_
  )
  cores_logical <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) NA_integer_
  )

  r_ver    <- paste(R.version$major, R.version$minor, sep = ".")
  platform <- R.version$platform

  # ── Verdict ───────────────────────────────────────────────────────────────
  # can_proceed = TRUE if we cannot measure RAM (unknown) or if there is enough
  can_proceed <- is.na(free_mb) || (overall_peak_mb <= free_mb)

  ratio <- if (!is.na(free_mb) && free_mb > 0) overall_peak_mb / free_mb else NA_real_

  verdict_line <- if (is.na(free_mb)) {
    "  [?] Available RAM unknown -- cannot verify, proceed with caution."
  } else if (!can_proceed) {
    sprintf("  [NO-GO] Peak %.1f MB EXCEEDS available %.1f MB (%.0f%% of free RAM).",
            overall_peak_mb, free_mb, ratio * 100)
  } else if (!is.na(ratio) && ratio > 0.8) {
    sprintf("  [CAUTION] Peak uses %.0f%% of free RAM -- monitor closely.",
            ratio * 100)
  } else {
    sprintf("  [GO] Peak uses %.0f%% of free RAM -- OK to proceed.",
            ratio * 100)
  }

  # ── Print report ──────────────────────────────────────────────────────────
  cat("============================================================\n")
  cat("        RAS Pipeline - System & Memory Check\n")
  cat("============================================================\n")
  cat("  R version  :", r_ver, "\n")
  cat("  Platform   :", platform, "\n")
  if (!is.na(cores_physical)) {
    cat(sprintf("  CPU cores  : %d physical / %d logical\n",
                cores_physical, cores_logical))
  } else {
    cat("  CPU cores  : could not be detected\n")
  }
  if (!is.na(total_mb)) {
    cat(sprintf("  Total RAM  : %8.1f MB  (%.1f GB)\n", total_mb, total_mb / 1024))
  }
  if (!is.na(free_mb)) {
    cat(sprintf("  Free  RAM  : %8.1f MB  (%.1f GB)\n", free_mb, free_mb / 1024))
  } else {
    cat("  Free  RAM  : could not be determined\n")
  }
  cat("------------------------------------------------------------\n")
  cat(sprintf("  Input : %d samples x %d SNPs\n", n_total, n_snps))
  cat(sprintf("  Split : %d train / %d hold-out\n", n_train, n_holdout))
  cat("------------------------------------------------------------\n")
  cat(sprintf("  geno matrix             : %10.1f MB\n", geno_mb))
  cat(sprintf("  pgs.mat matrix          : %10.1f MB\n", pgsmat_mb))
  cat(sprintf("  coef.mat (GWAS output)  : %10.1f MB\n", coefmat_mb))
  cat("------------------------------------------------------------\n")
  cat(sprintf("  Stage 1  compute_gwas_weights  peak : %8.1f MB\n", stage1_mb))
  cat(sprintf("  Stage 2  compute_pgs_matrix    peak : %8.1f MB  (geno copied)\n",
              stage2_mb))
  cat(sprintf("  Stage 3  screen_forward_max    peak : %8.1f MB\n", stage3_mb))
  cat("------------------------------------------------------------\n")
  cat(sprintf("  Overall estimated peak               : %8.1f MB\n", overall_peak_mb))
  cat(verdict_line, "\n")
  cat("============================================================\n")

  if (!can_proceed && abort) {
    stop(sprintf(
      "RAS pipeline aborted: estimated peak (%.1f MB) exceeds available RAM (%.1f MB).",
      overall_peak_mb, free_mb
    ), call. = FALSE)
  }

  invisible(list(
    memory = c(
      geno_mb         = geno_mb,
      pgsmat_mb       = pgsmat_mb,
      coefmat_mb      = coefmat_mb,
      stage1_peak_mb  = stage1_mb,
      stage2_peak_mb  = stage2_mb,
      stage3_peak_mb  = stage3_mb,
      overall_peak_mb = overall_peak_mb,
      available_mb    = free_mb,
      total_mb        = total_mb
    ),
    system = list(
      cores_physical = cores_physical,
      cores_logical  = cores_logical,
      r_version      = r_ver,
      platform       = platform,
      can_proceed    = can_proceed
    )
  ))
}


# Internal helper: query free physical memory in MB (cross-platform).
.get_free_memory_mb <- function() {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      raw <- system("wmic OS get FreePhysicalMemory /Value", intern = TRUE)
      val <- raw[grep("FreePhysicalMemory=", raw)]
      kb  <- as.numeric(sub("FreePhysicalMemory=", "", trimws(val)))
      kb / 1024
    } else if (file.exists("/proc/meminfo")) {
      lines <- readLines("/proc/meminfo")
      avail <- lines[grep("^MemAvailable:", lines)]
      kb    <- as.numeric(gsub("[^0-9]", "", avail))
      kb / 1024
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
}

# Internal helper: query total physical memory in MB (cross-platform).
.get_total_memory_mb <- function() {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      raw <- system("wmic OS get TotalVisibleMemorySize /Value", intern = TRUE)
      val <- raw[grep("TotalVisibleMemorySize=", raw)]
      kb  <- as.numeric(sub("TotalVisibleMemorySize=", "", trimws(val)))
      kb / 1024
    } else if (file.exists("/proc/meminfo")) {
      lines <- readLines("/proc/meminfo")
      total <- lines[grep("^MemTotal:", lines)]
      kb    <- as.numeric(gsub("[^0-9]", "", total))
      kb / 1024
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
}


#' Find the Local Maximum Within a Window
#'
#' Returns the index of the maximum value of \code{y} within a symmetric
#' window of half-width \code{window.size} centred at \code{x0}.
#'
#' @param y Numeric vector. Values to search (e.g., \eqn{-\log_{10}(p)}
#'   values from the RAS scan).
#' @param x0 Integer. Centre index of the search window.
#' @param window.size Integer. Half-width of the search window.  The search
#'   covers indices \code{max(1, x0 - window.size)} to
#'   \code{min(length(y), x0 + window.size)}.  Default \code{50}.
#'
#' @details
#' This function is used by \code{\link{ras_detect}}
#' after the sliding-window loop to snap each accepted changepoint index to
#' the nearest local peak in the \eqn{-\log_{10}(p)} profile.  Snapping to
#' the peak ensures that reported positions correspond to the most significant
#' SNP in the association region rather than to the mathematical breakpoint
#' of the piecewise linear fit, which may be slightly offset.
#'
#' @return
#' Integer. The index (into \code{y}) of the local maximum within the window
#' centred at \code{x0}.  If multiple positions tie for the maximum,
#' \code{\link[base]{which.max}} returns the first.
#'
#' @seealso
#' \code{\link{ras_detect}} which calls this
#' function to refine detected changepoint positions.
#'
#' @examples
#' y <- c(1, 3, 7, 5, 2, 8, 4, 1)
#'
#' ## Peak within a window of half-width 2 centred at index 3
#' get_local_maximum(y, x0 = 3, window.size = 2)  # returns 3 (value = 7)
#'
#' ## Peak within a window of half-width 3 centred at index 3
#' ## covers indices 1:6; max is at index 6 (value = 8)
#' get_local_maximum(y, x0 = 3, window.size = 3)  # returns 6
#' @export
get_local_maximum <- function(y, x0, window.size = 50) {
  lower.bound <- max(1, x0 - window.size)
  upper.bound <- min(length(y), x0 + window.size)
  x0 <- which.max(y[lower.bound:upper.bound]) + lower.bound - 1
  return(x0)
}
