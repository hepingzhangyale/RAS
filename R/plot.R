#' Plot a RAS Result Object
#'
#' S3 plot method for objects of class \code{"ras"} returned by
#' \code{\link{ras}}.  When \code{zoom = FALSE} (default) produces the
#' full-chromosome scan profile with detected changepoints marked.  When
#' \code{zoom = TRUE} produces the zoomed multi-panel figure around each
#' detected changepoint.
#'
#' @param x Object of class \code{"ras"} as returned by \code{\link{ras}}.
#' @param zoom Logical. \code{FALSE} (default) plots the full scan profile;
#'   \code{TRUE} plots zoomed panels around each detected changepoint.
#' @param device Character. Output device: \code{"pdf"} (default),
#'   \code{"png"}, or \code{"screen"}.
#' @param p.threshold Numeric. Significance reference line on the
#'   \eqn{-\log_{10}} scale.  Default \code{8}.
#' @param y_cap Numeric or \code{NULL}. Cap the y-axis (scan plot only).
#'   Default \code{NULL}.
#' @param min_display_p Numeric. Minimum \code{all.p.values} for a candidate
#'   to appear in the dual-axis overlay (scan plot only).  Default \code{1}.
#' @param xlim Numeric vector of length 2, or \code{NULL}. Restrict plot to
#'   this genomic range (scan plot only).  Default \code{NULL}.
#' @param zoom_half_width Numeric. Half-width of each zoom panel in SNP index
#'   units (zoom plot only).  Default \code{3000}.
#' @param ncol Integer. Columns in the zoom panel grid (zoom plot only).
#'   Default \code{3}.
#' @param min_signal Numeric. Low-signal colour boundary on the
#'   \eqn{-\log_{10}} scale: scan values below this threshold are drawn in
#'   blue; values at or above it are drawn in yellow or higher.  Should match
#'   the value passed to \code{\link{ras_validate}}.  Default \code{2.5}.
#' @param ... Currently unused.
#'
#' @return Invisibly returns \code{x}.  Called for its side effect of writing
#'   plot files or rendering to the active graphics device.
#'
#' @seealso \code{\link{ras}} for the function that produces the
#'   \code{"ras"} object. \code{\link{print.ras}} for the console summary.
#'
#' @examples
#' ## Build a minimal "ras" object by hand (see ras() for the full pipeline)
#' set.seed(7)
#' xg <- seq(1, 2000, by = 10)
#' yg <- c(seq(0, 9, length.out = length(xg) %/% 2),
#'         seq(9, 1, length.out = length(xg) - length(xg) %/% 2)) +
#'       rnorm(length(xg), sd = 0.4)
#' detection <- list(
#'   tau_hats         = xg[100],
#'   all.changepoints = xg[c(95, 100, 105)],
#'   all.p.values     = c(5, 12, 4),
#'   left.slopes      = 0.3,
#'   right.slopes     = -0.3
#' )
#' result <- structure(
#'   list(scan = list(x = xg, y = yg), detection = detection,
#'        chrom = 1, save_dir = tempdir()),
#'   class = "ras")
#'
#' plot(result, device = "screen")
#' plot(result, zoom = TRUE, device = "screen")
#' @export
plot.ras <- function(x, zoom = FALSE, device = "pdf", p.threshold = 8,
                     y_cap = NULL, min_display_p = 1, xlim = NULL,
                     zoom_half_width = 3000, ncol = 3, min_signal = 2.5, ...) {
  if (!zoom) {
    plot_ras_scan(
      x                = x$scan$x,
      y                = x$scan$y,
      detection.result = x$detection,
      this_chrom       = x$chrom,
      save.directory   = x$save_dir,
      p.threshold      = p.threshold,
      device           = device,
      y_cap            = y_cap,
      min_display_p    = min_display_p,
      xlim             = xlim,
      min_signal       = min_signal
    )
  } else {
    plot_ras_zoom_regions(
      x                = x$scan$x,
      y                = x$scan$y,
      detection.result = x$detection,
      this_chrom       = x$chrom,
      save.directory   = x$save_dir,
      p.threshold      = p.threshold,
      device           = device,
      zoom_half_width  = zoom_half_width,
      ncol             = ncol,
      min_signal       = min_signal
    )
  }
  invisible(x)
}


#' Print a RAS Result Object
#'
#' Prints a one-line summary of a \code{"ras"} object showing the chromosome
#' and detected changepoint positions.
#'
#' @param x Object of class \code{"ras"} as returned by \code{\link{ras}}.
#' @param ... Currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{ras}}, \code{\link{plot.ras}}.
#' @export
print.ras <- function(x, ...) {
  n <- length(x$detection$tau_hats)
  cat(sprintf(
    "RAS result  chr %d  |  %s\n",
    x$chrom,
    if (n == 0) "no changepoints detected"
    else sprintf("%d changepoint(s) at: %s",
                 n, paste(x$detection$tau_hats, collapse = ", "))
  ))
  invisible(x)
}


#' Plot RAS Scan Profile with Detected Changepoints
#'
#' Produces two diagnostic plots for one chromosome: the RAS scan profile
#' with detected changepoints marked, and a dual-axis overlay of the profile
#' with Davies test significance at each candidate position.
#'
#' @param x Numeric vector. SNP position index grid from the scan
#'   (\code{scan$x} returned by \code{\link{ras_scan}}).
#' @param y Numeric vector. Averaged \eqn{-\log_{10}(p)}-value profile
#'   (\code{scan$y} returned by \code{\link{ras_scan}}).
#' @param detection.result List. Output from \code{\link{ras_validate}},
#'   containing \code{tau_hats}, \code{all.changepoints}, \code{all.p.values},
#'   \code{left.slopes}, and \code{right.slopes}.
#' @param this_chrom Integer. Chromosome number used in output file names and
#'   panel titles.
#' @param save.directory Character. Directory path for saved files.  Ignored
#'   when \code{device = "screen"}.
#' @param p.threshold Numeric. Significance reference line drawn on both plots
#'   (on the \eqn{-\log_{10}} scale).  Default \code{8}.
#' @param device Character. Output device: \code{"pdf"} (default),
#'   \code{"png"}, or \code{"screen"} (renders to the active graphics device
#'   without writing files).
#' @param y_cap Numeric or \code{NULL}. If provided, the y-axis is capped at
#'   this value; positions above the cap are annotated with upward arrows and
#'   their true values.  Useful when boundary effects produce extreme outliers.
#'   Default \code{NULL} (no cap).
#' @param min_display_p Numeric. Minimum \code{all.p.values} required for a
#'   candidate changepoint to appear in Plot 2.  Lower-valued candidates are
#'   omitted to reduce overplotting.  Default \code{1}.
#' @param xlim Numeric vector of length 2, or \code{NULL}. If provided, both
#'   plots show only the genomic range \code{[xlim[1], xlim[2]]}.
#'   Default \code{NULL} (full range).
#' @param min_signal Numeric. Low-signal colour boundary on the
#'   \eqn{-\log_{10}} scale.  Values below this are drawn in blue; values at
#'   or above it switch to yellow (and higher tiers).  A grey reference line
#'   is drawn at this level.  Default \code{2.5}.
#'
#' @details
#' \strong{Plot 1 — Scan profile.}  The scan line is drawn segment-by-segment
#' with colour determined by the local \eqn{-\log_{10}(p)} value:
#' \itemize{
#'   \item Red (\code{#d73027}): \eqn{\ge} \strong{p.threshold}
#'   \item Orange (\code{#fc8d59}): \eqn{\ge} \code{p.threshold / 2}
#'   \item Yellow (\code{#fee090}): \eqn{\ge} \strong{min_signal}
#'   \item Blue (\code{#91bfdb}): below \strong{min_signal}
#' }
#' Detected changepoints (\code{tau_hats}) are marked with vertical dashed
#' red lines, shaded bands, and labelled with their position and Davies
#' \eqn{-\log_{10}(p)}.  Left and right slope lines (dark green and purple)
#' are overlaid at each changepoint; their half-width is fixed at 2.5\% of
#' the total x-range.  When \strong{y_cap} is set, positions above the cap are
#' shown as arrows with their true value annotated.
#'
#' \strong{Plot 2 — Dual-axis overlay.}  The scan profile (left y-axis, blue
#' line) is overlaid with triangle markers showing the Davies
#' \eqn{-\log_{10}(p)} for each candidate (right y-axis, coloured by whether
#' it exceeds \strong{p.threshold}).  Vertical dashed segments connect markers
#' to the baseline.
#'
#' Output files are named:
#' \code{chr-<chrom>-cp-plot.<ext>} and
#' \code{chr-<chrom>-cp-p-values-plot.<ext>}.
#'
#' @return
#' Invisibly returns \code{NULL}.  Called for its side effect of writing plot
#' files or rendering to the active graphics device.
#'
#' @seealso
#' \code{\link{plot_ras_zoom_regions}} for zoomed panels around each
#' changepoint.
#' \code{\link{ras}} which calls this function automatically.
#' \code{\link{ras_validate}} whose output is passed as
#' \code{detection.result}.
#'
#' @keywords internal
plot_ras_scan <- function(x, y, detection.result, this_chrom, save.directory,
                          p.threshold = 8, device = "pdf", y_cap = NULL,
                          min_display_p = 1, xlim = NULL, min_signal = 2.5) {

  # Restore the caller's graphics parameters on exit. Only device = "screen"
  # draws on the caller's active device; the pdf/png branches open and close
  # their own device, whose par() settings are discarded when it closes.
  if (device == "screen") {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  tau_hats         <- detection.result$tau_hats
  left.slopes      <- detection.result$left.slopes
  right.slopes     <- detection.result$right.slopes
  all.changepoints <- detection.result$all.changepoints
  all.p.values     <- detection.result$all.p.values

  # apply xlim subsetting if requested
  if (!is.null(xlim)) {
    keep      <- x >= xlim[1] & x <= xlim[2]
    x         <- x[keep]
    y         <- y[keep]
    tau_hats  <- tau_hats[tau_hats >= xlim[1] & tau_hats <= xlim[2]]
    cp_keep   <- all.changepoints >= xlim[1] & all.changepoints <= xlim[2]
    all.changepoints <- all.changepoints[cp_keep]
    all.p.values     <- all.p.values[cp_keep]
    # realign slopes to surviving tau_hats
    orig_tau   <- detection.result$tau_hats
    keep_k     <- which(orig_tau >= xlim[1] & orig_tau <= xlim[2])
    left.slopes  <- left.slopes[keep_k]
    right.slopes <- right.slopes[keep_k]
  }

  ext <- switch(device, png = ".png", pdf = ".pdf", screen = "")
  open_device <- function(stem) {
    if (device == "screen") return(invisible(NULL))
    path <- paste0(save.directory, "/chr-", this_chrom, stem, ext)
    if (device == "png") png(filename = path, width = 12 * 200, height = 6 * 200, res = 200)
    else                 pdf(file    = path, width = 10,  height = 5.5)
  }
  close_device <- function() {
    if (device != "screen") dev.off()
  }

  # ── colour scheme ─────────────────────────────────────────────────────────
  col_high   <- "#d73027"   # >= p.threshold
  col_mid    <- "#fc8d59"   # >= p.threshold/2
  col_low    <- "#fee090"   # >= 2.5
  col_bg     <- "#91bfdb"   # below 2.5
  col_scan   <- "#2c7bb6"
  col_shade  <- rgb(0.95, 0.80, 0.80, 0.35)

  sig_color <- function(v) {
    ifelse(v >= p.threshold,       col_high,
    ifelse(v >= p.threshold / 2,   col_mid,
    ifelse(v >= min_signal,        col_low,  col_bg)))
  }

  slope_hw <- diff(range(x)) * 0.025   # half-width for slope line segments

  # apply y_cap: clip display range, remember which positions are above cap
  y_display <- y
  clipped   <- logical(length(y))
  if (!is.null(y_cap)) {
    clipped             <- y > y_cap
    y_display[clipped]  <- y_cap
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Plot 1: scan profile
  # ══════════════════════════════════════════════════════════════════════════
  open_device("-cp-plot")
  par(mar = c(8.5, 4, 4, 2) + 0.1)

  y_max  <- max(y_display, na.rm = TRUE)
  y_ceil <- y_max * 1.18

  plot(x, y_display, type = "n",
       main  = paste0("Chr ", this_chrom,
                      "  RAS Scan  (p-threshold = 1e-",
                      round(p.threshold), ")"),
       xlab  = "SNP Index",
       ylab  = expression(-log[10](p)),
       ylim  = c(0, y_ceil))

  # significance reference lines
  abline(h = p.threshold, col = col_high, lty = 2, lwd = 0.9)
  abline(h = min_signal,  col = "grey55", lty = 3, lwd = 0.8)

  # shaded region around each detected changepoint
  if (length(tau_hats) > 0) {
    for (tau in tau_hats) {
      rect(tau - slope_hw, 0, tau + slope_hw, y_ceil,
           col = col_shade, border = NA)
    }
  }

  # colour-coded scan line (using clipped values)
  seg_cols <- sig_color(y_display)
  for (i in seq_len(length(x) - 1L)) {
    lines(x[i:(i + 1L)], y_display[i:(i + 1L)], col = seg_cols[i], lwd = 1.3)
  }

  # arrow + label for any positions clipped above y_cap
  if (any(clipped)) {
    clip_x <- x[clipped]
    clip_y <- y[clipped]
    arrows(clip_x, y_cap * 0.92, clip_x, y_cap * 0.99,
           length = 0.07, col = col_high, lwd = 1.5)
    text(clip_x, y_cap * 0.88,
         sprintf("%.1f", clip_y),
         col = col_high, cex = 0.65, font = 2)
  }

  # left / right slope lines at each changepoint
  if (length(tau_hats) > 0 && !is.null(left.slopes) && length(left.slopes) > 0) {
    for (k in seq_along(tau_hats)) {
      tau <- tau_hats[k]
      y0  <- y_display[which.min(abs(x - tau))]
      lines(c(tau - slope_hw, tau),
            c(y0 - left.slopes[k]  * slope_hw, y0),
            col = "darkgreen", lwd = 2.2)
      lines(c(tau, tau + slope_hw),
            c(y0, y0 + right.slopes[k] * slope_hw),
            col = "purple",    lwd = 2.2)
    }
  }

  # vertical changepoint lines
  abline(v = tau_hats, col = "red", lty = 2, lwd = 1.3)

  # text labels at changepoints — stagger vertically to avoid overlap
  if (length(tau_hats) > 0) {
    # compute stagger levels: group changepoints that are too close on x-axis
    min_x_gap  <- diff(range(x)) * 0.06
    stagger_y  <- numeric(length(tau_hats))
    stagger_levels <- c(0.97, 0.87, 0.77)
    last_x     <- -Inf
    level_idx  <- 1L
    for (k in seq_along(tau_hats)) {
      tau <- tau_hats[k]
      if ((tau - last_x) < min_x_gap) {
        level_idx <- (level_idx %% length(stagger_levels)) + 1L
      } else {
        level_idx <- 1L
      }
      stagger_y[k] <- stagger_levels[level_idx]
      last_x       <- tau
    }
    for (k in seq_along(tau_hats)) {
      tau   <- tau_hats[k]
      cp_pv <- if (length(all.changepoints) > 0) {
        all.p.values[which.min(abs(all.changepoints - tau))]
      } else NA
      lbl <- if (!is.na(cp_pv) && cp_pv > 0) {
        sprintf("pos %d\n1e-%.1f", tau, cp_pv)
      } else {
        sprintf("pos %d", tau)
      }
      text(tau, y_ceil * stagger_y[k], lbl,
           col = "red", cex = 0.65, adj = c(0.5, 1), font = 2)
    }
  }

  legend("bottomleft", inset = c(0, -0.45), xpd = TRUE,
         bty = "o", bg = "white", box.col = "grey80",
         cex = 0.72,
         legend = c(
           sprintf("-log10(p) >= %.0f",        p.threshold),
           sprintf("-log10(p) >= %.0f",        p.threshold / 2),
           sprintf("-log10(p) >= %.1f",        min_signal),
           sprintf("below %.1f",                   min_signal),
           "Detected changepoint",
           "Left slope",
           "Right slope"
         ),
         col  = c(col_high, col_mid, col_low, col_bg, "red", "darkgreen", "purple"),
         lty  = c(1, 1, 1, 1, 2, 1, 1),
         lwd  = c(2, 2, 2, 2, 1.3, 2.2, 2.2))

  close_device()

  # ══════════════════════════════════════════════════════════════════════════
  # Plot 2: dual-axis overlay — RAS profile + changepoint significance
  # ══════════════════════════════════════════════════════════════════════════
  open_device("-cp-p-values-plot")
  par(mar = c(8, 4, 4, 5) + 0.1)

  # use clipped y for the scan profile in Plot 2 as well
  y_range <- c(0, max(y_display, na.rm = TRUE) * 1.2)
  # only show candidates above min_display_p to avoid overplotting
  show_idx <- all.p.values >= min_display_p
  vis_cp   <- all.changepoints[show_idx]
  vis_pv   <- all.p.values[show_idx]

  p_range <- if (length(vis_pv) > 0 && max(vis_pv, na.rm = TRUE) > 0) {
    c(0, max(vis_pv, na.rm = TRUE) * 1.5)
  } else c(0, 1)

  plot(x, y_display, type = "l", col = col_scan, lwd = 1.6,
       main  = paste0("Chr ", this_chrom,
                      "  -  RAS Profile & Changepoint Significance"),
       xlim  = range(x),
       ylab  = expression(-log[10](p) ~ "(RAS scan)"),
       xlab  = "SNP Index",
       ylim  = y_range)

  abline(h = min_signal, col = "grey60", lty = 3, lwd = 0.8)
  text(x[1], min_signal, sprintf("min signal = %.1f", min_signal),
       adj = c(0, -0.4), col = "grey40", cex = 0.75)
  abline(v = tau_hats,  col = "red",    lty = 2, lwd = 1.2)

  # mark clipped positions with upward arrows
  if (any(clipped)) {
    arrows(x[clipped], y_cap * 0.88, x[clipped], y_cap * 0.98,
           length = 0.07, col = col_high, lwd = 1.5)
  }

  if (length(vis_cp) > 0) {
    pt_cols <- ifelse(vis_pv >= p.threshold, col_high, col_mid)

    par(new = TRUE)
    plot(vis_cp, vis_pv,
         col = pt_cols, pch = 17, cex = 0.85,
         ylab = "", xlab = "", xaxt = "n", yaxt = "n",
         ylim = p_range, xlim = range(x))

    segments(x0 = vis_cp, y0 = vis_pv,
             x1 = vis_cp, y1 = 0,
             col = pt_cols, lty = 2, lwd = 0.8)

    abline(h = p.threshold, lty = 2, col = col_high, lwd = 1)
    text(x[1], p.threshold, sprintf("threshold: p = 1e-%g", p.threshold),
         adj = c(0, -0.4), col = col_high, cex = 0.75)
  }

  axis(4)
  mtext(expression(-log[10](p) ~ "(changepoint test)"), side = 4, line = 3)

  legend("bottomright", inset = c(0, -0.35), xpd = TRUE,
         bty = "o", bg = "white", box.col = "grey80",
         cex = 0.72,
         legend = c(
           "RAS scan profile",
           sprintf("Candidate (p >= 1e-%.0f)", p.threshold),
           "Candidate (weaker)",
           "Detected changepoint",
           "Significance threshold"
         ),
         col  = c(col_scan, col_high, col_mid, "red", col_high),
         pch  = c(NA, 17, 17, NA, NA),
         lty  = c(1, NA, NA, 2, 2),
         lwd  = c(1.6, NA, NA, 1.2, 1))

  close_device()

  invisible(NULL)
}


#' Zoom-In Plots Around Each Detected Changepoint
#'
#' Generates a multi-panel figure with one panel per detected changepoint,
#' each showing a zoomed view of the RAS scan profile centred on that
#' position.
#'
#' @param x Numeric vector. SNP position index grid used in the scan.
#' @param y Numeric vector. Averaged \eqn{-\log_{10}(p)}-value profile.
#' @param detection.result List. Output from \code{\link{ras_validate}},
#'   the same object passed to \code{\link{plot_ras_scan}}.
#' @param this_chrom Integer. Chromosome number used in file names and panel
#'   titles.
#' @param save.directory Character. Directory for output files.  Ignored when
#'   \code{device = "screen"}.
#' @param p.threshold Numeric. Significance reference line drawn on each panel
#'   (on the \eqn{-\log_{10}} scale).  Default \code{8}.
#' @param device Character. \code{"pdf"} (default), \code{"png"}, or
#'   \code{"screen"}.
#' @param zoom_half_width Numeric. Half-width of each zoom window in the same
#'   units as \code{x}.  Default \code{3000}.
#' @param ncol Integer. Number of columns in the panel grid.  Default
#'   \code{3}.
#' @param min_signal Numeric. Low-signal colour boundary on the
#'   \eqn{-\log_{10}} scale; also drawn as a grey reference line on each
#'   panel.  Default \code{2.5}.
#'
#' @details
#' Each panel shows a \eqn{[\tau - \texttt{zoom\_half\_width},\;
#' \tau + \texttt{zoom\_half\_width}]} window around changepoint \eqn{\tau}.
#' The scan line is drawn in blue with colour-coded overlays for regions above
#' the significance thresholds (\strong{min_signal}, \code{p.threshold / 2},
#' \strong{p.threshold}).  Left and right slope lines are plotted in dark green
#' and purple respectively, with their horizontal extent scaled to span
#' approximately 45\% of the local y-range (minimum 15\% of
#' \code{zoom_half_width}).
#'
#' If \code{detection.result$tau_hats} is empty the function returns
#' invisibly with a message.  Blank panels are added to complete the grid if
#' the number of changepoints is not a multiple of \code{ncol}.
#'
#' The output file is named \code{chr-<chrom>-zoom.<ext>}.
#'
#' @return
#' Invisibly returns \code{NULL}.  Called for its side effect of writing a
#' multi-panel plot file or rendering to the active graphics device.
#'
#' @seealso
#' \code{\link{plot_ras_scan}} for the full-chromosome overview plots.
#' \code{\link{ras}} which calls this function automatically.
#'
#' @keywords internal
plot_ras_zoom_regions <- function(x, y, detection.result,
                                  this_chrom, save.directory,
                                  p.threshold    = 8,
                                  device         = "pdf",
                                  zoom_half_width = 3000,
                                  ncol           = 3,
                                  min_signal     = 2.5) {

  tau_hats         <- unique(detection.result$tau_hats)
  left.slopes      <- detection.result$left.slopes
  right.slopes     <- detection.result$right.slopes
  all.changepoints <- detection.result$all.changepoints
  all.p.values     <- detection.result$all.p.values

  n_tau <- length(tau_hats)
  if (n_tau == 0) { message("No changepoints to zoom into."); return(invisible(NULL)) }

  nrow_panels <- ceiling(n_tau / ncol)

  # Restore the caller's graphics parameters on exit. Only device = "screen"
  # draws on the caller's active device; the pdf/png branches open and close
  # their own device, whose par() settings are discarded when it closes.
  if (device == "screen") {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  # ── colour helpers ────────────────────────────────────────────────────────
  col_high  <- "#d73027"
  col_mid   <- "#fc8d59"
  col_low   <- "#fee090"
  col_bg    <- "#91bfdb"
  col_shade <- rgb(0.95, 0.80, 0.80, 0.40)
  sig_color <- function(v) {
    ifelse(v >= p.threshold,     col_high,
    ifelse(v >= p.threshold / 2, col_mid,
    ifelse(v >= min_signal,      col_low, col_bg)))
  }

  # ── open output device ────────────────────────────────────────────────────
  ext <- switch(device, png = ".png", pdf = ".pdf", screen = "")
  if (device != "screen") {
    path <- paste0(save.directory, "/chr-", this_chrom, "-zoom", ext)
    pw   <- ncol * 3.2;  ph <- nrow_panels * 2.8
    if (device == "png")
      png(filename = path, width = round(pw * 200), height = round(ph * 200), res = 200)
    else
      pdf(file = path, width = pw, height = ph)
  }

  par(mfrow = c(nrow_panels, ncol),
      mar   = c(3, 3, 2.2, 1) + 0.1,
      mgp   = c(1.8, 0.5, 0))

  for (k in seq_along(tau_hats)) {
    tau  <- tau_hats[k]
    x_lo <- max(tau - zoom_half_width, min(x))
    x_hi <- min(tau + zoom_half_width, max(x))

    # subset scan to window
    idx <- x >= x_lo & x <= x_hi
    xz  <- x[idx];  yz <- y[idx]
    if (length(xz) < 2) { plot.new(); next }

    y_local <- max(yz, na.rm = TRUE)
    y_ceil  <- y_local * 1.22

    # nearby candidates
    cp_in   <- all.changepoints >= x_lo & all.changepoints <= x_hi
    vis_cp  <- all.changepoints[cp_in]
    vis_pv  <- all.p.values[cp_in]

    # title p-value label
    nearest <- which.min(abs(all.changepoints - tau))
    cp_pv   <- all.p.values[nearest]
    title_lbl <- if (length(cp_pv) > 0 && cp_pv > 0)
      sprintf("Chr%d  pos %d  (1e-%.1f)", this_chrom, tau, cp_pv)
    else
      sprintf("Chr%d  pos %d", this_chrom, tau)

    # ── draw panel ──────────────────────────────────────────────────────────
    plot(xz, yz, type = "n",
         main = title_lbl, cex.main = 0.80,
         xlab = "Position", ylab = expression(-log[10](p)),
         xlim = c(x_lo, x_hi), ylim = c(0, y_ceil),
         cex.axis = 0.72, cex.lab = 0.78)

    abline(h = p.threshold, col = col_high, lty = 2, lwd = 0.8)
    abline(h = min_signal,  col = "grey55", lty = 3, lwd = 0.7)

    # shading around tau
    shade_hw <- zoom_half_width * 0.04
    rect(tau - shade_hw, 0, tau + shade_hw, y_ceil,
         col = col_shade, border = NA)

    # base scan line: smooth single colour (same blue as Plot 2)
    lines(xz, yz, col = "#2c7bb6", lwd = 1.4)

    # overlay highlighted sections above thresholds
    for (thr in c(min_signal, p.threshold / 2, p.threshold)) {
      over <- yz >= thr
      if (!any(over)) next
      col_over <- sig_color(rep(thr, 1))
      run_start <- which(diff(c(FALSE, over)) == 1)
      run_end   <- which(diff(c(over, FALSE)) == -1)
      for (r in seq_along(run_start)) {
        ri <- run_start[r]:run_end[r]
        lines(xz[ri], yz[ri], col = col_over, lwd = 1.8)
      }
    }

    # slope lines — hw scaled to ~45% of local y range, with a minimum of
    # 15% of the zoom window so steep-slope lines are always visible.
    # Lines that run off the y-axis are clipped automatically by R.
    if (!is.null(left.slopes) && k <= length(left.slopes) &&
        left.slopes[k] != 0 && right.slopes[k] != 0) {
      y0    <- yz[which.min(abs(xz - tau))]
      hw_min <- zoom_half_width * 0.15
      hw_l  <- max(min(abs(y_local * 0.45 / left.slopes[k]),  zoom_half_width * 0.9), hw_min)
      hw_r  <- max(min(abs(y_local * 0.45 / right.slopes[k]), zoom_half_width * 0.9), hw_min)
      lines(c(tau - hw_l, tau),
            c(y0 - left.slopes[k]  * hw_l, y0),
            col = "darkgreen", lwd = 2)
      lines(c(tau, tau + hw_r),
            c(y0, y0 + right.slopes[k] * hw_r),
            col = "purple",    lwd = 2)
    }

    # detected changepoint line
    abline(v = tau, col = "red", lty = 2, lwd = 1.2)

  }

  # blank panels to fill grid
  n_blank <- nrow_panels * ncol - n_tau
  for (i in seq_len(n_blank)) plot.new()

  if (device != "screen") dev.off()
  invisible(NULL)
}
