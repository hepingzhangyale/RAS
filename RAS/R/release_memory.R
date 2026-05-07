#' Release Memory Back to the OS
#'
#' Runs full garbage collection and, on Linux/glibc, calls
#' \code{malloc_trim(0)} to return free heap pages to the operating system.
#' Useful after large temporary matrices are removed in memory-heavy RAS
#' pipeline stages.
#'
#' On non-Linux platforms the C call is skipped and \code{NA} is returned
#' silently; no error is raised.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints a one-line message
#'   with the \code{malloc_trim} return value so RSS changes can be monitored
#'   in pipeline logs.
#'
#' @return Invisibly returns the \code{malloc_trim(0)} result:
#'   \code{1} if heap pages were returned to the OS,
#'   \code{0} if nothing was returned,
#'   \code{NA_integer_} on non-Linux platforms.
#' @export
release_memory <- function(verbose = TRUE) {
  invisible(gc(full = TRUE))

  out <- tryCatch(
    .Call("RAS_malloc_trim", PACKAGE = "RAS"),
    error = function(e) NA_integer_
  )

  if (isTRUE(verbose)) {
    message("gc(full = TRUE) complete; malloc_trim returned: ", out)
  }

  invisible(out)
}
