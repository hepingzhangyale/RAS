# RAS 1.0.3

Resubmission addressing the CRAN reviewer's comments on the 1.0.0 submission,
with a bug fix folded in.

* Removed the `\dontrun{}` examples that called the unexported internal
  functions `plot_ras_scan()` and `plot_ras_zoom_regions()` through `:::`.
  These functions are internal (`@keywords internal`) and are already exercised
  by the runnable example of the exported `plot()` method for `"ras"` objects.
* Made the `ras_memory()` "abort" example runnable by wrapping the intentional
  `stop()` in `try()`, instead of hiding it in `\dontrun{}`.
* The plotting functions now capture the caller's graphics parameters and
  restore them via an immediate `on.exit(par(oldpar))`, so plotting a `"ras"`
  object no longer leaves the user's `par()` settings modified.
* Fixed a crash in `screen_forward_max_region()` (and therefore `ras()`) when
  covariates contained missing values. Incomplete cases are now dropped before
  the model matrix is built, matching the documented "removes incomplete cases"
  behaviour, on all three scan paths (continuous exact-FWL, binary score, and
  the legacy glm path).
