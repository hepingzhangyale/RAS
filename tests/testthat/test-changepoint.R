make_peak_series <- function(n = 200, peak_pos = 100, sd = 0.3) {
  set.seed(42)
  y <- c(seq(0, 5, length.out = peak_pos),
         seq(5, 0, length.out = n - peak_pos)) + rnorm(n, sd = sd)
  y
}

# --- slope_test ---

test_that("slope_test returns value between 0 and 1", {
  set.seed(1)
  x <- 1:30
  y <- 2 * x + rnorm(30, sd = 0.5)
  p <- slope_test(x, y, lower.tail = FALSE)
  expect_true(p >= 0 && p <= 1)
})

test_that("slope_test: positive slope gives small p when lower.tail=FALSE", {
  # pt(t_stat, df, lower.tail=FALSE) = P(T > t_stat)
  # strong positive slope -> large positive t_stat -> P(T > large) ≈ 0
  set.seed(1)
  x <- 1:50
  y <- 5 * x + rnorm(50, sd = 0.5)
  p <- slope_test(x, y, lower.tail = FALSE)
  expect_lt(p, 0.05)
})

test_that("slope_test: negative slope gives small p when lower.tail=TRUE", {
  # pt(t_stat, df, lower.tail=TRUE) = P(T <= t_stat)
  # strong negative slope -> large negative t_stat -> P(T <= large_negative) ≈ 0
  set.seed(1)
  x <- 1:50
  y <- -5 * x + rnorm(50, sd = 0.5)
  p <- slope_test(x, y, lower.tail = TRUE)
  expect_lt(p, 0.05)
})

# --- get_break_points ---

test_that("get_break_points returns list with required elements", {
  set.seed(42)
  x      <- 1:100
  y      <- make_peak_series(100, peak_pos = 50, sd = 0.2)
  result <- get_break_points(x, y, t = 100)
  expect_true(is.list(result))
  expect_true(all(c("break.points", "p.values", "slope.left", "slope.right") %in% names(result)))
})

test_that("get_break_points p.value is between 0 and 1", {
  set.seed(42)
  x      <- 1:100
  y      <- make_peak_series(100, peak_pos = 50, sd = 0.2)
  result <- get_break_points(x, y, t = 100)
  expect_true(result$p.values >= 0 && result$p.values <= 1)
})

test_that("get_break_points returns p.value=1 on flat data", {
  x      <- 1:30
  y      <- rep(1, 30)
  result <- get_break_points(x, y, t = 30)
  expect_equal(result$p.values, 1)
})

# --- ras_detect ---

test_that("ras_detect returns list with required elements", {
  set.seed(42)
  x      <- 1:200
  y      <- make_peak_series(200, peak_pos = 100)
  result <- ras_detect(x, y, window_size = 100, skip = 10,
                       slope.p.values.threshold.left  = 1e-3,
                       slope.p.values.threshold.right = 1e-3)
  expected_names <- c("tau_hats", "p.values", "slope.left", "slope.right",
                      "all.changepoints", "all.p.values", "slope.angle", "previous_tau_hats")
  expect_true(all(expected_names %in% names(result)))
})

test_that("ras_detect output vectors are consistent length", {
  set.seed(42)
  x      <- 1:200
  y      <- make_peak_series(200, peak_pos = 100)
  result <- ras_detect(x, y, window_size = 100, skip = 10,
                       slope.p.values.threshold.left  = 1e-3,
                       slope.p.values.threshold.right = 1e-3)
  expect_equal(length(result$tau_hats), length(result$p.values))
  expect_equal(length(result$tau_hats), length(result$slope.left))
  expect_equal(length(result$tau_hats), length(result$slope.right))
})

# --- ras_validate ---

test_that("ras_validate returns list with required elements", {
  set.seed(42)
  n <- 300
  x <- 1:n
  y <- c(seq(0, 20, length.out = 150), seq(20, 0, length.out = 150)) + rnorm(n, sd = 0.05)

  # construct mock first-pass result with a known changepoint at position 150
  mock_first <- list(
    tau_hats         = 150L,
    p.values         = 0.001,
    slope.left       = 0.13,
    slope.right      = -0.13,
    all.changepoints = c(148L, 149L, 150L, 151L, 152L),
    all.p.values     = c(0.01, 0.005, 0.001, 0.005, 0.01),
    slope.angle      = 120,
    previous_tau_hats = 150L
  )

  result <- ras_validate(mock_first, x = x, y = y,
                         second_window_size = 30, p.value.threshold = 0.5)
  expected_names <- c("all.changepoints", "tau_hats", "all.p.values",
                      "left.slopes", "right.slopes")
  expect_true(all(expected_names %in% names(result)))
})
