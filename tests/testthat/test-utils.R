test_that("get_local_maximum returns index of peak within window", {
  y <- c(1, 2, 5, 3, 1)
  expect_equal(get_local_maximum(y, x0 = 3, window.size = 2), 3)
})

test_that("get_local_maximum respects window boundary", {
  y <- c(1, 2, 5, 3, 8)
  # window around x0=3 with size=1 covers indices 2:4, max is index 3
  expect_equal(get_local_maximum(y, x0 = 3, window.size = 1), 3)
})

test_that("get_local_maximum handles edge: x0 near start", {
  y <- c(9, 2, 3, 4, 5)
  expect_equal(get_local_maximum(y, x0 = 1, window.size = 2), 1)
})

test_that("get_local_maximum handles edge: x0 near end", {
  y <- c(1, 2, 3, 4, 9)
  expect_equal(get_local_maximum(y, x0 = 5, window.size = 2), 5)
})
