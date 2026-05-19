test_that("release_memory() runs without error", {
  expect_no_error(release_memory(verbose = FALSE))
})

test_that("release_memory() returns an integer (or NA_integer_)", {
  out <- release_memory(verbose = FALSE)
  expect_true(is.integer(out) || is.na(out))
})

test_that("release_memory() returns NA on non-Linux platforms", {
  skip_if(identical(.Platform$OS.type, "unix") && grepl("linux", R.version$os),
          "skipped on Linux — returns 0 or 1 there")
  out <- release_memory(verbose = FALSE)
  expect_true(is.na(out))
})

test_that("release_memory() returns 0 or 1 on Linux", {
  skip_if_not(identical(.Platform$OS.type, "unix") && grepl("linux", R.version$os),
              "only runs on Linux")
  out <- release_memory(verbose = FALSE)
  expect_true(out %in% c(0L, 1L))
})

test_that("release_memory() emits a message when verbose = TRUE", {
  expect_message(release_memory(verbose = TRUE), "malloc_trim returned")
})

test_that("release_memory() emits no message when verbose = FALSE", {
  expect_no_message(release_memory(verbose = FALSE))
})

test_that("release_memory() return value is invisible", {
  # withVisible() reports whether the return was auto-printed
  rv <- withVisible(release_memory(verbose = FALSE))
  expect_false(rv$visible)
})
