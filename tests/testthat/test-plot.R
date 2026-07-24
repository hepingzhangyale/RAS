test_that("plot_ras_scan creates two PDF files", {
  set.seed(42)
  x    <- 1:50
  y    <- c(seq(0, 3, length.out = 25), seq(3, 0, length.out = 25))
  mock_result <- list(
    tau_hats         = 25L,
    all.changepoints = c(20L, 25L, 30L),
    all.p.values     = c(5, 10, 4)
  )

  tmp_dir <- tempdir()
  plot_ras_scan(x, y, mock_result, this_chrom = 1, save.directory = tmp_dir)

  expect_true(file.exists(file.path(tmp_dir, "chr-1-cp-plot.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "chr-1-cp-p-values-plot.pdf")))
})
