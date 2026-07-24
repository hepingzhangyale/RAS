make_screening_data <- function(N_sample = 40, N_snp = 50) {
  set.seed(42)
  geno    <- matrix(sample(0:2, N_sample * N_snp, replace = TRUE,
                           prob = c(0.6, 0.3, 0.1)),
                   nrow = N_sample, ncol = N_snp)
  weights <- rnorm(N_snp)
  leftout <- 1:N_sample
  pgs.mat <- compute_pgs_matrix(geno, leftout, weights)

  df <- data.frame(
    phenotype2  = rnorm(N_sample),
    sex         = sample(c("Male", "Female"), N_sample, replace = TRUE),
    age         = sample(20:60, N_sample, replace = TRUE)
  )
  df$age_squared <- df$age^2
  df$age_sex     <- df$age * as.numeric(df$sex == "Male")
  for (i in 1:10) df[[paste0("pc", i)]] <- rnorm(N_sample)

  list(geno = geno, pgs.mat = pgs.mat, df = df, N_snp = N_snp)
}

test_that("screen_forward_max_region returns numeric vector of correct length", {
  d      <- make_screening_data()
  result <- screen_forward_max_region(
    geno              = d$geno,
    pgs.mat           = d$pgs.mat,
    this.df           = d$df,
    num_signals       = 0,
    isPlot            = FALSE,
    skip1             = 10,
    skip2             = 5,
    min_window_size   = 5,
    max_window_size   = 20,
    is_continuous     = TRUE,
    covariate_formula = "sex + age"
  )
  expect_true(is.numeric(result))
  expected_len <- length(seq(1, d$N_snp, by = 10))
  expect_equal(length(result), expected_len)
})

test_that("binary score test agrees with glm and returns valid output", {
  d <- make_screening_data()
  d$df$phenotype2 <- rbinom(nrow(d$df), 1, 0.5)  # binary outcome

  args <- list(
    geno = d$geno, pgs.mat = d$pgs.mat, this.df = d$df, num_signals = 0,
    isPlot = FALSE, skip1 = 10, skip2 = 5,
    min_window_size = 5, max_window_size = 20,
    is_continuous = FALSE, covariate_formula = "sex + age"
  )
  p_glm   <- do.call(screen_forward_max_region, c(args, scan_test = "glm"))
  p_score <- do.call(screen_forward_max_region, c(args, scan_test = "score"))

  expect_equal(length(p_score), length(seq(1, d$N_snp, by = 10)))
  expect_true(all(p_score >= 0 | is.nan(p_score) | is.infinite(p_score)))
  # score is asymptotically equivalent to the glm Wald test: profiles correlate
  ok <- is.finite(p_glm) & is.finite(p_score)
  if (sum(ok) >= 3) {
    expect_gt(cor(p_glm[ok], p_score[ok]), 0.9)
  }
})

test_that("incomplete covariate rows are dropped, matching the same data pre-dropped", {
  # Regression guard. The fast paths once resolved completeness on the OUTPUT of
  # model.matrix(), which has already applied na.action and dropped the NA rows.
  # `keep` then no longer indexed this.df, so it desynchronised from the
  # full-length this.pgs and the run died with "(subscript) logical subscript
  # too long". Every test above uses complete data, which is why it survived to
  # a release. Blanking a covariate must give exactly the profile obtained by
  # deleting those rows up front -- that is what "removes incomplete cases"
  # means, and it must hold on all three scan paths.
  d       <- make_screening_data()
  na_rows <- c(3, 11, 27, 34)

  common <- list(
    num_signals = 0, isPlot = FALSE, skip1 = 10, skip2 = 5,
    min_window_size = 5, max_window_size = 20,
    covariate_formula = "sex + age + pc1"
  )

  expect_same <- function(df, is_cont, test) {
    df_na <- df
    df_na$pc1[na_rows] <- NA

    with_na <- do.call(screen_forward_max_region, c(list(
      geno = d$geno, pgs.mat = d$pgs.mat, this.df = df_na,
      is_continuous = is_cont, scan_test = test), common))

    pre_dropped <- do.call(screen_forward_max_region, c(list(
      geno    = d$geno[-na_rows, , drop = FALSE],
      pgs.mat = d$pgs.mat[-na_rows, , drop = FALSE],
      this.df = df[-na_rows, , drop = FALSE],
      is_continuous = is_cont, scan_test = test), common))

    expect_equal(with_na, pre_dropped)
  }

  expect_same(d$df, TRUE, "glm")        # continuous -> exact-FWL path

  df_bin <- d$df
  df_bin$phenotype2 <- rbinom(nrow(df_bin), 1, 0.5)
  expect_same(df_bin, FALSE, "score")   # binary -> Rao score path
  expect_same(df_bin, FALSE, "glm")     # binary -> legacy glm path
})

test_that("screen_forward_max_region returns non-negative values", {
  d      <- make_screening_data()
  result <- screen_forward_max_region(
    geno              = d$geno,
    pgs.mat           = d$pgs.mat,
    this.df           = d$df,
    num_signals       = 0,
    isPlot            = FALSE,
    skip1             = 10,
    skip2             = 5,
    min_window_size   = 5,
    max_window_size   = 20,
    is_continuous     = TRUE,
    covariate_formula = "sex + age"
  )
  expect_true(all(result >= 0 | is.nan(result) | is.infinite(result)))
})
