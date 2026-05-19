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
