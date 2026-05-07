make_test_data <- function(N_sample = 60, N_snp = 10) {
  set.seed(42)
  geno <- matrix(sample(0:2, N_sample * N_snp, replace = TRUE,
                        prob = c(0.6, 0.3, 0.1)),
                 nrow = N_sample, ncol = N_snp)
  df <- data.frame(
    sex         = sample(c("Male", "Female"), N_sample, replace = TRUE),
    age         = sample(20:60, N_sample, replace = TRUE)
  )
  df$age_squared <- df$age^2
  df$age_sex     <- df$age * as.numeric(df$sex == "Male")
  for (i in 1:10) df[[paste0("pc", i)]] <- rnorm(N_sample)
  list(geno = geno, df = df, N_sample = N_sample, N_snp = N_snp)
}

test_that("compute_gwas_weights returns matrix with correct dimensions", {
  d <- make_test_data()
  phenotype1 <- rnorm(d$N_sample)
  result <- compute_gwas_weights(d$geno, phenotype1, 1:d$N_sample,
                                 d$df, is_continuous = TRUE)
  expect_equal(dim(result), c(d$N_snp, 4))
  expect_equal(colnames(result), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
})

test_that("compute_gwas_weights works for binary trait", {
  d <- make_test_data()
  phenotype1 <- rbinom(d$N_sample, 1, 0.3)
  covariate_formula <- "this.x + sex + age"
  result <- compute_gwas_weights(d$geno, phenotype1, 1:d$N_sample,
                                 d$df, is_continuous = FALSE,
                                 covariate_formula = covariate_formula)
  expect_equal(dim(result), c(d$N_snp, 4))
})

test_that("compute_pgs_matrix returns correct dimensions", {
  set.seed(42)
  geno    <- matrix(sample(0:2, 50 * 20, replace = TRUE), nrow = 50, ncol = 20)
  weights <- rnorm(20)
  leftout <- 1:15
  result  <- compute_pgs_matrix(geno, leftout, weights)
  expect_equal(dim(result), c(15, 20))
})

test_that("compute_pgs_matrix imputes NA to 0", {
  geno        <- matrix(c(1, NA, 2, 1), nrow = 2, ncol = 2)
  weights     <- c(1, 1)
  result      <- compute_pgs_matrix(geno, 1:2, weights)
  expect_false(any(is.na(result)))
})

test_that("compute_pgs_matrix values equal genotype times weight", {
  geno    <- matrix(c(2, 0, 1, 2), nrow = 2, ncol = 2)
  weights <- c(3, 2)
  result  <- compute_pgs_matrix(geno, 1:2, weights)
  expect_equal(result[1, ], c(6, 2))
  expect_equal(result[2, ], c(0, 4))
})
