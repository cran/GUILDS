context("logLikelihood.Guilds")

test_that("logLikelihood.Guilds: use", {
  set.seed(666)
  J <- 1000
  theta <- 100
  alpha_x <- 0.1
  alpha_y <- alpha_x

  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)
  LL1 <- logLikelihood.Guilds(parameters = c(theta, alpha_x, alpha_y),
                              model = "D0",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)
  LL2 <- logLikelihood.Guilds(parameters = c(20, 0.5, 0.25),
                              model = "D0",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)

  LL3 <- logLikelihood.Guilds(parameters = c(theta, alpha_x, alpha_x * 0.1),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)

  LL4 <- logLikelihood.Guilds(parameters = c(theta, alpha_x, alpha_x),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)
  testthat::expect_gt(LL1[[1]], LL2[[1]])
  testthat::expect_gt(LL1[[1]], LL3[[1]])
  testthat::expect_equal(LL1[[1]], LL4[[1]])

  alpha_y <- 0.01

  simul_data <- generate.Guilds(theta, alpha_x, alpha_y, J)
  LL1 <- logLikelihood.Guilds(parameters = c(theta, alpha_x, alpha_y),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)
  LL2 <- logLikelihood.Guilds(parameters = c(20, 0.99, 0.1),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)

  LL3 <- logLikelihood.Guilds(parameters = c(theta, alpha_x, alpha_x),
                              model = "D0",
                              simul_data$guildX, simul_data$guildY,
                              verbose = FALSE)


  testthat::expect_gt(LL1[[1]], LL2[[1]])
  testthat::expect_gt(LL1[[1]], LL3[[1]])

  testthat::expect_output(
    LL1 <- logLikelihood.Guilds(parameters = c(theta, alpha_x, alpha_y),
                              model = "D1",
                              simul_data$guildX, simul_data$guildY,
                              verbose = TRUE)
)
})
