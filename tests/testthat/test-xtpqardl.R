make_panel <- function() {
  set.seed(3)
  do.call(rbind, lapply(1:8, function(i) {
    T <- 40
    x <- cumsum(rnorm(T))
    y <- numeric(T)
    for (t in 2:T) y[t] <- y[t - 1] - 0.4 * (y[t - 1] - 0.8 * x[t - 1]) +
      0.3 * (x[t] - x[t - 1]) + rnorm(1, 0, 0.3)
    data.frame(id = i, t = 1:T, y = y, x = x,
               d_y = c(NA, diff(y)), d_x = c(NA, diff(x)),
               L_y = c(NA, y[-T]), L_x = c(NA, x[-T]))
  }))
}

test_that("MG equals the average of the panel estimates", {
  d <- make_panel()
  f <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5,
                model = "mg")
  expect_equal(as.numeric(f$rho_mg), mean(f$rho_all[, 1]))
  expect_equal(as.numeric(f$beta_mg), mean(f$beta_all[, 1]))
  expect_equal(f$n_obs, 8 * 39)
})

test_that("PMG pools the long run and differs from MG", {
  d <- make_panel()
  mg  <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "mg")
  pmg <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "pmg")
  expect_false(isTRUE(all.equal(as.numeric(mg$rho_mg), as.numeric(pmg$rho_mg))))
  expect_true(abs(pmg$beta_mg[1] - 0.8) < 0.2)
  expect_true(is.list(pmg$hausman))
})

test_that("half-life uses ln(0.5)/ln(1 + rho)", {
  d <- make_panel()
  f <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "mg")
  expect_equal(as.numeric(f$halflife_mg), log(0.5) / log(1 + as.numeric(f$rho_mg)))
})
