## Regression tests for 1.0.4

make_panel_k <- function(N = 6, TT = 17, K = 3, seed = 11) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) {
    X <- apply(matrix(rnorm(TT * K), TT, K), 2, cumsum)
    y <- numeric(TT)
    for (t in 2:TT) y[t] <- y[t - 1] - 0.4 * (y[t - 1] - sum(0.5 * X[t - 1, ])) +
      sum(0.2 * (X[t, ] - X[t - 1, ])) + rnorm(1, 0, 0.3)
    out <- data.frame(id = i, t = seq_len(TT), d_y = c(NA, diff(y)), L_y = c(NA, y[-TT]))
    for (k in seq_len(K)) {
      out[[paste0("d_x", k)]] <- c(NA, diff(X[, k]))
      out[[paste0("L_x", k)]] <- c(NA, X[-TT, k])
    }
    out
  }))
}

test_that("wald_test(joint = FALSE) prints for DFE instead of erroring", {
  d <- make_panel_k(K = 2)
  fit <- xtpqardl(d_y ~ d_x1 + d_x2, d, "id", "t", lr = c("L_y", "L_x1", "L_x2"),
                  tau = c(.25, .5, .75), model = "dfe")
  w <- wald_test(fit, joint = FALSE)
  expect_length(w$beta_test, 2)
  expect_output(print(w), "Long-Run Coefficients")
  wj <- wald_test(fit, joint = TRUE)
  expect_output(print(wj), "ECT Speed")
})

test_that("DFE has cross-quantile covariances and Koenker-Bassett tests", {
  d <- make_panel_k(K = 2)
  ## 1.0.5: the anova reproduction needs the "nid" estimator and all = TRUE
  fit <- xtpqardl(d_y ~ d_x1 + d_x2, d, "id", "t", lr = c("L_y", "L_x1", "L_x2"),
                  tau = c(.25, .5, .75), model = "dfe", dfe_cov = "nid")
  expect_null(fit$cross_quantile_cov_error)
  expect_true(fit$cross_quantile_cov)
  expect_false(anyNA(fit$beta_V))
  expect_false(anyNA(fit$rho_V))
  expect_true(isSymmetric(unname(fit$rho_V), tol = 1e-8))
  w <- wald_test(fit, all = TRUE)
  expect_true(is.finite(w$beta_test$statistic))
  expect_true(is.finite(w$rho_test$statistic))
  expect_false(is.null(w$kb_test))
  expect_true(is.finite(w$kb_test$joint[1, "pvalue"]))
  ## kb_test with all = TRUE reproduces quantreg::anova on the equivalent pooled regression
  dd <- d[!is.na(d$d_y), ]
  fq <- d_y ~ L_y + L_x1 + L_x2 + d_x1 + d_x2 + factor(id)
  fits <- lapply(c(.25, .5, .75), function(tt) quantreg::rq(fq, tau = tt, data = dd))
  ref <- suppressWarnings(anova(fits[[1]], fits[[2]], fits[[3]], joint = TRUE, se = "nid"))$table
  expect_equal(w$kb_test$joint[1, "Tn"], ref[1, "Tn"], tolerance = 1e-6)
  expect_equal(w$kb_test$joint[1, "pvalue"], ref[1, "pvalue"], tolerance = 1e-6)
})

test_that("PMG on very short panels: fallback covariance with warning, or informative error", {
  d <- make_panel_k(N = 6, TT = 17, K = 6)   # n - k = 2 per panel
  res <- tryCatch(
    withCallingHandlers(
      xtpqardl(d_y ~ d_x1 + d_x2 + d_x3 + d_x4 + d_x5 + d_x6, d, "id", "t",
               lr = c("L_y", paste0("L_x", 1:6)), tau = c(.25, .5, .75), model = "pmg"),
      ## 1.0.5: one warning names every estimator switch (summary.rq se = "nid"
      ## or the iid formula with a kernel sparsity estimate)
      warning = function(w) { expect_match(conditionMessage(w), "kernel sparsity|estimator used instead|dropped|Hausman"); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res), "Per-panel status")
    expect_match(conditionMessage(res), "pmg_weights")
  } else {
    expect_s3_class(res$panel_cov_status, "data.frame")
    expect_true(all(c("usable", "reason", "se") %in% names(res$panel_cov_status)))
    expect_true(all(is.finite(res$beta_mg)))
  }
  mg <- suppressWarnings(xtpqardl(d_y ~ d_x1 + d_x2 + d_x3 + d_x4 + d_x5 + d_x6, d, "id", "t",
                 lr = c("L_y", paste0("L_x", 1:6)), tau = c(.25, .5, .75), model = "mg"))
  expect_s3_class(mg$panel_cov_status, "data.frame")
  expect_equal(nrow(mg$panel_cov_status), 6 * 3)
  pe <- suppressWarnings(xtpqardl(d_y ~ d_x1 + d_x2 + d_x3 + d_x4 + d_x5 + d_x6, d, "id", "t",
                 lr = c("L_y", paste0("L_x", 1:6)), tau = c(.25, .5, .75), model = "pmg",
                 pmg_weights = "equal"))
  expect_equal(as.numeric(pe$beta_mg), as.numeric(mg$beta_mg), tolerance = 1e-10)
  expect_equal(pe$pmg_weights, "equal")
})

test_that("inverse-variance PMG is unchanged when covariances are estimable", {
  d <- make_panel_k(N = 6, TT = 40, K = 1, seed = 2)
  pmg <- xtpqardl(d_y ~ d_x1, d, "id", "t", lr = c("L_y", "L_x1"), tau = 0.5, model = "pmg")
  expect_true(all(pmg$panel_cov_status$usable))
  expect_true(all(pmg$panel_cov_status$se == "iid"))
  expect_true(is.list(pmg$hausman))
})

test_that("input checks on lr give clear errors", {
  d <- make_panel_k(K = 1)
  expect_error(xtpqardl(d_y ~ d_x1, d, "id", "t", lr = c("d_y", "L_x1"), model = "dfe"),
               "LAGGED LEVEL")
  expect_error(xtpqardl(d_y ~ d_x1, d, "id", "t", lr = c("L_y", "d_x1"), model = "dfe"),
               "both in 'formula'")
  expect_error(xtpqardl(d_y ~ d_x1, d, "id", "t", lr = c("L_y", "L_x1", "L_x1"), model = "dfe"),
               "duplicated")
})
