## Tests for the 1.0.5 corrections (audit of 1.0.4). Every pinned value is a
## hand computation written out in the test (rq.fit, qr, lm) or an oracle
## from quantreg; nothing compares the package to itself.

sim_panel <- function(N = 4, TT = 20, seed = 7) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) {
    x <- cumsum(rnorm(TT))
    y <- numeric(TT)
    for (t in 2:TT) y[t] <- y[t - 1] - 0.4 * (y[t - 1] - 0.8 * x[t - 1]) +
      0.3 * (x[t] - x[t - 1]) + rnorm(1, 0, 0.3)
    data.frame(id = i, t = seq_len(TT), d_y = c(NA, diff(y)), d_x = c(NA, diff(x)),
               L_y = c(NA, y[-TT]), L_x = c(NA, x[-TT]))
  }))
}

## quantreg's iid sparsity recipe (summary.rq, se = "iid"), written out
hand_sparsity <- function(resid, n, p, tau) {
  eps <- .Machine$double.eps^(1/2)
  pz <- sum(abs(resid) < eps)
  h <- max(p + 1, ceiling(n * quantreg::bandwidth.rq(tau, n, hs = TRUE)))
  ir <- (pz + 1):(h + pz + 1)
  ord.resid <- sort(resid[order(abs(resid))][ir])
  xt <- ir / (n - p)
  unname(quantreg::rq(ord.resid ~ xt)$coef[2])
}

test_that("DFE joint covariance is the Koenker-Bassett (1982) Theorem 3.1 iid form, hand-built", {
  d <- sim_panel()
  taus <- c(0.25, 0.75)
  fit <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = taus, model = "dfe")
  expect_equal(fit$dfe_cov, "iid")
  ## hand design: slopes, intercept, N - 1 contrasts (panel 1 is the reference)
  dd <- d[complete.cases(d), ]
  XX <- cbind(dd$L_y, dd$L_x, dd$d_x, 1, outer(dd$id, 2:4, "=="))
  n <- nrow(XX); p <- ncol(XX)
  XtXinv <- solve(crossprod(XX))
  b <- s <- vector("list", 2)
  for (a in 1:2) {
    f <- quantreg::rq.fit(XX, dd$d_y, tau = taus[a], method = "br")
    b[[a]] <- f$coefficients
    s[[a]] <- hand_sparsity(f$residuals, n, p, taus[a])
  }
  V <- matrix(NA_real_, 2 * p, 2 * p)
  for (a in 1:2) for (bb in 1:2) {
    w <- min(taus[a], taus[bb]) - taus[a] * taus[bb]
    V[(a - 1) * p + 1:p, (bb - 1) * p + 1:p] <- w * s[[a]] * s[[bb]] * XtXinv
  }
  expect_equal(unname(fit$dfe_coef), cbind(b[[1]], b[[2]]), tolerance = 1e-10)
  expect_equal(unname(fit$dfe_coef_V), V, tolerance = 1e-10)
  ## the diagonal blocks equal quantreg's summary.rq(se = "iid") covariance
  sm <- summary(fit$rq_fits[[1]], se = "iid", covariance = TRUE)$cov
  expect_equal(unname(fit$dfe_coef_V[1:p, 1:p]), unname(sm), tolerance = 1e-10)
  ## rho_V is the (L_y, L_y) element of every block; beta_V the delta transform
  rho_V <- matrix(c(V[1, 1], V[1, p + 1], V[p + 1, 1], V[p + 1, p + 1]), 2)
  expect_equal(unname(fit$rho_V), rho_V, tolerance = 1e-10)
  G <- lapply(1:2, function(a) { r <- b[[a]][1]; th <- b[[a]][2]; cbind(th / r^2, -1 / r) })
  beta_V <- matrix(NA_real_, 2, 2)
  for (a in 1:2) for (bb in 1:2) {
    Cab <- V[(a - 1) * p + 1:2, (bb - 1) * p + 1:2]
    beta_V[a, bb] <- G[[a]] %*% Cab %*% t(G[[bb]])
  }
  expect_equal(unname(fit$beta_V), beta_V, tolerance = 1e-10)
  expect_equal(as.numeric(fit$beta_mg), -c(b[[1]][2] / b[[1]][1], b[[2]][2] / b[[2]][1]), tolerance = 1e-10)
  ## never mixed: the cross blocks use the same sparsities as the diagonal
  expect_equal(fit$rho_V[1, 2], (0.25 - 0.25 * 0.75) * s[[1]] * s[[2]] * XtXinv[1, 1], tolerance = 1e-10)
})

test_that("nid and ker DFE covariances are quantreg's H^{-1} J H^{-1} in every block", {
  d <- sim_panel()
  taus <- c(0.25, 0.5, 0.75)
  for (type in c("nid", "ker")) {
    fit <- suppressWarnings(xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = taus,
                                     model = "dfe", dfe_cov = type))
    expect_equal(fit$dfe_cov, type)
    sm <- lapply(fit$rq_fits, function(x) suppressWarnings(summary(x, se = type, covariance = TRUE)))
    p <- nrow(fit$dfe_coef)
    for (a in 1:3) for (b in 1:3) {
      w <- min(taus[a], taus[b]) - taus[a] * taus[b]
      ref <- w * sm[[a]]$Hinv %*% sm[[a]]$J %*% sm[[b]]$Hinv
      expect_equal(unname(fit$dfe_coef_V[(a - 1) * p + 1:p, (b - 1) * p + 1:p]), unname(ref), tolerance = 1e-10)
    }
    expect_equal(unname(fit$dfe_coef_V[1:p, 1:p]), unname(sm[[1]]$cov), tolerance = 1e-10)
  }
})

test_that("Koenker-Bassett joint test restricts only the economic slopes by default", {
  d <- sim_panel(N = 5, TT = 25)
  taus <- c(0.25, 0.5, 0.75)
  fit <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = taus, model = "dfe")
  w <- wald_test(fit)
  kb <- w$kb_test
  expect_equal(kb$restricted, c("L_y", "L_x", "d_x"))
  expect_equal(kb$joint[1, "ndf"], 3 * (length(taus) - 1))
  expect_equal(kb$joint[1, "ddf"], fit$n_obs * length(taus) - 3 * (length(taus) - 1))
  expect_equal(nrow(kb$individual), 3)
  expect_false(kb$all)
  ## hand: Tn/ndf with consecutive differences on the three slopes
  B <- fit$dfe_coef; V <- fit$dfe_coef_V; p <- nrow(B); m <- 3
  E <- cbind(diag(3), matrix(0, 3, p - 3))
  D <- kronecker(diff(diag(m)), E)
  Db <- D %*% as.numeric(B)
  Tn <- as.numeric(t(Db) %*% solve(D %*% V %*% t(D), Db)) / nrow(D)
  expect_equal(kb$joint[1, "Tn"], Tn, tolerance = 1e-10)
  expect_equal(kb$joint[1, "pvalue"], pf(Tn, nrow(D), fit$n_obs * m - nrow(D), lower.tail = FALSE), tolerance = 1e-10)
  ## all = TRUE adds the N - 1 panel contrasts
  wa <- wald_test(fit, all = TRUE)
  expect_equal(wa$kb_test$joint[1, "ndf"], (3 + 4) * (m - 1))
  expect_equal(length(wa$kb_test$restricted), 7)
  expect_true(wa$kb_test$all)
  expect_output(print(w), "economic slopes")
  expect_output(print(wa), "panel contrasts")
  expect_false(any(grepl("=>", capture.output(print(w)), fixed = TRUE)))
})

test_that("lags are built by the time index: panel with a gap", {
  d <- sim_panel(N = 4, TT = 30)
  dg <- d[!(d$id == 1 & d$t == 15), ]
  expect_error(xtpqardl(d_y ~ d_x, dg, "id", "t", lr = c("L_y", "L_x"), tau = 0.5,
                        model = "dfe", q = 2), "gap in 't'.*allow_gaps")
  fit <- suppressWarnings(xtpqardl(d_y ~ d_x, dg, "id", "t", lr = c("L_y", "L_x"), tau = 0.5,
                                   model = "dfe", q = 2, allow_gaps = TRUE))
  ## hand: the lag of d_x at time t is the value at time t - 1 of the same panel
  dg$L1 <- NA_real_
  for (i in unique(dg$id)) {
    s <- dg$id == i
    dg$L1[s] <- dg$d_x[s][match(dg$t[s] - 1, dg$t[s])]
  }
  dd <- dg[complete.cases(dg[, c("d_y", "L_y", "L_x", "d_x", "L1")]), ]
  XX <- cbind(dd$L_y, dd$L_x, dd$d_x, dd$L1, 1, outer(dd$id, 2:4, "=="))
  hb <- quantreg::rq.fit(XX, dd$d_y, tau = 0.5, method = "br")$coefficients
  expect_equal(fit$n_obs, nrow(dd))
  expect_equal(as.numeric(fit$rho_mg), hb[1], tolerance = 1e-8)
  expect_equal(as.numeric(fit$sr_mg), hb[3:4], tolerance = 1e-8)
  ## row-position lags would use a value two periods old after the gap
  expect_false(isTRUE(all.equal(dg$L1[dg$id == 1 & dg$t == 16], dg$d_x[dg$id == 1 & dg$t == 14])))
  ## the same for MG
  mg <- suppressWarnings(xtpqardl(d_y ~ d_x, dg, "id", "t", lr = c("L_y", "L_x"), tau = 0.5,
                                  model = "mg", q = 2, allow_gaps = TRUE))
  d1 <- dd[dd$id == 1, ]
  h1 <- quantreg::rq.fit(cbind(d1$L_y, d1$L_x, d1$d_x, d1$L1, 1), d1$d_y, tau = 0.5, method = "br")$coefficients
  expect_equal(mg$rho_all[1, 1], h1[1], tolerance = 1e-8)
  ## duplicated (id, time) rows are an error
  expect_error(xtpqardl(d_y ~ d_x, rbind(d, d[7, ]), "id", "t", lr = c("L_y", "L_x"), tau = 0.5,
                        model = "dfe"), "duplicated")
  ## non-numeric time is an error
  dc <- d; dc$t <- as.character(dc$t)
  expect_error(xtpqardl(d_y ~ d_x, dc, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "dfe"),
               "must be numeric")
})

test_that("collinear regressors, transformed formulas, and NA covariance blocks", {
  d <- sim_panel()
  d2 <- d; d2$d_x2 <- 2 * d$d_x
  expect_error(xtpqardl(d_y ~ d_x + d_x2, d2, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "dfe"),
               "collinear.*d_x2")
  expect_error(xtpqardl(d_y ~ d_x + d_x2, d2, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "mg"),
               "collinear.*panel 1.*d_x2")
  ## log() and I() terms are expanded with model.matrix
  d$xpos <- exp(d$d_x / 2)
  fl <- xtpqardl(d_y ~ log(xpos) + I(d_x^2), d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "mg")
  expect_equal(fl$indepvars, c("log(xpos)", "I(d_x^2)"))
  dd <- d[complete.cases(d), ]; d1 <- dd[dd$id == 1, ]
  h1 <- quantreg::rq.fit(cbind(d1$L_y, d1$L_x, log(d1$xpos), d1$d_x^2, 1), d1$d_y, tau = 0.5, method = "br")$coefficients
  expect_equal(fl$rho_all[1, 1], h1[1], tolerance = 1e-8)
  expect_equal(unname(fl$sr_all[1, ]), unname(h1[3:4]), tolerance = 1e-8)
  ## a covariance that cannot be estimated gives NA blocks with a warning, not a
  ## substitute: y an exact linear function of the design (all residuals zero)
  de <- d; de$d_y <- -0.4 * de$L_y + 0.3 * de$L_x + 0.2 * de$d_x + de$id
  expect_warning(fe <- xtpqardl(d_y ~ d_x, de, "id", "t", lr = c("L_y", "L_x"), tau = c(0.25, 0.75), model = "dfe"),
                 "could not be estimated")
  expect_true(all(is.na(fe$rho_V)))
  expect_false(fe$cross_quantile_cov)
  expect_false(all(fe$panel_cov_status$ok))
  w <- suppressWarnings(wald_test(fe))
  expect_true(is.na(w$rho_test$statistic))
  expect_true(is.na(w$kb_test$joint[1, "Tn"]))
})

test_that("Hausman test is NA with a warning when the difference matrix is not positive definite", {
  d <- sim_panel(N = 6, TT = 30)
  ## equal weights: V(PMG) = V(MG), so the difference is the zero matrix
  expect_warning(pe <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = c(0.25, 0.75),
                                model = "pmg", pmg_weights = "equal"), "not positive definite")
  expect_s3_class(pe$hausman, "data.frame")
  expect_true(all(is.na(pe$hausman$statistic)))
  expect_equal(pe$hausman$method, c("strict", "strict"))
  expect_output(print(summary(pe)), "Hausman test")
  expect_output(print(summary(pe)), "not defined")
  ## ginv only on request; on the zero matrix there is no positive eigenvalue
  ## to invert, so the statistic stays NA and the method is recorded
  pg <- suppressWarnings(xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = c(0.25, 0.75),
                                  model = "pmg", pmg_weights = "equal", hausman = "ginv"))
  expect_equal(pg$hausman$method, c("ginv", "ginv"))
  expect_true(all(is.na(pg$hausman$statistic)))
  ## ginv on a hand-built indefinite difference: eigenvalues (2, -1), d = (1, 1)
  expect_warning(h <- xtpqardl:::.hausman_one(c(1, 1), diag(c(2, -1)), 0.5, method = "ginv"), "discarded")
  expect_equal(h$statistic, 1 / 2)
  expect_equal(h$df, 1L)
  expect_warning(hs <- xtpqardl:::.hausman_one(c(1, 1), diag(c(2, -1)), 0.5, method = "strict"), "not positive definite")
  expect_true(is.na(hs$statistic))
  ## a positive definite case reproduces the hand quadratic form
  pm <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "pmg")
  mg <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "mg")
  ## the difference V(MG) - V(PMG) is positive for this seeded panel, so the
  ## statistic must be defined (asserted, not conditional)
  expect_gt(mg$beta_V[1, 1] - pm$beta_V[1, 1], 0)
  expect_false(is.na(pm$hausman$statistic[1]))
  dlt <- as.numeric(mg$beta_mg) - as.numeric(pm$beta_mg)
  h <- dlt^2 / (mg$beta_V[1, 1] - pm$beta_V[1, 1])
  expect_equal(pm$hausman$statistic[1], h, tolerance = 1e-10)
  expect_equal(pm$hausman$p.value[1], pchisq(h, 1, lower.tail = FALSE), tolerance = 1e-10)
  expect_output(print(summary(pm)), "Hausman")
})

test_that("MG with fewer than three panels warns", {
  d <- sim_panel(N = 2, TT = 30)
  expect_warning(xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = 0.5, model = "mg"),
                 "Only 2 panel")
})

test_that("lag selection penalises all estimated coefficients and uses least squares", {
  d <- sim_panel(N = 4, TT = 30)
  prep <- xtpqardl:::.prepare_data(d_y ~ d_x, d, "id", "t", allow_gaps = FALSE)
  sel <- xtpqardl:::.select_lag_order(prep$data, prep$depvar, prep$indepvars, c("L_y", "L_x"),
                                       "id", "t", pmax = 2, qmax = 2, criterion = "bic")
  tab <- sel$table
  ## p = 1, q = 1: 4 coefficients per panel (L_y, L_x, d_x, constant), 4 panels
  r11 <- tab[tab$p == 1 & tab$q == 1, ]
  expect_equal(r11$K, 4 * 4)
  ## hand: pooled RSS and n from lm on each panel
  dd <- d[complete.cases(d), ]
  rss <- 0; n <- 0
  for (i in 1:4) {
    di <- dd[dd$id == i, ]
    fit <- lm(d_y ~ L_y + L_x + d_x, data = di)
    rss <- rss + sum(residuals(fit)^2); n <- n + nrow(di)
  }
  expect_equal(r11$ic, n * log(rss / n) + 16 * log(n), tolerance = 1e-10)
  ## p = 2, q = 2: 6 coefficients per panel
  expect_equal(tab[tab$p == 2 & tab$q == 2, "K"], 6 * 4)
})

test_that("covariance estimator switches are reported with a warning and recorded", {
  d <- sim_panel(N = 4, TT = 20)
  fit <- xtpqardl(d_y ~ d_x, d, "id", "t", lr = c("L_y", "L_x"), tau = c(0.25, 0.75), model = "dfe")
  expect_equal(fit$panel_cov_status$se, c("iid", "iid"))
  expect_output(print(fit), "dfe_cov = \"iid\"")
  expect_output(print(summary(fit)), "Koenker and Bassett")
  ## the printed output has no text arrows
  out <- c(capture.output(print(fit)), capture.output(print(summary(fit))),
           capture.output(print(wald_test(fit))))
  expect_false(any(grepl("=>|->", out)))
})
