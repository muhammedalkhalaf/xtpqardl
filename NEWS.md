# xtpqardl 1.0.5

## Breaking changes

* `model = "dfe"`: the standard errors of rho and beta and the cross-quantile blocks of `rho_V` and `beta_V` come from one joint covariance matrix, the Koenker and Bassett (1982, Theorem 3.1) form with the scalar sparsity of `summary.rq(se = "iid")` (new argument `dfe_cov`, default `"iid"`). In 1.0.4 the diagonal blocks were the Powell kernel sandwich (`se = "ker"`) and the off-diagonal blocks the Hendricks-Koenker form (`se = "nid"`). The estimator is stored in `dfe_cov` and printed by `print()` and `summary()`. `dfe_cov = "nid"` and `"ker"` use the corresponding `summary.rq` estimator for every block; the man page documents the Monte Carlo size of each choice.
* `wald_test()`: the Koenker-Bassett tests in `kb_test` restrict only the economic slopes (rho, theta, short run) and are computed from the joint matrix of the fit; in 1.0.4 they were `quantreg::anova(..., se = "nid")` and also restricted the N - 1 panel-intercept contrasts. New argument `all = TRUE` restricts every non-intercept coefficient; with `dfe_cov = "nid"` this reproduces the 1.0.4 values when the quantile regression solutions are unique (when they are not, the ordering of the panel contrasts can change the Hendricks-Koenker densities and so the values). `kb_test` now holds `restricted`, `se` and `all`; the printed labels state which coefficients are restricted.
* Lags of the differenced variables are built by the time index (the value at `time - j` within the panel), not by row position. Duplicated (`id`, `time`) rows are an error. A gap in a panel's time index is an error unless the new argument `allow_gaps = TRUE` is given; rows whose lagged period is missing are then dropped. The `time` variable must be numeric.
* `lagsel`: the information criteria count the coefficients of all panel regressions (N times the per-panel count); in 1.0.4 only one panel's coefficients were counted. The man page states that the selection is by least squares.
* `model = "pmg"`: `hausman` is a data frame with one row per quantile (statistic, df, p-value, number of non-positive eigenvalues, method). The test is `NA` with a warning when V(MG) - V(PMG) is not positive definite (new argument `hausman = "strict"`, the default); `hausman = "ginv"` uses a Moore-Penrose inverse with a warning. In 1.0.4 one statistic over all quantiles was computed with the PMG cross-quantile blocks set to zero and negative eigenvalues dropped silently.
* `summary()` prints the Hausman table; in 1.0.4 the printing code was never reached.
* `rq_fits` (DFE): 1.0.4 stored `rq(y ~ XX)` fits with coefficients named `(Intercept)`, `XXL_y`, ..., and panel-contrast levels in character order. 1.0.5 stores `rq(y ~ XX - 1)` fits whose coefficients are the slopes, then `(Intercept)`, then the panel contrasts, named by the id variable followed by the panel level (for example `countryCountry_B`), in panel order. Code that indexes `coef(fit$rq_fits[[1]])` by name or position must be updated.
* The warning for fewer than three panels (`N < 3`) is issued for `model = "pmg"` as well as `model = "mg"`, since both use the same panel code.

## Other changes

* Exactly collinear regressors stop with an error naming the columns (and the panel for MG and PMG).
* Formulas with `log()`, `I()` and interactions are expanded with `model.matrix` before the lags are built; the column names of the expanded design are the regressor names.
* Warnings are emitted once when the per-panel covariance estimator differs from `summary.rq(se = "iid")` (MG, PMG), when a panel is dropped from the PMG pooling, when a DFE covariance block cannot be estimated (the block is `NA`), and when `summary.rq` reports non-positive density estimates. `panel_cov_status` is also returned for DFE.
* `model = "mg"` with fewer than three panels warns.
* New return elements `dfe_coef`, `dfe_coef_V`, `dfe_slopes`, `n_dropped`, `time_step`; `summary()` prints the number of rows dropped for missing values.
* The rank-one property of the asymptotic covariance of theta (Cho, Kim and Shin 2015, Theorem 1(ii)) and its consequence for joint tests on several theta coefficients are documented in `xtpqardl` and `wald_test`.
* The Koenker and Bassett (1982) reference is "Robust Tests for Heteroscedasticity Based on Regression Quantiles", Econometrica 50(1), 43-61, doi:10.2307/1912528 (1.0.4 cited "Tests of Linear Hypotheses and L1 Estimation", doi:10.2307/1913398); Hausman (1978) added to the references.
* `print.wald_test()` no longer writes text arrows.
* Tests: the 1.0.4 test of the `quantreg::anova` reproduction now uses `dfe_cov = "nid"` and `all = TRUE`; the PMG short-panel test accepts the new switch warning. New tests pin the iid joint covariance, the "nid" and "ker" blocks, the restriction count of the joint test, time-aligned lags on a panel with a gap, the collinearity and formula handling, the `NA` covariance path, the Hausman rule and the lag-selection penalty to hand computations.

# xtpqardl 1.0.4

* Bug fix: `wald_test(fit, joint = FALSE)` failed inside `print()` with "$ operator is invalid for atomic vectors" whenever the covariance matrix had `NA` blocks (every DFE and PMG fit in 1.0.3), because the helper returned one flat list instead of one list per variable. Fixed, and `print()` now guards against malformed elements.
* New: for `model = "dfe"` the cross-quantile covariance blocks of `beta_V` and `rho_V` are estimated with the Koenker and Bassett (1982) formula (all quantiles share the pooled design), so the cross-quantile Wald tests of `wald_test()` are now available for DFE. `wald_test()` additionally returns `kb_test`, the `quantreg::anova(..., se = "nid")` equality tests (joint and per coefficient) on the pooled quantile-regression coefficients; these reproduce `anova(rq_tau1, rq_tau2, ...)` on the equivalent regression `d_y ~ L_y + L_x + d_x + factor(id)`. The `rq` fits are kept in `rq_fits`.
* PMG with short panels: the per-panel covariance used for the inverse-variance weights comes from `summary.rq`, which fails when the degrees of freedom per panel are very few (the study data had n - k = 2); in 1.0.3 every panel was then silently dropped and PMG stopped with "needs at least two panels with a usable long-run covariance". The function now tries the other `summary.rq` estimators and an iid formula with a kernel sparsity estimate before giving up, reports the status and reason for every panel and quantile (`panel_cov_status`, also for MG), explains the degrees-of-freedom problem in the error message, and offers `pmg_weights = "equal"` (unweighted pooling).
* `wald_test()`: the singularity check `det(RVR') < 1e-15` declared the test "Singular variance matrix" whenever the dependent variable had a small scale (the determinant of a matrix of small variances is tiny even when it is well conditioned); it now uses the reciprocal condition number (`rcond < 1e-12`).
* Input checks: `lr` must not contain the dependent variable of `formula`, nor any of its regressors, and all variables must be numeric; the previous behaviour was "subscript out of bounds" or a silent run with a wrong error correction term.
* Tests added for all of the above.

# xtpqardl 1.0.3

* Bug fix: the long-run variables in `lr` were lagged a second time, so with `lr = c("L_y", ...)` the error correction term used y(t-2). They now enter exactly as supplied.
* Bug fix: `model = "pmg"` returned the Mean Group estimator. PMG now pools the long-run coefficients by minimum distance (inverse delta-method covariance weights) and re-estimates each panel with the pooled long run imposed; a Hausman test of long-run homogeneity (MG against PMG) is returned in `hausman`.
* Bug fix: `model = "dfe"` returned fixed placeholder variances (0.01 on the diagonal). Its standard errors now come from the Powell kernel sandwich covariance of the pooled quantile regression, transformed by the delta method.
* Bug fix: the half-life was ln(2)/|rho|; it is now the exact ln(0.5)/ln(1 + rho), defined for -1 < rho < 0.
* Cross-quantile covariance blocks that are not estimated (PMG and DFE) are now `NA`, so `wald_test()` reports them as unavailable instead of assuming independence.
* The MG and DFE point estimates agree with the Stata command xtpqardl (SSC, v1.0.4) on the same data (MG: rho = -0.5123, beta = 0.5946; DFE: rho = -0.4702, beta = 0.5642).
* Added unit tests.

# xtpqardl 1.0.2

* Corrected the DOI of Cho, Kim and Shin (2015) to 10.1016/j.jeconom.2015.05.003 (all occurrences).
* Removed a DOI attached to Bildirici and Kayikci (2022) that could not be verified in CrossRef; the citation text is unchanged.
* Fixed a stray brace in xtpqardl-package.Rd; Authors@R updated.

# xtpqardl 1.0.1

## New Features

* Initial CRAN release
* Panel Quantile ARDL estimation with PMG, MG, and DFE estimators
* Support for multiple quantiles
* Long-run cointegrating parameters estimation
* ECT speed of adjustment computation
* Half-life of adjustment
* Wald tests for parameter equality across quantiles
* Impulse response function computation
* BIC/AIC automatic lag selection

## References

* Based on Pesaran, Shin, and Smith (1999) panel ARDL methodology
* Implements quantile cointegration framework of Cho, Kim, and Shin (2015)
