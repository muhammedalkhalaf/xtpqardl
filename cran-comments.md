## xtpqardl 1.0.3

This release corrects the computations below; the 1.0.2 submission (reference metadata only) should be discarded in favour of this one.

* Bug fix: the long-run variables in `lr` were lagged a second time, so with `lr = c("L_y", ...)` the error correction term used y(t-2). They now enter exactly as supplied.
* Bug fix: `model = "pmg"` returned the Mean Group estimator. PMG now pools the long-run coefficients by minimum distance (inverse delta-method covariance weights) and re-estimates each panel with the pooled long run imposed; a Hausman test of long-run homogeneity (MG against PMG) is returned in `hausman`.
* Bug fix: `model = "dfe"` returned fixed placeholder variances (0.01 on the diagonal). Its standard errors now come from the Powell kernel sandwich covariance of the pooled quantile regression, transformed by the delta method.
* Bug fix: the half-life was ln(2)/|rho|; it is now the exact ln(0.5)/ln(1 + rho), defined for -1 < rho < 0.
* Cross-quantile covariance blocks that are not estimated (PMG and DFE) are now `NA`, so `wald_test()` reports them as unavailable instead of assuming independence.
* The MG and DFE point estimates agree with the Stata command xtpqardl (SSC, v1.0.4) on the same data (MG: rho = -0.5123, beta = 0.5946; DFE: rho = -0.4702, beta = 0.5642).
* Added unit tests.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
