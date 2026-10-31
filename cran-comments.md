## xtpqardl 1.0.5

CRAN version: 1.0.4.

### Changes

* DFE: one joint covariance estimator for every block (new argument `dfe_cov`, default the Koenker and Bassett (1982) iid form); estimator stored and printed.
* `wald_test()`: Koenker-Bassett tests restrict the economic slopes only; `all = TRUE` restores the full restriction set.
* Lags built by the time index; duplicated (id, time) rows and time gaps are errors (`allow_gaps = TRUE` keeps panels with gaps).
* Lag selection penalises all estimated coefficients.
* Hausman test per quantile, NA with a warning when the difference matrix is not positive definite (`hausman = "ginv"` on request); `summary()` prints it.
* Collinear regressors and non-numeric time stop with informative errors; `log()`, `I()` and interaction terms are expanded with `model.matrix`.
* Warnings for every covariance estimator switch, dropped panel and NA block.
* Corrected the Koenker and Bassett (1982) reference (doi:10.2307/1912528).
* Removed text arrows from printed output.
* New tests pinned to hand computations.

### Test environments

* R 4.3.3 on Linux (x86_64), quantreg 6.1. No R-devel was available in the build container; the package was checked on R 4.3.3 only.

### R CMD check results

`R CMD check --as-cran --no-manual` on R 4.3.3: 0 errors, 0 warnings, 1 NOTE. The NOTE ("checking CRAN incoming feasibility": the GitHub URLs in DESCRIPTION returned HTTP 403) comes from the build container's outbound proxy; the URLs resolve in a browser.
