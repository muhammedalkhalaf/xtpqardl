# xtpqardl

Panel Quantile Autoregressive Distributed Lag Model for R

## Overview

The `xtpqardl` package provides functions for estimating Panel Quantile ARDL 
(PQARDL) models. It combines the panel ARDL methodology of Pesaran, Shin, and 
Smith (1999) with quantile regression to allow for heterogeneous effects across 
the conditional distribution of the response variable.

## Installation

```r
# Install from CRAN (when available)
install.packages("xtpqardl")

# Or install development version from GitHub
```

## Usage

```r
library(xtpqardl)

# Load example data
data(pqardl_sample)

# Estimate PQARDL model at multiple quantiles
fit <- xtpqardl(
  formula = d_y ~ d_x1 + d_x2,
  data = pqardl_sample,
  id = "country",
  time = "year",
  lr = c("L_y", "x1", "x2"),
  tau = c(0.25, 0.50, 0.75),
  model = "pmg"
)

# View results
summary(fit)

# Test parameter equality across quantiles
wald_test(fit)

# Compute impulse response function
irf <- compute_irf(fit, horizon = 20)
print(irf)
```

## Key Features

- **Multiple estimators**: Pooled Mean Group (PMG), Mean Group (MG), and 
  Dynamic Fixed Effects (DFE)
- **Multiple quantiles**: Estimate at any set of quantiles simultaneously
- **Long-run parameters**: Compute cointegrating coefficients β(τ)
- **Error correction**: Speed of adjustment ρ(τ) with convergence diagnostics
- **Half-life**: Time to close 50% of disequilibrium
- **Wald tests**: Test for parameter equality across quantiles
- **IRF**: Impulse response function by quantile
- **Lag selection**: Automatic BIC/AIC lag order selection

## References

- Pesaran MH, Shin Y, Smith RP (1999). "Pooled Mean Group Estimation of Dynamic 
  Heterogeneous Panels." *Journal of the American Statistical Association*, 
  94(446), 621-634. [doi:10.1080/01621459.1999.10474156](https://doi.org/10.1080/01621459.1999.10474156)

- Cho JS, Kim TH, Shin Y (2015). "Quantile Cointegration in the Autoregressive 
  Distributed-Lag Modeling Framework." *Journal of Econometrics*, 188(1), 281-300. 
  [doi:10.1016/j.jeconom.2015.02.030](https://doi.org/10.1016/j.jeconom.2015.02.030)

- Bildirici M, Kayikci F (2022). "Uncertainty, Renewable Energy, and CO2 
  Emissions in Top Renewable Energy Countries: A Panel Quantile Regression 
  Approach." *Energy*, 247, 124303. 
  [doi:10.1016/j.energy.2022.124303](https://doi.org/10.1016/j.energy.2022.124303)

## Author

## License
GPL-3
