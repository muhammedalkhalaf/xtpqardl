#' Panel Quantile Autoregressive Distributed Lag Model
#'
#' Estimate Panel Quantile ARDL (PQARDL) models that combine panel ARDL
#' methodology with quantile regression. Supports Pooled Mean Group (PMG),
#' Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators.
#'
#' @param formula A formula specifying the model. The response variable should
#'   be in first differences (e.g., \code{d.y ~ d.x1 + d.x2}), where short-run
#'   dynamics are estimated. Transformations such as \code{log(x)},
#'   \code{I(x^2)} and interactions \code{a:b} are expanded with
#'   \code{model.matrix} before the lags are built, so the lag orders in
#'   \code{q} refer to the columns of the expanded design (in the order of
#'   \code{colnames(model.matrix(formula, data))} without the intercept).
#' @param data A data frame containing panel data with variables specified in
#'   the formula and \code{lr} argument.
#' @param id Character string specifying the panel (cross-section) identifier
#'   variable name.
#' @param time Character string specifying the time variable name. It must be
#'   numeric (for example the year as an integer). Lags are built by matching
#'   \code{time - j} within each panel, never by row position; see Details.
#' @param lr Character vector of long-run level variable names. The first
#'   element should be the lagged dependent variable level (for the error
#'   correction term), and remaining elements are the long-run explanatory
#'   variables. The variables enter the regression exactly as supplied (no
#'   further lag is applied); supply lagged regressors, for example
#'   \code{L_x1}, if the error correction term should use \eqn{x_{t-1}}.
#'   The names must be columns of \code{data} (no transformations).
#' @param tau Numeric vector of quantiles to estimate, each in (0,1).
#'   Default is \code{c(0.25, 0.50, 0.75)}.
#' @param p Integer specifying the autoregressive lag order for the dependent
#'   variable. Default is 1.
#' @param q Integer or integer vector specifying the distributed lag order(s)
#'   for explanatory variables. If a single integer, the same lag order is
#'   applied to all variables. Default is 1.
#' @param model Character string specifying the estimation method:
#'   \code{"pmg"} for Pooled Mean Group (default), \code{"mg"} for Mean Group,
#'   or \code{"dfe"} for Dynamic Fixed Effects.
#' @param lagsel Character string for automatic lag selection. If \code{"bic"}
#'   or \code{"aic"}, the lag orders are selected by least squares on the
#'   panel regressions, see Details. Default is \code{NULL} (no automatic
#'   selection).
#' @param pmax Maximum p to consider in lag selection. Default is 4.
#' @param qmax Maximum q to consider in lag selection. Default is 4.
#' @param constant Logical. Include a constant term? Default is \code{TRUE}.
#' @param pmg_weights Character, for \code{model = "pmg"} only.
#'   \code{"inverse_variance"} (default) weights each panel's long-run
#'   estimate by the inverse of its delta-method covariance, as described in
#'   Details; it requires that covariance to be estimable for at least two
#'   panels per quantile. \code{"equal"} pools with equal weights (the pooled
#'   long run is then the simple average of the panel estimates and its
#'   covariance the mean-group covariance); use it when the panels are too
#'   short for the covariances to be estimated.
#' @param dfe_cov Character, for \code{model = "dfe"} only: the estimator of
#'   the joint (all quantiles) covariance matrix of the pooled quantile
#'   regression coefficients. \code{"iid"} (default) is the Koenker and
#'   Bassett (1982, Theorem 3.1) form with the scalar sparsity estimate of
#'   \code{summary.rq(se = "iid")}; \code{"nid"} uses the Hendricks-Koenker
#'   local densities of \code{summary.rq(se = "nid")}; \code{"ker"} the
#'   Powell kernel densities of \code{summary.rq(se = "ker")}. One estimator
#'   is used for every block of the matrix; see Details for the choice.
#' @param allow_gaps Logical. If \code{FALSE} (default) the function stops
#'   when a panel has a gap in \code{time}. If \code{TRUE} the panels with
#'   gaps are kept: lags are matched by the time index, and a row whose
#'   lagged period is missing is dropped from that panel's regression.
#' @param hausman Character, for \code{model = "pmg"} only: \code{"strict"}
#'   (default) returns \code{NA} with a warning when the difference
#'   \eqn{V(MG) - V(PMG)} is not positive definite; \code{"ginv"} inverts it
#'   with a Moore-Penrose inverse over the eigenvalues larger than
#'   \code{1e-8} times the largest, with degrees of freedom equal to that
#'   rank, and warns.
#'
#' @return An object of class \code{"xtpqardl"} containing:
#' \describe{
#'   \item{beta_mg}{Matrix of mean group long-run coefficients across quantiles}
#'   \item{rho_mg}{Vector of mean group ECT speed of adjustment by quantile}
#'   \item{halflife_mg}{Half-life \eqn{\ln(0.5)/\ln(1+\rho)} of the mean
#'     group speed of adjustment, by quantile}
#'   \item{sr_mg}{Matrix of mean group short-run coefficients}
#'   \item{phi_mg}{Matrix of mean group AR coefficients (if p > 1)}
#'   \item{hausman}{For \code{model = "pmg"}, a data frame with one row per
#'     quantile: the Hausman test of long-run homogeneity (MG against PMG)
#'     at that quantile, with \code{statistic}, \code{df}, \code{p.value},
#'     the number of non-positive eigenvalues of the difference matrix and
#'     the method used}
#'   \item{beta_V}{Variance-covariance matrix for beta_mg. For MG and DFE the
#'     full matrix is estimated (for DFE every block comes from the estimator
#'     named in \code{dfe_cov}, see Details); for PMG the cross-quantile
#'     blocks are \code{NA}}
#'   \item{rho_V}{Variance-covariance matrix for rho_mg}
#'   \item{dfe_cov}{For \code{model = "dfe"}, the covariance estimator used
#'     (\code{"iid"}, \code{"nid"} or \code{"ker"})}
#'   \item{dfe_coef}{For \code{model = "dfe"}, the matrix of pooled
#'     quantile-regression coefficients (rows: coefficients, columns: tau)}
#'   \item{dfe_coef_V}{For \code{model = "dfe"}, the joint covariance matrix
#'     of \code{dfe_coef} stacked by quantile, used by
#'     \code{\link{wald_test}} for the Koenker-Bassett tests}
#'   \item{dfe_slopes}{For \code{model = "dfe"}, the names of the economic
#'     slopes (speed of adjustment, long-run levels, short run), as opposed
#'     to the intercept and the panel-intercept contrasts}
#'   \item{rq_fits}{For \code{model = "dfe"}, the list of \code{rq} fits (one
#'     per quantile) of the pooled regression}
#'   \item{panel_cov_status}{A data frame with one row per panel and quantile
#'     (MG, PMG) or per quantile (DFE) reporting whether the covariance
#'     could be estimated, which estimator was used and why not}
#'   \item{beta_all}{Matrix of per-panel long-run coefficients}
#'   \item{rho_all}{Matrix of per-panel ECT coefficients}
#'   \item{halflife_all}{Matrix of per-panel half-life values}
#'   \item{tau}{Vector of estimated quantiles}
#'   \item{p}{AR lag order used}
#'   \item{q}{Distributed lag order(s) used}
#'   \item{model}{Estimation method used}
#'   \item{n_obs}{Total number of observations}
#'   \item{n_dropped}{Number of rows dropped because a regressor, the
#'     response or a lag was missing}
#'   \item{n_panels}{Number of panels}
#'   \item{valid_panels}{Number of successfully estimated panels}
#'   \item{depvar}{Dependent variable name}
#'   \item{lrvars}{Long-run variable names}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' The PQARDL model extends the standard panel ARDL framework to allow for
#' heterogeneous effects across the conditional distribution of the response
#' variable. The error correction representation is:
#'
#' \deqn{\Delta y_{it} = \rho_i(\tau) \cdot ECT_{i,t-1} + \sum_{j=1}^{p-1}
#' \phi_{ij} \Delta y_{i,t-j} + \sum_{m=0}^{q-1} \theta_{im} \Delta x_{i,t-m}
#' + \varepsilon_{it}(\tau)}
#'
#' where \eqn{ECT = y_{i,t-1} - \beta(\tau)' X} is the error correction term
#' built from the \code{lr} variables as supplied, \eqn{\rho(\tau)} is the
#' speed of adjustment (negative for convergence), and \eqn{\beta(\tau)} are
#' the long-run parameters, \eqn{\beta = -\theta/\rho} (the plug-in estimator
#' of Cho, Kim and Shin 2015, equation (3)).
#'
#' \strong{Lags, gaps and missing values.} Within each panel the lag
#' \eqn{j} of a variable at time \eqn{t} is the value at time \eqn{t - j},
#' where the unit step is the smallest positive time difference in the data.
#' Duplicated (\code{id}, \code{time}) pairs are an error. A gap in a panel's
#' time index is an error unless \code{allow_gaps = TRUE}, in which case the
#' rows whose lagged period is missing are dropped. Rows with a missing value
#' in the response, a regressor or a lag are dropped panel by panel (the
#' leading rows lost to the user's own differencing and lagging are the
#' usual case); the number of dropped rows is returned in \code{n_dropped}.
#' Exactly collinear regressors stop with an error naming the columns.
#'
#' \code{model = "mg"} averages the panel estimates (Pesaran and Smith
#' covariance \eqn{\sum_i (b_i - \bar b)(b_i - \bar b)'/(N(N-1))}). With fewer
#' than three panels this covariance is unreliable (with N = 2 the audit's
#' Monte Carlo rejected 26 to 31 percent at nominal 5 percent) and the
#' function warns. \code{model = "pmg"} pools the long run by minimum
#' distance, weighting each panel by the inverse of its delta-method
#' covariance, then re-estimates every panel with the pooled long run
#' imposed; the speed of adjustment and the short run are the mean group of
#' that second stage. A Hausman test of long-run homogeneity (MG against PMG)
#' is computed quantile by quantile; the statistic is \eqn{d'(V_{MG} -
#' V_{PMG})^{-1} d} with \eqn{d} the difference of the long-run estimates at
#' that quantile. The quadratic form is defined only when the difference
#' matrix is positive definite (Hausman 1978 gives non-negative definiteness
#' asymptotically); see \code{hausman}.
#'
#' \strong{DFE covariance.} \code{model = "dfe"} pools all panels in one
#' quantile regression per quantile with an intercept and \eqn{N - 1}
#' panel-intercept contrasts. All quantiles share the design \eqn{X}, so the
#' joint covariance of the coefficient vectors \eqn{b(\tau_a)}, \eqn{b(\tau_b)}
#' is estimated by one formula for every block, with the estimator named in
#' \code{dfe_cov}:
#' \deqn{Cov(b(\tau_a), b(\tau_b)) = (\min(\tau_a,\tau_b) - \tau_a\tau_b)
#' H_a^{-1} (X'X) H_b^{-1}.}
#' With \code{"iid"} (default) \eqn{H_a^{-1} = s(\tau_a) (X'X)^{-1}}, where
#' \eqn{s(\tau) = 1/f(F^{-1}(\tau))} is the sparsity estimated as in
#' \code{summary.rq(se = "iid")}; this is Theorem 3.1 of Koenker and Bassett
#' (1982), whose covariance has typical element
#' \eqn{[\min(u_i,u_j) - u_i u_j]/(f(Q(u_i)) f(Q(u_j)))} times
#' \eqn{D^{-1} = \lim (X'X/n)^{-1}}, under i.i.d. errors. With \code{"nid"}
#' and \code{"ker"}, \eqn{H_\tau = X' diag(f_i(\tau)) X} with the local density
#' estimates of \code{summary.rq(se = "nid")} (Hendricks-Koenker) or
#' \code{summary.rq(se = "ker")} (Powell kernel), including quantreg's
#' bandwidth rule (the bandwidth is halved while \eqn{\tau \pm h} leaves
#' (0, 1)). The diagonal blocks equal \code{summary.rq(se = dfe_cov)$cov}.
#' The default is \code{"iid"} because in the audit of version 1.0.4 (panel
#' designs with N = 5 to 20, T = 20 to 40, i.i.d. errors) the
#' Hendricks-Koenker densities under-stated the variance of the speed of
#' adjustment by 23 to 55 percent in the presence of the panel intercepts
#' (Monte Carlo variance 0.00161 against a mean "nid" estimate 0.00072 at
#' N = 20, T = 40), so that the "nid" equality tests rejected 9 to 23 percent
#' at nominal 5 percent for N at most 10 and 47 to 76 percent at N = 20,
#' whereas the "iid" form rejected 4.3, 5.0 and 7.0 percent in the three
#' designs checked. For the default in the N = 10, T = 40, k = 1 design (iid
#' errors, nominal 5 percent) \code{rho_test} rejected between 6 and 9 percent
#' in two independent Monte Carlo batches (R = 200 and R = 500); a third batch
#' (R = 300) gave 5.7 percent. The test is therefore mildly over-sized in that
#' cell, against 1.2 percent with the 1.0.4 mixed matrix. \code{"nid"} and \code{"ker"} remain available for users
#' who need them; they are not recommended for the cross-quantile tests of
#' panels with many intercepts. The same joint matrix is the basis of the
#' Koenker-Bassett tests of \code{\link{wald_test}}; if the estimator fails
#' for a quantile (for example an extreme quantile in a small sample), the
#' affected blocks are \code{NA} and a warning names the quantile and the
#' reason. The standard errors of \eqn{\rho} and \eqn{\beta} come from the
#' diagonal blocks and the delta method.
#'
#' \strong{Rank of the covariance of theta.} Cho, Kim and Shin (2015, Theorem
#' 1(ii)) show that under their normalisation the asymptotic covariance of
#' the long-run level coefficients \eqn{\gamma(\tau)} (\code{theta} here, the
#' coefficients of \code{lr[-1]}) has rank one, and they restrict their Wald
#' tests on \eqn{\gamma} to restriction matrices of rank one, because a
#' higher rank gives a degenerate null distribution. The finite-sample
#' sandwich blocks for \eqn{\theta} used here are of full rank, so a joint
#' test that restricts several \eqn{\theta} coefficients at once (as the
#' Koenker-Bassett joint test of \code{\link{wald_test}} does when
#' \eqn{k_x > 1}) relies on finite-sample full rank of a matrix whose limit
#' is singular under their assumptions. See \code{\link{wald_test}} for the
#' size numbers. Tests on a single \eqn{\theta_j} across quantiles (the
#' individual rows) are of rank \eqn{m - 1} in the quantile direction.
#'
#' \strong{Lag selection.} \code{lagsel} fits every (p, q) combination with
#' \eqn{p \le} \code{pmax}, \eqn{q \le} \code{qmax} (the same q for all
#' regressors) by least squares (\code{lm.fit}) panel by panel, sums the
#' residual sums of squares and the sample sizes over panels, and penalises
#' the total number of estimated coefficients (the sum over the panels of the
#' columns of their regressions): \eqn{BIC = n \log(RSS/n) + K \log n},
#' \eqn{AIC = n \log(RSS/n) + 2K}. The criterion is least squares, not the
#' quantile check loss, and applies to the mean regression.
#'
#' \strong{Short panels.} Each panel regression has \eqn{1 + k_x + (p-1) +
#' \sum q + 1} coefficients. When \eqn{T} is close to that number the
#' covariance of the panel estimates cannot be estimated by
#' \code{summary.rq(se = "iid")} (its sparsity estimate needs the quantiles at
#' \eqn{\tau \pm h}); the function then tries the other \code{summary.rq}
#' estimators and finally an iid formula with a kernel density estimate of
#' the sparsity, and warns once naming how many panel-quantile pairs switched
#' and to which estimator (details in \code{panel_cov_status}). If no
#' covariance is available for at least two panels, PMG stops with a message
#' listing the status of every panel; \code{pmg_weights = "equal"} or a
#' smaller model are the alternatives. Panels dropped from the PMG pooling
#' because their covariance is not invertible are named in a warning.
#'
#' @references
#' Pesaran MH, Shin Y, Smith RP (1999). "Pooled Mean Group Estimation of
#' Dynamic Heterogeneous Panels." \emph{Journal of the American Statistical
#' Association}, 94(446), 621-634. \doi{10.1080/01621459.1999.10474156}
#'
#' Cho JS, Kim TH, Shin Y (2015). "Quantile Cointegration in the
#' Autoregressive Distributed-Lag Modeling Framework." \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.05.003}
#'
#' Bildirici M, Kayikci F (2022). "Uncertainty, Renewable Energy, and CO2
#' Emissions in Top Renewable Energy Countries: A Panel Quantile Regression
#' Approach." \emph{Energy}, 247, 124303.
#'
#' Hausman JA (1978). "Specification Tests in Econometrics."
#' \emph{Econometrica}, 46(6), 1251-1271. \doi{10.2307/1913827}
#'
#' Koenker R, Bassett G (1978). "Regression Quantiles." \emph{Econometrica},
#' 46(1), 33-50. \doi{10.2307/1913643}
#'
#' Koenker R, Bassett G (1982). "Robust Tests for Heteroscedasticity Based on
#' Regression Quantiles." \emph{Econometrica}, 50(1), 43-61.
#' \doi{10.2307/1912528}
#'
#' @examples
#' \donttest{
#' # Load example panel data
#' data(pqardl_sample)
#'
#' # Estimate PQARDL model at 25th, 50th, and 75th quantiles
#' fit <- xtpqardl(
#'   formula = d_y ~ d_x1 + d_x2,
#'   data = pqardl_sample,
#'   id = "country",
#'   time = "year",
#'   lr = c("L_y", "x1", "x2"),
#'   tau = c(0.25, 0.50, 0.75),
#'   model = "pmg"
#' )
#'
#' # View results
#' summary(fit)
#'
#' # Wald test for parameter equality across quantiles
#' wald_test(fit)
#' }
#'
#' @importFrom stats model.frame model.matrix model.response var cov na.omit
#'   pnorm pchisq pf formula terms complete.cases as.formula coef lm.fit
#'   na.pass
#' @importFrom quantreg rq
#' @export
xtpqardl <- function(formula, data, id, time, lr,
                     tau = c(0.25, 0.50, 0.75),
                     p = 1, q = 1,
                     model = c("pmg", "mg", "dfe"),
                     lagsel = NULL,
                     pmax = 4, qmax = 4,
                     constant = TRUE,
                     pmg_weights = c("inverse_variance", "equal"),
                     dfe_cov = c("iid", "nid", "ker"),
                     allow_gaps = FALSE,
                     hausman = c("strict", "ginv")) {
  pmg_weights <- match.arg(pmg_weights)
  dfe_cov <- match.arg(dfe_cov)
  hausman <- match.arg(hausman)

  # Match arguments
  model <- match.arg(model)
  call <- match.call()

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!id %in% names(data)) {
    stop("Panel identifier '", id, "' not found in data")
  }
  if (!time %in% names(data)) {
    stop("Time variable '", time, "' not found in data")
  }
  if (!is.character(lr)) {
    stop("'lr' must be a character vector of variable names: the lagged level ",
         "of the dependent variable first, then the long-run regressors")
  }
  if (length(lr) < 2) {
    stop("'lr' must specify at least 2 variables (lagged y and at least one x)")
  }
  if (!all(lr %in% names(data))) {
    missing_vars <- lr[!lr %in% names(data)]
    stop("Long-run variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (!is.numeric(tau) || !all(tau > 0 & tau < 1)) {
    stop("All quantiles in 'tau' must be between 0 and 1 (exclusive)")
  }
  if (anyDuplicated(tau)) stop("'tau' contains duplicated quantiles")
  if (p < 1) {
    stop("'p' must be at least 1")
  }
  ## Input checks added in 1.0.4: formerly a misspecified 'lr' produced
  ## "subscript out of bounds" or ran silently with a wrong error correction term.
  fvars <- all.vars(formula)
  if (fvars[1] %in% lr) {
    stop("The dependent variable of 'formula' (", fvars[1], ") also appears in 'lr'. ",
         "'formula' must be written in first differences (d_y ~ d_x ...) and lr[1] ",
         "must be the LAGGED LEVEL of the dependent variable (for example 'L_y').")
  }
  if (any(lr %in% fvars[-1])) {
    bad <- lr[lr %in% fvars[-1]]
    stop("Variable(s) ", paste(bad, collapse = ", "), " appear both in 'formula' ",
         "(short-run differences) and in 'lr' (long-run levels). Use the ",
         "differenced variables in 'formula' and the (lagged) levels in 'lr'.")
  }
  if (anyDuplicated(lr)) {
    stop("'lr' contains duplicated variable names")
  }
  for (v in c(fvars, lr)) {
    if (!v %in% names(data)) stop("Variable '", v, "' not found in data")
    if (!is.numeric(data[[v]])) stop("Variable '", v, "' must be numeric")
  }

  ## Sort by panel and time, check duplicates and gaps, expand the formula
  ## (log(), I(), interactions) into numeric columns (1.0.5)
  prep <- .prepare_data(formula, data, id, time, allow_gaps)
  data <- prep$data
  depvar <- prep$depvar
  indepvars <- prep$indepvars
  k <- length(indepvars)

  # Get panel information
  panels <- unique(data[[id]])
  n_panels <- length(panels)

  # Long-run variables
  lr_y <- lr[1]  # Lagged y level (for ECT)
  lr_x <- lr[-1]  # Long-run x variables
  k_lr <- length(lr_x)

  # Parse q (lag orders for each x variable)
  if (length(q) == 1) {
    qlags <- rep(q, k)
  } else if (length(q) == k) {
    qlags <- q
  } else {
    stop("'q' must be a single integer or a vector of length ", k,
         " (one per column of the expanded design: ",
         paste(indepvars, collapse = ", "), ")")
  }
  if (any(qlags < 1)) stop("All lag orders in 'q' must be at least 1")

  # Number of quantiles
  ntau <- length(tau)

  # Automatic lag selection
  if (!is.null(lagsel)) {
    lagsel <- match.arg(lagsel, c("aic", "bic"))
    message("Performing ", toupper(lagsel), " lag selection (least squares)...")

    lag_result <- .select_lag_order(data, depvar, indepvars, lr, id, time,
                                     pmax, qmax, lagsel)
    p <- lag_result$p
    qlags <- lag_result$q

    message("  Selected lag order: p = ", p, ", q = ",
            paste(qlags, collapse = ","))
  }

  # Build ARDL order string
  ardl_order <- paste0("PQARDL(", p, ",", paste(qlags, collapse = ","), ")")

  n_sr <- sum(qlags)  # Total short-run coefficients
  n_ar <- max(0, p - 1)  # AR lag coefficients

  est <- if (model == "dfe") {
    .estimate_dfe(data, depvar, indepvars, lr, id, time, tau, p, qlags,
                  constant, dfe_cov = dfe_cov)
  } else {
    .estimate_mg_pmg(data, depvar, indepvars, lr, id, time, tau, p, qlags,
                     constant, pool = (model == "pmg"), pmg_weights = pmg_weights,
                     hausman_method = hausman)
  }

  rho_mg <- est$rho
  halflife_mg <- .half_life(rho_mg)

  result <- list(
    beta_mg = matrix(est$beta, nrow = 1),
    rho_mg = matrix(rho_mg, nrow = 1),
    halflife_mg = matrix(halflife_mg, nrow = 1),
    sr_mg = matrix(est$sr, nrow = 1),
    phi_mg = if (n_ar > 0) matrix(est$phi, nrow = 1) else NULL,
    beta_V = est$beta_V,
    rho_V = est$rho_V,
    beta_all = est$beta_all,
    rho_all = est$rho_all,
    halflife_all = if (is.null(est$rho_all)) NULL else
      matrix(.half_life(est$rho_all), nrow = nrow(est$rho_all)),
    phi_all = if (n_ar > 0) est$phi_all else NULL,
    sr_all = est$sr_all,
    hausman = est$hausman,
    rq_fits = est$rq_fits,
    dfe_cov = if (model == "dfe") dfe_cov else NULL,
    dfe_coef = est$dfe_coef,
    dfe_coef_V = est$dfe_coef_V,
    dfe_slopes = est$dfe_slopes,
    cross_quantile_cov = isTRUE(est$cross_quantile_cov),
    cross_quantile_cov_error = est$cross_quantile_cov_error,
    panel_cov_status = est$panel_cov_status,
    pmg_weights = if (model == "pmg") pmg_weights else NULL,
    tau = tau,
    p = p,
    q = qlags,
    model = model,
    ardl_order = ardl_order,
    n_obs = est$n_obs,
    n_dropped = est$n_dropped,
    n_panels = n_panels,
    valid_panels = est$valid_panels,
    depvar = depvar,
    indepvars = indepvars,
    lrvars = lr,
    lr_y = lr_y,
    lr_x = lr_x,
    k_lr = k_lr,
    constant = constant,
    allow_gaps = allow_gaps,
    time_step = prep$step,
    call = call
  )

  class(result) <- "xtpqardl"
  return(result)
}


#' Exact half-life ln(0.5)/ln(1 + rho), defined for -1 < rho < 0
#' @keywords internal
#' @noRd
.half_life <- function(rho) {
  out <- rho
  out[] <- NA_real_
  ok <- !is.na(rho) & rho < 0 & rho > -1
  out[ok] <- log(0.5) / log(1 + rho[ok])
  out
}


#' Sort, check the time index, expand the formula
#'
#' Returns the sorted data with one numeric column per column of the
#' expanded model matrix (names as in \code{model.matrix}) and the response,
#' plus the integer time index \code{.tidx} used to build lags by time.
#' @keywords internal
#' @noRd
.prepare_data <- function(formula, data, id, time, allow_gaps) {
  tvec <- data[[time]]
  if (!is.numeric(tvec)) {
    stop("Time variable '", time, "' must be numeric (for example the year as ",
         "an integer); found class ", class(tvec)[1])
  }
  if (anyNA(tvec) || anyNA(data[[id]])) {
    stop("'", id, "' and '", time, "' must not contain missing values (",
         sum(is.na(data[[id]])), " missing in '", id, "', ", sum(is.na(tvec)),
         " missing in '", time, "')")
  }
  data <- data[order(data[[id]], tvec), , drop = FALSE]
  dup <- duplicated(data[, c(id, time), drop = FALSE])
  if (any(dup)) {
    ex <- utils::head(which(dup), 5)
    stop(sum(dup), " duplicated (", id, ", ", time, ") pair(s), for example ",
         paste(sprintf("(%s, %s)", data[[id]][ex], data[[time]][ex]), collapse = ", "),
         ". Each panel must have one row per period.")
  }
  tt <- data[[time]]
  dts <- unlist(lapply(split(tt, data[[id]]), diff), use.names = FALSE)
  dts <- dts[is.finite(dts) & dts > 0]
  step <- if (length(dts)) min(dts) else 1
  tidx <- round((tt - min(tt)) / step)
  if (any(abs((tt - min(tt)) / step - tidx) > 1e-6)) {
    stop("The time variable '", time, "' is not on a regular grid: the unit ",
         "step was taken as ", step, " (the smallest positive difference) and ",
         "some time values are not integer multiples of it")
  }
  gaps <- lapply(split(tidx, data[[id]]), function(ti) {
    d <- diff(ti)
    w <- which(d > 1)
    if (length(w)) cbind(from = ti[w], to = ti[w + 1]) else NULL
  })
  gaps <- gaps[!vapply(gaps, is.null, TRUE)]
  if (length(gaps)) {
    ex <- utils::head(names(gaps), 5)
    txt <- vapply(ex, function(pid) {
      g <- gaps[[pid]]
      paste0("panel ", pid, ": ", paste(sprintf("%s to %s", min(tt) + g[, 1] * step,
                                                min(tt) + g[, 2] * step), collapse = ", "))
    }, "")
    msg <- paste0(length(gaps), " panel(s) have a gap in '", time, "' (", paste(txt, collapse = "; "),
                  if (length(gaps) > 5) "; ..." else "", ").")
    if (!allow_gaps) {
      stop(msg, " Lags are built by the time index; set allow_gaps = TRUE to keep ",
           "these panels (rows whose lagged period is missing are dropped).")
    }
  }
  ## Expand the formula: transformations and interactions become numeric
  ## columns named as in model.matrix; NA rows are kept (na.pass) so that the
  ## time alignment is preserved and the complete-case step happens per panel.
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(formula, data = mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (!ncol(X)) stop("'formula' must have at least one short-run regressor")
  depvar <- paste(deparse(formula[[2]]), collapse = "")
  data[[depvar]] <- y
  for (cn in colnames(X)) data[[cn]] <- as.numeric(X[, cn])
  data[[".tidx"]] <- tidx
  list(data = data, depvar = depvar, indepvars = colnames(X), step = step,
       n_gap_panels = length(gaps))
}


#' Lag a vector by j periods of the time index (NA where the period is missing)
#' @keywords internal
#' @noRd
.lag_by_time <- function(v, tidx, j) {
  if (j == 0) return(v)
  v[match(tidx - j, tidx)]
}


#' Stop when the design is exactly collinear, naming the columns
#' @keywords internal
#' @noRd
.check_collinear <- function(X, where) {
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[(qx$rank + 1):ncol(X)]]
    stop("Exactly collinear regressors in ", where, ": ",
         paste(bad, collapse = ", "), " (rank ", qx$rank, " of ", ncol(X),
         " columns). Remove or respecify these columns.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Build the per-panel PQARDL regression
#'
#' Columns: the long-run variables exactly as supplied in \code{lr} (the
#' lagged dependent level first), the lagged differences of the dependent
#' variable (p - 1 of them), and each short-run regressor with lags
#' 0, ..., q - 1; a constant is appended if requested. Lags are matched by
#' the time index (\code{.tidx}), never by row position.
#' @keywords internal
#' @noRd
.build_panel_regression <- function(panel_data, depvar, indepvars, lr,
                                    p, qlags, time, constant) {
  required_vars <- c(depvar, indepvars, lr, ".tidx")
  if (!all(required_vars %in% names(panel_data))) {
    return(NULL)
  }
  n <- nrow(panel_data)
  tidx <- panel_data[[".tidx"]]
  y <- panel_data[[depvar]]
  X_list <- list()
  col_names <- character(0)

  for (lv in lr) {
    X_list[[length(X_list) + 1]] <- panel_data[[lv]]
    col_names <- c(col_names, lv)
  }
  if (p > 1) {
    for (lag in 1:(p - 1)) {
      X_list[[length(X_list) + 1]] <- .lag_by_time(panel_data[[depvar]], tidx, lag)
      col_names <- c(col_names, paste0("L", lag, ".", depvar))
    }
  }
  for (j in seq_along(indepvars)) {
    xvar <- indepvars[j]
    for (lag in 0:(qlags[j] - 1)) {
      X_list[[length(X_list) + 1]] <- .lag_by_time(panel_data[[xvar]], tidx, lag)
      col_names <- c(col_names, if (lag == 0) xvar else paste0("L", lag, ".", xvar))
    }
  }
  if (constant) {
    X_list[[length(X_list) + 1]] <- rep(1, n)
    col_names <- c(col_names, "constant")
  }
  X <- do.call(cbind, X_list)
  colnames(X) <- col_names
  ok <- stats::complete.cases(cbind(y, X))
  y <- y[ok]
  X <- X[ok, , drop = FALSE]
  if (length(y) < ncol(X) + 1) return(NULL)
  list(y = y, X = X, n_dropped = sum(!ok))
}


#' Quantile fit with its covariance matrix (MG and PMG panels)
#' @keywords internal
#' @noRd
.qr_fit_cov <- function(y, X, tau, se = "iid") {
  fit <- tryCatch(quantreg::rq(y ~ X - 1, tau = tau, method = "br"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- stats::coef(fit)
  names(b) <- colnames(X)
  ## Covariance: try the requested estimator, then the other summary.rq
  ## estimators. summary.rq() fails when the degrees of freedom are very few
  ## (its sparsity estimate needs quantiles at tau +/- h to be estimable), so
  ## the last resort is the iid formula with a kernel sparsity estimate that
  ## only needs the residuals. The caller warns when a switch happened.
  tried <- unique(c(se, "nid", "iid", "ker"))
  V <- NULL; se_used <- NA_character_; err <- NULL
  for (s in tried) {
    V <- tryCatch(suppressWarnings(summary(fit, se = s, covariance = TRUE)$cov),
                  error = function(e) { err <<- conditionMessage(e); NULL })
    if (!is.null(V) && all(is.finite(V))) { se_used <- s; break }
    V <- NULL
  }
  if (is.null(V)) {
    V <- tryCatch(.iid_cov_fallback(fit, X, tau), error = function(e) NULL)
    if (!is.null(V)) se_used <- "iid_kernel"
  }
  if (!is.null(V)) dimnames(V) <- list(colnames(X), colnames(X))
  list(b = b, V = V, fit = fit, se_used = se_used, cov_error = err,
       n = length(y), k = ncol(X))
}


#' iid quantile-regression covariance with a kernel sparsity estimate
#'
#' \eqn{V = \tau(1-\tau) s(\tau)^2 (X'X)^{-1}}, where the sparsity
#' \eqn{s(\tau) = 1/f(0)} is estimated from the residual density at zero with
#' a Gaussian kernel (Silverman bandwidth on the non-zero residuals). Used only
#' when \code{summary.rq} cannot deliver a covariance (very few degrees of
#' freedom per panel).
#' @keywords internal
#' @noRd
.iid_cov_fallback <- function(fit, X, tau) {
  r <- stats::residuals(fit)
  r_nz <- r[abs(r) > sqrt(.Machine$double.eps)]
  if (length(r_nz) < 3L) return(NULL)
  h <- stats::bw.nrd0(r_nz)
  f0 <- mean(stats::dnorm(r_nz / h)) / h
  if (!is.finite(f0) || f0 <= 0) return(NULL)
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NULL)
  tau * (1 - tau) / f0^2 * XtX_inv
}


#' Hausman test of long-run homogeneity for one quantile
#'
#' \eqn{h = d'(V_{MG} - V_{PMG})^{-1} d}. With \code{method = "strict"} the
#' statistic is \code{NA} (with a warning) unless the difference matrix is
#' positive definite; with \code{"ginv"} a Moore-Penrose inverse over the
#' eigenvalues larger than \code{tol} times the largest absolute eigenvalue
#' is used and the degrees of freedom equal that rank.
#' @keywords internal
#' @noRd
.hausman_one <- function(d, Vd, tau, method = "strict", tol = 1e-8) {
  out <- list(tau = tau, statistic = NA_real_, df = NA_integer_, p.value = NA_real_,
              n_nonpositive = NA_integer_, method = method)
  if (anyNA(d) || anyNA(Vd)) {
    out$note <- "estimates or covariances missing"
    return(out)
  }
  Vd <- (Vd + t(Vd)) / 2
  e <- eigen(Vd, symmetric = TRUE)
  thr <- tol * max(abs(e$values))
  nonpos <- sum(e$values <= thr)
  out$n_nonpositive <- as.integer(nonpos)
  if (method == "strict") {
    if (nonpos > 0) {
      warning(sprintf(paste("Hausman test at tau = %g not defined: V(MG) - V(PMG) is not",
                            "positive definite (%d eigenvalue(s) not positive). The statistic is NA;",
                            "hausman = \"ginv\" uses a Moore-Penrose inverse instead."),
                      tau, nonpos), call. = FALSE)
      out$note <- "difference matrix not positive definite"
      return(out)
    }
    stat <- as.numeric(t(d) %*% solve(Vd, d))
    df <- length(d)
  } else {
    keep <- e$values > thr
    if (!any(keep)) {
      out$note <- "difference matrix has no positive eigenvalue"
      return(out)
    }
    if (nonpos > 0) {
      warning(sprintf(paste("Hausman test at tau = %g: %d of %d eigenvalue(s) of V(MG) - V(PMG)",
                            "discarded by the generalised inverse; the statistic can be inflated."),
                      tau, nonpos, length(d)), call. = FALSE)
    }
    Vinv <- e$vectors[, keep, drop = FALSE] %*%
      diag(1 / e$values[keep], sum(keep)) %*% t(e$vectors[, keep, drop = FALSE])
    stat <- as.numeric(t(d) %*% Vinv %*% d)
    df <- sum(keep)
  }
  out$statistic <- stat
  out$df <- as.integer(df)
  out$p.value <- stats::pchisq(stat, df, lower.tail = FALSE)
  out
}


#' Mean Group and Pooled Mean Group estimation
#'
#' Stage 1 fits the quantile ARDL panel by panel. MG averages the panel
#' coefficients, with the Pesaran-Smith covariance
#' \eqn{\sum_i (b_i - \bar b)(b_i - \bar b)' / (N(N-1))}. PMG pools the
#' long-run coefficients by minimum distance,
#' \eqn{\beta_P = (\sum_i W_i)^{-1} \sum_i W_i \beta_i} with
#' \eqn{W_i} the inverse delta-method covariance of \eqn{\beta_i}, and then
#' re-estimates every panel with \eqn{ECT = lr_y - \beta_P' lr_x} imposed;
#' the speed of adjustment and the short run are the mean group of this
#' second stage.
#' @keywords internal
#' @noRd
.estimate_mg_pmg <- function(data, depvar, indepvars, lr, id, time, tau, p,
                             qlags, constant, pool,
                             pmg_weights = "inverse_variance",
                             hausman_method = "strict") {
  lr_y <- lr[1]
  lr_x <- lr[-1]
  kx <- length(lr_x)
  ntau <- length(tau)
  n_ar <- max(0, p - 1)
  n_sr <- sum(qlags)
  panels <- unique(data[[id]])
  N <- length(panels)

  regs <- lapply(panels, function(pid) {
    .build_panel_regression(data[data[[id]] == pid, , drop = FALSE], depvar,
                            indepvars, lr, p, qlags, time, constant)
  })
  built <- !vapply(regs, is.null, TRUE)
  if (!any(built)) {
    stop("No panel has enough complete observations for the number of ",
         "coefficients (", length(lr) + n_ar + n_sr + as.integer(constant), ")")
  }
  for (i in which(built)) {
    .check_collinear(regs[[i]]$X, paste0("panel ", panels[i]))
  }
  rest_names <- setdiff(colnames(regs[[which(built)[1]]]$X), c(lr, "constant"))
  n_dropped <- sum(vapply(regs[built], function(r) r$n_dropped, numeric(1)))

  if (N < 3) {
    warning(sprintf(paste("Only %d panel(s): the Pesaran-Smith (mean group) covariance needs",
                          "several panels; with N = 2 it rejected 26 to 31 percent at nominal",
                          "5 percent in the Monte Carlo of the 1.0.4 audit, and with N = 1 it",
                          "is NA. Standard errors and Wald tests are unreliable."), N),
            call. = FALSE)
  }

  rho_all <- matrix(NA_real_, N, ntau)
  beta_all <- matrix(NA_real_, N, kx * ntau)
  rest_all <- matrix(NA_real_, N, length(rest_names) * ntau)
  beta_Vi <- vector("list", N * ntau)
  n_obs <- 0
  valid <- rep(FALSE, N)
  cov_status <- data.frame(panel = rep(as.character(panels), each = ntau),
                           tau = rep(tau, times = N), n = NA_integer_,
                           k = NA_integer_, se = NA_character_,
                           usable = FALSE, reason = NA_character_,
                           stringsAsFactors = FALSE)

  for (i in seq_len(N)) {
    rd <- regs[[i]]
    if (is.null(rd)) {
      cov_status$reason[cov_status$panel == as.character(panels[i])] <-
        "regression not built (too few observations for the number of coefficients)"
      next
    }
    for (ti in seq_len(ntau)) {
      row <- (i - 1) * ntau + ti
      f <- .qr_fit_cov(rd$y, rd$X, tau[ti])
      if (is.null(f)) { cov_status$reason[row] <- "quantile regression failed"; next }
      cov_status$n[row] <- f$n; cov_status$k[row] <- f$k; cov_status$se[row] <- f$se_used
      rho <- unname(f$b[lr_y])
      rho_all[i, ti] <- rho
      if (abs(rho) > 1e-10) {
        th <- f$b[lr_x]
        beta_all[i, (ti - 1) * kx + seq_len(kx)] <- -th / rho
        if (!is.null(f$V)) {
          # delta method for beta = -theta / rho
          G <- cbind(th / rho^2, diag(-1 / rho, kx))
          Vsub <- f$V[c(lr_y, lr_x), c(lr_y, lr_x)]
          Vb <- G %*% Vsub %*% t(G)
          beta_Vi[[row]] <- Vb
          ok <- all(is.finite(Vb)) && all(diag(Vb) > 0) &&
            !inherits(tryCatch(solve(Vb), error = function(e) e), "error")
          cov_status$usable[row] <- ok
          cov_status$reason[row] <- if (ok) "ok" else "delta-method covariance not invertible"
        } else {
          cov_status$reason[row] <- paste0("covariance not estimable (n - k = ", f$n - f$k,
                                           if (!is.null(f$cov_error)) paste0("; summary.rq: ", f$cov_error) else "", ")")
        }
      } else {
        cov_status$reason[row] <- "speed of adjustment numerically zero"
      }
      if (length(rest_names)) {
        rest_all[i, (ti - 1) * length(rest_names) + seq_along(rest_names)] <-
          f$b[rest_names]
      }
    }
    valid[i] <- TRUE
    n_obs <- n_obs + length(rd$y)
  }

  ## Warn once when the per-panel covariance estimator differed from the
  ## requested summary.rq(se = "iid") (silent in 1.0.4)
  sw <- cov_status$se[!is.na(cov_status$se) & cov_status$se != "iid"]
  if (length(sw)) {
    tb <- table(sw)
    lab <- paste(sprintf("%s (%d)", ifelse(names(tb) == "iid_kernel",
                                           "iid formula with a kernel sparsity estimate",
                                           paste0("summary.rq se = \"", names(tb), "\"")),
                         as.integer(tb)), collapse = ", ")
    warning(sprintf(paste("summary.rq(se = \"iid\") could not estimate the panel covariance for",
                          "%d of %d panel-quantile pairs; the estimator used instead: %s.",
                          "See panel_cov_status."),
                    length(sw), sum(!is.na(cov_status$se)), lab), call. = FALSE)
  }

  mg_mean <- function(M) colMeans(M, na.rm = TRUE)
  mg_V <- function(M) {
    M <- M[stats::complete.cases(M), , drop = FALSE]
    n <- nrow(M)
    if (n < 2) return(matrix(NA_real_, ncol(M), ncol(M)))
    D <- sweep(M, 2, colMeans(M))
    crossprod(D) / (n * (n - 1))
  }

  beta_mg <- mg_mean(beta_all)
  beta_V_mg <- mg_V(beta_all)
  out <- list(beta_all = beta_all, n_obs = n_obs, n_dropped = n_dropped,
              valid_panels = sum(valid), hausman = NULL,
              panel_cov_status = cov_status, rq_fits = NULL,
              dfe_coef = NULL, dfe_coef_V = NULL, dfe_slopes = NULL)

  if (!pool) {
    out$beta <- beta_mg
    out$beta_V <- beta_V_mg
    out$rho_all <- rho_all
  } else {
    beta_P <- rep(NA_real_, kx * ntau)
    # cross-quantile covariances of the pooled estimates are not estimated
    beta_VP <- matrix(NA_real_, kx * ntau, kx * ntau)
    dropped <- character(0)
    for (ti in seq_len(ntau)) {
      SW <- matrix(0, kx, kx)
      SWb <- rep(0, kx)
      np <- 0
      for (i in seq_len(N)) {
        row <- (i - 1) * ntau + ti
        bi <- beta_all[i, (ti - 1) * kx + seq_len(kx)]
        if (anyNA(bi)) next
        if (pmg_weights == "equal") {
          Wi <- diag(kx)
        } else {
          if (!isTRUE(cov_status$usable[row])) {
            dropped <- c(dropped, sprintf("%s at tau = %g (%s)", cov_status$panel[row],
                                          tau[ti], cov_status$reason[row]))
            next
          }
          Wi <- solve(beta_Vi[[row]])
        }
        SW <- SW + Wi
        SWb <- SWb + Wi %*% bi
        np <- np + 1
      }
      if (np < 2) {
        st <- cov_status[cov_status$tau == tau[ti], ]
        detail <- paste(sprintf("%s: %s", st$panel, st$reason), collapse = "; ")
        stop("PMG pooling at tau = ", tau[ti], " needs at least two panels with a usable ",
             "long-run covariance; found ", np, ". Per-panel status: ", detail, ". ",
             "With T close to the number of coefficients per panel (here n - k = ",
             paste(unique(st$n - st$k), collapse = "/"), ") the panel covariances cannot ",
             "be estimated. Options: reduce the long-run or short-run variables, use ",
             "model = \"mg\" or \"dfe\", or pmg_weights = \"equal\" (unweighted pooling, ",
             "no covariance-based pooled standard errors).", call. = FALSE)
      }
      VP <- solve(SW)
      idx <- (ti - 1) * kx + seq_len(kx)
      beta_P[idx] <- as.numeric(VP %*% SWb)
      beta_VP[idx, idx] <- if (pmg_weights == "equal") {
        ## equal weights: pooled estimate is the simple average, so its
        ## covariance is the mean-group (Pesaran-Smith) covariance
        beta_V_mg[idx, idx]
      } else VP
    }
    if (length(dropped)) {
      warning("PMG pooling dropped ", length(dropped), " panel-quantile pair(s) whose ",
              "long-run covariance is not usable: ", paste(dropped, collapse = "; "),
              ". See panel_cov_status.", call. = FALSE)
    }

    # Stage 2: impose the pooled long run, re-estimate rho and the short run
    rho_all2 <- matrix(NA_real_, N, ntau)
    rest_all <- matrix(NA_real_, N, length(rest_names) * ntau)
    for (i in seq_len(N)) {
      rd <- regs[[i]]
      if (is.null(rd)) next
      for (ti in seq_len(ntau)) {
        if (is.na(rho_all[i, ti])) next
        bP <- beta_P[(ti - 1) * kx + seq_len(kx)]
        ect <- rd$X[, lr_y] - as.numeric(rd$X[, lr_x, drop = FALSE] %*% bP)
        X2 <- cbind(ECT = ect, rd$X[, setdiff(colnames(rd$X), lr), drop = FALSE])
        f <- tryCatch(quantreg::rq.fit(X2, rd$y, tau = tau[ti], method = "br")$coefficients,
                      error = function(e) NULL)
        if (is.null(f)) next
        names(f) <- colnames(X2)
        rho_all2[i, ti] <- f[["ECT"]]
        if (length(rest_names)) {
          rest_all[i, (ti - 1) * length(rest_names) + seq_along(rest_names)] <-
            f[rest_names]
        }
      }
    }
    rho_all <- rho_all2

    # Hausman test of long-run homogeneity (MG versus PMG), quantile by
    # quantile (the cross-quantile PMG covariances are not estimated)
    H <- lapply(seq_len(ntau), function(ti) {
      idx <- (ti - 1) * kx + seq_len(kx)
      .hausman_one(beta_mg[idx] - beta_P[idx], beta_V_mg[idx, idx, drop = FALSE] -
                     beta_VP[idx, idx, drop = FALSE], tau[ti], method = hausman_method)
    })
    H <- data.frame(tau = tau,
                    statistic = vapply(H, function(h) h$statistic, numeric(1)),
                    df = vapply(H, function(h) as.integer(h$df), integer(1)),
                    p.value = vapply(H, function(h) h$p.value, numeric(1)),
                    n_nonpositive = vapply(H, function(h) as.integer(h$n_nonpositive), integer(1)),
                    method = hausman_method, stringsAsFactors = FALSE)

    out$beta <- beta_P
    out$beta_V <- beta_VP
    out$rho_all <- rho_all
    out$hausman <- H
  }

  out$rho <- mg_mean(out$rho_all)
  out$rho_V <- mg_V(out$rho_all)
  nr <- length(rest_names)
  ar_idx <- if (n_ar > 0) seq_len(n_ar) else integer(0)
  sr_idx <- n_ar + seq_len(n_sr)
  pick <- function(cols) {
    if (!length(cols)) return(NULL)
    unlist(lapply(seq_len(ntau), function(ti) (ti - 1) * nr + cols))
  }
  out$phi_all <- if (n_ar > 0) rest_all[, pick(ar_idx), drop = FALSE] else NULL
  out$sr_all <- rest_all[, pick(sr_idx), drop = FALSE]
  out$phi <- if (n_ar > 0) mg_mean(out$phi_all) else numeric(0)
  out$sr <- mg_mean(out$sr_all)
  out
}


#' Dynamic Fixed Effects estimation
#'
#' Pools all panels in one quantile regression per quantile with an
#' intercept and N - 1 panel-intercept contrasts. The joint covariance of
#' the stacked coefficient vectors comes from one estimator
#' (\code{dfe_cov}) for every block, see \code{.dfe_joint_cov}; rho_V and
#' beta_V are its (rho, theta) blocks mapped by the delta method.
#' @keywords internal
#' @noRd
.estimate_dfe <- function(data, depvar, indepvars, lr, id, time, tau, p,
                          qlags, constant, dfe_cov = "iid") {
  lr_y <- lr[1]
  lr_x <- lr[-1]
  kx <- length(lr_x)
  ntau <- length(tau)
  n_ar <- max(0, p - 1)
  n_sr <- sum(qlags)
  panels <- unique(data[[id]])

  Ys <- list(); Xs <- list(); ids <- list(); n_dropped <- 0
  for (pid in panels) {
    rd <- .build_panel_regression(data[data[[id]] == pid, , drop = FALSE], depvar,
                                  indepvars, lr, p, qlags, time,
                                  constant = FALSE)
    if (is.null(rd)) next
    Ys[[length(Ys) + 1]] <- rd$y
    Xs[[length(Xs) + 1]] <- rd$X
    ids[[length(ids) + 1]] <- rep(as.character(pid), length(rd$y))
    n_dropped <- n_dropped + rd$n_dropped
  }
  if (!length(Ys)) stop("No panel has enough complete observations for the number of ",
                        "coefficients (", length(lr) + n_ar + n_sr, ")")
  y <- unlist(Ys)
  X <- do.call(rbind, Xs)
  g <- factor(unlist(ids), levels = as.character(panels)[as.character(panels) %in% unlist(ids)])
  ## intercept + N - 1 contrasts (first panel is the reference), same column
  ## space as one dummy per panel
  D <- stats::model.matrix(~ g)
  colnames(D) <- c("(Intercept)", paste0(id, levels(g)[-1L]))
  XX <- cbind(X, D)
  .check_collinear(XX, "the pooled DFE regression")
  slopes <- colnames(X)
  rest_names <- setdiff(slopes, lr)
  pn <- ncol(XX); n <- nrow(XX)

  fits <- vector("list", ntau)
  B <- matrix(NA_real_, pn, ntau, dimnames = list(colnames(XX), paste0("tau", tau)))
  for (ti in seq_len(ntau)) {
    f <- tryCatch(quantreg::rq(y ~ XX - 1, tau = tau[ti], method = "br"),
                  error = function(e) stop("quantile regression at tau = ", tau[ti],
                                           " failed: ", conditionMessage(e), call. = FALSE))
    names(f$coefficients) <- colnames(XX)
    fits[[ti]] <- f
    B[, ti] <- f$coefficients
  }

  ## One estimator for every block of the joint covariance (1.0.5)
  jc <- .dfe_joint_cov(fits, XX, tau, type = dfe_cov)
  V <- jc$V
  bad <- jc$status[!jc$status$ok, , drop = FALSE]
  if (nrow(bad)) {
    warning(sprintf(paste("DFE covariance (dfe_cov = \"%s\") could not be estimated at tau = %s",
                          "(%s); the affected blocks of rho_V, beta_V and dfe_coef_V are NA."),
                    dfe_cov, paste(bad$tau, collapse = ", "),
                    paste(unique(bad$reason), collapse = "; ")), call. = FALSE)
  }
  if (length(jc$notes)) {
    warning(sprintf("summary.rq(se = \"%s\") reported: %s", dfe_cov,
                    paste(jc$notes, collapse = "; ")), call. = FALSE)
  }

  rho <- as.numeric(B[lr_y, ])
  beta <- as.numeric(-B[lr_x, , drop = FALSE] / rep(rho, each = kx))
  rest <- matrix(as.numeric(B[rest_names, , drop = FALSE]), nrow = 1)
  rho_V <- matrix(NA_real_, ntau, ntau)
  beta_V <- matrix(NA_real_, kx * ntau, kx * ntau)
  pos <- match(c(lr_y, lr_x), colnames(XX))
  Gs <- lapply(seq_len(ntau), function(ti) {
    r <- rho[ti]; th <- B[lr_x, ti]
    if (!is.finite(r) || abs(r) < 1e-10) return(NULL)
    cbind(th / r^2, diag(-1 / r, kx))
  })
  for (a in seq_len(ntau)) for (b in seq_len(ntau)) {
    ia <- (a - 1) * pn + pos; ib <- (b - 1) * pn + pos
    Cab <- V[ia, ib, drop = FALSE]
    rho_V[a, b] <- Cab[1, 1]
    if (!is.null(Gs[[a]]) && !is.null(Gs[[b]]) && !anyNA(Cab)) {
      beta_V[(a - 1) * kx + seq_len(kx), (b - 1) * kx + seq_len(kx)] <-
        Gs[[a]] %*% Cab %*% t(Gs[[b]])
    }
  }
  nr <- length(rest_names)
  pick <- function(cols) unlist(lapply(seq_len(ntau), function(ti) (ti - 1) * nr + cols))
  list(rho = rho, beta = beta, rho_V = rho_V, beta_V = beta_V,
       phi = if (n_ar > 0) rest[, pick(seq_len(n_ar))] else numeric(0),
       sr = rest[, pick(n_ar + seq_len(n_sr))],
       rho_all = NULL, beta_all = NULL, phi_all = NULL, sr_all = NULL,
       n_obs = n, n_dropped = n_dropped, valid_panels = length(Ys), hausman = NULL,
       rq_fits = fits, dfe_coef = B, dfe_coef_V = V, dfe_slopes = slopes,
       panel_cov_status = jc$status,
       cross_quantile_cov = all(jc$status$ok),
       cross_quantile_cov_error = if (nrow(bad)) paste(unique(bad$reason), collapse = "; ") else NULL)
}


#' Joint covariance of quantile regression fits sharing one design
#'
#' Block \eqn{(a, b)} is \eqn{(\min(\tau_a,\tau_b) - \tau_a\tau_b) H_a^{-1}
#' (X'X) H_b^{-1}}. For \code{type = "iid"} (Koenker and Bassett 1982,
#' Theorem 3.1) \eqn{H_a^{-1} = s_a (X'X)^{-1}} with the sparsity \eqn{s_a}
#' of \code{summary.rq(se = "iid")}; for \code{"nid"} and \code{"ker"}
#' \eqn{H_a^{-1}} is \code{summary.rq(se = type)$Hinv}. The diagonal blocks
#' therefore equal \code{summary.rq(se = type)$cov}.
#'
#' @param fits List of \code{rq} fits (same design \code{X}), one per tau.
#' @param X Design matrix.
#' @param tau Quantiles.
#' @param type One of "iid", "nid", "ker".
#' @return List: \code{V} (pm x pm, NA where an estimate failed),
#'   \code{status} (data frame per tau), \code{notes} (warnings issued by
#'   summary.rq).
#' @keywords internal
#' @noRd
.dfe_joint_cov <- function(fits, X, tau, type = "iid") {
  n <- nrow(X); pn <- ncol(X); m <- length(tau)
  J <- crossprod(X)
  XtXinv <- chol2inv(chol(J))
  Hinv <- vector("list", m)
  status <- data.frame(tau = tau, se = type, n = n, k = pn, ok = FALSE,
                       reason = NA_character_, stringsAsFactors = FALSE)
  notes <- character(0)
  for (a in seq_len(m)) {
    msgs <- character(0)
    sm <- tryCatch(withCallingHandlers(
      summary(fits[[a]], se = type, covariance = TRUE),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) e)
    ## quantreg's generic LP message about the sparsity or tau +/- h fits is
    ## not a covariance failure; the density messages ("non-positive fis") are kept
    msgs <- unique(msgs[!grepl("nonunique", msgs)])
    if (length(msgs)) notes <- c(notes, sprintf("tau = %g: %s", tau[a], paste(msgs, collapse = ", ")))
    if (inherits(sm, "error")) {
      status$reason[a] <- paste0("summary.rq failed: ", conditionMessage(sm))
      next
    }
    Ha <- if (type == "iid") {
      s <- 1 / sm$scale
      if (!is.finite(s) || s <= 0) NULL else s * XtXinv
    } else sm$Hinv
    if (is.null(Ha) || !all(is.finite(Ha))) {
      status$reason[a] <- "density or sparsity estimate not finite"
      next
    }
    Hinv[[a]] <- Ha
    status$ok[a] <- TRUE
    status$reason[a] <- "ok"
  }
  cn <- paste0(rep(colnames(X), m), "_tau", rep(tau, each = pn))
  V <- matrix(NA_real_, pn * m, pn * m, dimnames = list(cn, cn))
  for (a in seq_len(m)) for (b in seq_len(m)) {
    if (!status$ok[a] || !status$ok[b]) next
    w <- min(tau[a], tau[b]) - tau[a] * tau[b]
    V[(a - 1) * pn + seq_len(pn), (b - 1) * pn + seq_len(pn)] <- w * Hinv[[a]] %*% J %*% Hinv[[b]]
  }
  list(V = V, status = status, notes = notes)
}


#' Lag selection by least squares on the panel regressions
#'
#' Pooled RSS over panels; the penalty counts the coefficients of all
#' panel regressions (N times the per-panel count when all panels are used).
#' @keywords internal
#' @noRd
.select_lag_order <- function(data, depvar, indepvars, lr, id, time,
                                pmax, qmax, criterion) {
  best_ic <- Inf
  best_p <- 1
  best_q <- rep(1, length(indepvars))

  panels <- unique(data[[id]])
  table_ic <- NULL

  for (ip in 1:pmax) {
    for (iq in 1:qmax) {
      qlags <- rep(iq, length(indepvars))

      total_rss <- 0
      total_n <- 0
      total_k <- 0

      for (panel_id in panels) {
        pdata <- data[data[[id]] == panel_id, , drop = FALSE]
        reg_data <- .build_panel_regression(pdata, depvar, indepvars, lr,
                                              ip, qlags, time, constant = TRUE)

        if (!is.null(reg_data) && nrow(reg_data$X) > ncol(reg_data$X)) {
          fit <- tryCatch({
            lm.fit(reg_data$X, reg_data$y)
          }, error = function(e) NULL)

          if (!is.null(fit)) {
            resid <- fit$residuals
            total_rss <- total_rss + sum(resid^2)
            total_n <- total_n + length(resid)
            ## all estimated coefficients: one set per panel regression
            total_k <- total_k + ncol(reg_data$X)
          }
        }
      }

      if (total_n > total_k && total_rss > 0) {
        if (criterion == "bic") {
          ic <- total_n * log(total_rss / total_n) + total_k * log(total_n)
        } else {
          ic <- total_n * log(total_rss / total_n) + 2 * total_k
        }
        table_ic <- rbind(table_ic, data.frame(p = ip, q = iq, n = total_n, K = total_k,
                                               rss = total_rss, ic = ic))

        if (ic < best_ic) {
          best_ic <- ic
          best_p <- ip
          best_q <- qlags
        }
      }
    }
  }

  list(p = best_p, q = best_q, ic = best_ic, table = table_ic)
}
