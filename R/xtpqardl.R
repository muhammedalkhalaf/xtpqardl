#' Panel Quantile Autoregressive Distributed Lag Model
#'
#' Estimate Panel Quantile ARDL (PQARDL) models that combine panel ARDL 
#' methodology with quantile regression. Supports Pooled Mean Group (PMG), 
#' Mean Group (MG), and Dynamic Fixed Effects (DFE) estimators.
#'
#' @param formula A formula specifying the model. The response variable should 
#'   be in first differences (e.g., \code{d.y ~ d.x1 + d.x2}), where short-run 
#'   dynamics are estimated.
#' @param data A data frame containing panel data with variables specified in 
#'   the formula and \code{lr} argument.
#' @param id Character string specifying the panel (cross-section) identifier 
#'   variable name.
#' @param time Character string specifying the time variable name.
#' @param lr Character vector of long-run level variable names. The first 
#'   element should be the lagged dependent variable level (for the error 
#'   correction term), and remaining elements are the long-run explanatory 
#'   variables. The variables enter the regression exactly as supplied (no
#'   further lag is applied); supply lagged regressors, for example
#'   \code{L_x1}, if the error correction term should use \eqn{x_{t-1}}.
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
#'   or \code{"aic"}, optimal lag orders are selected using the specified 
#'   criterion. Default is \code{NULL} (no automatic selection).
#' @param pmax Maximum p to consider in lag selection. Default is 4.
#' @param qmax Maximum q to consider in lag selection. Default is 4.
#' @param constant Logical. Include a constant term? Default is \code{TRUE}.
#'
#' @return An object of class \code{"xtpqardl"} containing:
#' \describe{
#'   \item{beta_mg}{Matrix of mean group long-run coefficients across quantiles}
#'   \item{rho_mg}{Vector of mean group ECT speed of adjustment by quantile}
#'   \item{halflife_mg}{Half-life \eqn{\ln(0.5)/\ln(1+\rho)} of the mean
#'     group speed of adjustment, by quantile}
#'   \item{sr_mg}{Matrix of mean group short-run coefficients}
#'   \item{phi_mg}{Matrix of mean group AR coefficients (if p > 1)}
#'   \item{hausman}{For \code{model = "pmg"}, the Hausman test of long-run
#'     homogeneity (MG against PMG): statistic, df and p-value}
#'   \item{beta_V}{Variance-covariance matrix for beta_mg (cross-quantile
#'     blocks are \code{NA} for PMG and DFE, where they are not estimated)}
#'   \item{rho_V}{Variance-covariance matrix for rho_mg}
#'   \item{beta_all}{Matrix of per-panel long-run coefficients}
#'   \item{rho_all}{Matrix of per-panel ECT coefficients}
#'   \item{halflife_all}{Matrix of per-panel half-life values}
#'   \item{tau}{Vector of estimated quantiles}
#'   \item{p}{AR lag order used}
#'   \item{q}{Distributed lag order(s) used}
#'   \item{model}{Estimation method used}
#'   \item{n_obs}{Total number of observations}
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
#' the long-run parameters, \eqn{\beta = -\theta/\rho}.
#'
#' \code{model = "mg"} averages the panel estimates (Pesaran and Smith
#' covariance). \code{model = "pmg"} pools the long run by minimum distance,
#' weighting each panel by the inverse of its delta-method covariance, then
#' re-estimates every panel with the pooled long run imposed; the speed of
#' adjustment and the short run are the mean group of that second stage, and
#' a Hausman test of long-run homogeneity (MG against PMG) is returned.
#' \code{model = "dfe"} pools all panels with panel intercepts; its standard
#' errors come from the Powell kernel sandwich of the quantile regression.
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
#' Koenker R, Bassett G (1978). "Regression Quantiles." \emph{Econometrica}, 
#' 46(1), 33-50. \doi{10.2307/1913643}
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
#'   pnorm pchisq formula terms complete.cases as.formula coef lm.fit
#' @importFrom quantreg rq
#' @export
xtpqardl <- function(formula, data, id, time, lr,
                     tau = c(0.25, 0.50, 0.75),
                     p = 1, q = 1,
                     model = c("pmg", "mg", "dfe"),
                     lagsel = NULL,
                     pmax = 4, qmax = 4,
                     constant = TRUE) {
  

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
  if (length(lr) < 2) {
    stop("'lr' must specify at least 2 variables (lagged y and at least one x)")
  }
  if (!all(lr %in% names(data))) {
    missing_vars <- lr[!lr %in% names(data)]
    stop("Long-run variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (!all(tau > 0 & tau < 1)) {
    stop("All quantiles in 'tau' must be between 0 and 1 (exclusive)")
  }
  if (p < 1) {
    stop("'p' must be at least 1")
  }
  
  # Sort data by panel and time
  data <- data[order(data[[id]], data[[time]]), ]
  
  # Get panel information
  panels <- unique(data[[id]])
  n_panels <- length(panels)
  
  # Parse formula
  mf <- model.frame(formula, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)
  
  # Remove intercept from X if present (will add back separately)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  
  depvar <- all.vars(formula)[1]
  indepvars <- colnames(X)
  k <- length(indepvars)
  
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
    stop("'q' must be a single integer or a vector of length ", k)
  }
  
  # Number of quantiles
  ntau <- length(tau)
  
  # Automatic lag selection
  if (!is.null(lagsel)) {
    lagsel <- match.arg(lagsel, c("aic", "bic"))
    message("Performing ", toupper(lagsel), " lag selection...")
    
    lag_result <- .select_lag_order(data, depvar, indepvars, lr, id, time,
                                     pmax, qmax, lagsel)
    p <- lag_result$p
    qlags <- lag_result$q
    
    message("  Optimal lag order: p = ", p, ", q = ", 
            paste(qlags, collapse = ","))
  }
  
  # Build ARDL order string
  ardl_order <- paste0("PQARDL(", p, ",", paste(qlags, collapse = ","), ")")
  
  n_sr <- sum(qlags)  # Total short-run coefficients
  n_ar <- max(0, p - 1)  # AR lag coefficients

  est <- if (model == "dfe") {
    .estimate_dfe(data, depvar, indepvars, lr, id, time, tau, p, qlags,
                  constant)
  } else {
    .estimate_mg_pmg(data, depvar, indepvars, lr, id, time, tau, p, qlags,
                     constant, pool = (model == "pmg"))
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
    tau = tau,
    p = p,
    q = qlags,
    model = model,
    ardl_order = ardl_order,
    n_obs = est$n_obs,
    n_panels = n_panels,
    valid_panels = est$valid_panels,
    depvar = depvar,
    indepvars = indepvars,
    lrvars = lr,
    lr_y = lr_y,
    lr_x = lr_x,
    k_lr = k_lr,
    constant = constant,
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


#' Lag a vector by k periods (NA padded)
#' @keywords internal
#' @noRd
.lagk <- function(v, k) {
  if (k == 0) return(v)
  n <- length(v)
  if (k >= n) return(rep(NA_real_, n))
  c(rep(NA_real_, k), v[1:(n - k)])
}


#' Build the per-panel PQARDL regression
#'
#' Columns: the long-run variables exactly as supplied in \code{lr} (the
#' lagged dependent level first), the lagged differences of the dependent
#' variable (p - 1 of them), and each short-run regressor with lags
#' 0, ..., q - 1; a constant is appended if requested.
#' @keywords internal
#' @noRd
.build_panel_regression <- function(panel_data, depvar, indepvars, lr,
                                    p, qlags, time, constant) {
  required_vars <- c(depvar, indepvars, lr)
  if (!all(required_vars %in% names(panel_data))) {
    return(NULL)
  }
  n <- nrow(panel_data)
  y <- panel_data[[depvar]]
  X_list <- list()
  col_names <- character(0)

  for (lv in lr) {
    X_list[[length(X_list) + 1]] <- panel_data[[lv]]
    col_names <- c(col_names, lv)
  }
  if (p > 1) {
    for (lag in 1:(p - 1)) {
      X_list[[length(X_list) + 1]] <- .lagk(panel_data[[depvar]], lag)
      col_names <- c(col_names, paste0("L", lag, ".", depvar))
    }
  }
  for (j in seq_along(indepvars)) {
    xvar <- indepvars[j]
    for (lag in 0:(qlags[j] - 1)) {
      X_list[[length(X_list) + 1]] <- .lagk(panel_data[[xvar]], lag)
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
  list(y = y, X = X)
}


#' Quantile fit with its iid covariance matrix
#' @keywords internal
#' @noRd
.qr_fit_cov <- function(y, X, tau, se = "iid") {
  fit <- tryCatch(quantreg::rq(y ~ X - 1, tau = tau, method = "br"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- stats::coef(fit)
  names(b) <- colnames(X)
  V <- tryCatch(suppressWarnings(
         summary(fit, se = se, covariance = TRUE)$cov),
       error = function(e) NULL)
  if (!is.null(V)) dimnames(V) <- list(colnames(X), colnames(X))
  list(b = b, V = V)
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
                             qlags, constant, pool) {
  lr_y <- lr[1]
  lr_x <- lr[-1]
  kx <- length(lr_x)
  ntau <- length(tau)
  n_ar <- max(0, p - 1)
  n_sr <- sum(qlags)
  panels <- unique(data[[id]])
  N <- length(panels)

  regs <- lapply(panels, function(pid) {
    .build_panel_regression(data[data[[id]] == pid, ], depvar, indepvars,
                            lr, p, qlags, time, constant)
  })
  rest_names <- setdiff(colnames(regs[[which(!vapply(regs, is.null, TRUE))[1]]]$X),
                        c(lr, "constant"))

  rho_all <- matrix(NA_real_, N, ntau)
  beta_all <- matrix(NA_real_, N, kx * ntau)
  rest_all <- matrix(NA_real_, N, length(rest_names) * ntau)
  beta_Vi <- vector("list", N * ntau)
  n_obs <- 0
  valid <- rep(FALSE, N)

  for (i in seq_len(N)) {
    rd <- regs[[i]]
    if (is.null(rd)) next
    for (ti in seq_len(ntau)) {
      f <- .qr_fit_cov(rd$y, rd$X, tau[ti])
      if (is.null(f)) next
      rho <- unname(f$b[lr_y])
      rho_all[i, ti] <- rho
      if (abs(rho) > 1e-10) {
        th <- f$b[lr_x]
        beta_all[i, (ti - 1) * kx + seq_len(kx)] <- -th / rho
        if (!is.null(f$V)) {
          # delta method for beta = -theta / rho
          G <- cbind(th / rho^2, diag(-1 / rho, kx))
          Vsub <- f$V[c(lr_y, lr_x), c(lr_y, lr_x)]
          beta_Vi[[(i - 1) * ntau + ti]] <- G %*% Vsub %*% t(G)
        }
      }
      if (length(rest_names)) {
        rest_all[i, (ti - 1) * length(rest_names) + seq_along(rest_names)] <-
          f$b[rest_names]
      }
    }
    valid[i] <- TRUE
    n_obs <- n_obs + length(rd$y)
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
  out <- list(beta_all = beta_all, n_obs = n_obs,
              valid_panels = sum(valid), hausman = NULL)

  if (!pool) {
    out$beta <- beta_mg
    out$beta_V <- beta_V_mg
    out$rho_all <- rho_all
  } else {
    beta_P <- rep(NA_real_, kx * ntau)
    # cross-quantile covariances of the pooled estimates are not estimated
    beta_VP <- matrix(NA_real_, kx * ntau, kx * ntau)
    for (ti in seq_len(ntau)) {
      SW <- matrix(0, kx, kx)
      SWb <- rep(0, kx)
      np <- 0
      for (i in seq_len(N)) {
        Vi <- beta_Vi[[(i - 1) * ntau + ti]]
        bi <- beta_all[i, (ti - 1) * kx + seq_len(kx)]
        if (is.null(Vi) || anyNA(bi) || anyNA(Vi) || any(diag(Vi) <= 0)) next
        Wi <- tryCatch(solve(Vi), error = function(e) NULL)
        if (is.null(Wi)) next
        SW <- SW + Wi
        SWb <- SWb + Wi %*% bi
        np <- np + 1
      }
      if (np < 2) stop("PMG pooling needs at least two panels with a usable long-run covariance; use model = \"mg\".")
      VP <- solve(SW)
      idx <- (ti - 1) * kx + seq_len(kx)
      beta_P[idx] <- as.numeric(VP %*% SWb)
      beta_VP[idx, idx] <- VP
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

    # Hausman test of long-run homogeneity (MG versus PMG)
    d <- beta_mg - beta_P
    VP0 <- beta_VP
    VP0[is.na(VP0)] <- 0
    Vd <- beta_V_mg - VP0
    H <- tryCatch({
      e <- eigen((Vd + t(Vd)) / 2, symmetric = TRUE)
      keep <- e$values > 1e-12 * max(abs(e$values))
      Vinv <- e$vectors[, keep, drop = FALSE] %*%
        diag(1 / e$values[keep], sum(keep)) %*% t(e$vectors[, keep, drop = FALSE])
      stat <- as.numeric(t(d) %*% Vinv %*% d)
      list(statistic = stat, df = sum(keep),
           p.value = stats::pchisq(stat, sum(keep), lower.tail = FALSE))
    }, error = function(e) NULL)

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
#' Pools all panels in one quantile regression with panel-specific
#' intercepts; the covariance of (rho, beta) is the delta-method transform
#' of the Powell kernel sandwich covariance of the quantile regression.
#' @keywords internal
#' @noRd
.estimate_dfe <- function(data, depvar, indepvars, lr, id, time, tau, p,
                          qlags, constant) {
  lr_y <- lr[1]
  lr_x <- lr[-1]
  kx <- length(lr_x)
  ntau <- length(tau)
  n_ar <- max(0, p - 1)
  n_sr <- sum(qlags)
  panels <- unique(data[[id]])

  Ys <- list(); Xs <- list(); ids <- list()
  for (pid in panels) {
    rd <- .build_panel_regression(data[data[[id]] == pid, ], depvar,
                                  indepvars, lr, p, qlags, time,
                                  constant = FALSE)
    if (is.null(rd)) next
    Ys[[length(Ys) + 1]] <- rd$y
    Xs[[length(Xs) + 1]] <- rd$X
    ids[[length(ids) + 1]] <- rep(as.character(pid), length(rd$y))
  }
  if (!length(Ys)) stop("No valid panel data for DFE estimation")
  y <- unlist(Ys)
  X <- do.call(rbind, Xs)
  g <- factor(unlist(ids))
  D <- stats::model.matrix(~ g - 1)
  XD <- cbind(X, D)
  rest_names <- setdiff(colnames(X), lr)

  rho <- rep(NA_real_, ntau)
  beta <- rep(NA_real_, kx * ntau)
  rest <- matrix(NA_real_, 1, length(rest_names) * ntau)
  rho_V <- matrix(NA_real_, ntau, ntau)
  # cross-quantile covariances are not estimated for DFE
  beta_V <- matrix(NA_real_, kx * ntau, kx * ntau)
  for (ti in seq_len(ntau)) {
    f <- .qr_fit_cov(y, XD, tau[ti], se = "ker")
    if (is.null(f)) next
    r <- unname(f$b[lr_y])
    rho[ti] <- r
    th <- f$b[lr_x]
    idx <- (ti - 1) * kx + seq_len(kx)
    beta[idx] <- -th / r
    rest[1, (ti - 1) * length(rest_names) + seq_along(rest_names)] <- f$b[rest_names]
    if (!is.null(f$V)) {
      rho_V[ti, ti] <- f$V[lr_y, lr_y]
      G <- cbind(th / r^2, diag(-1 / r, kx))
      beta_V[idx, idx] <- G %*% f$V[c(lr_y, lr_x), c(lr_y, lr_x)] %*% t(G)
    }
  }
  nr <- length(rest_names)
  pick <- function(cols) unlist(lapply(seq_len(ntau), function(ti) (ti - 1) * nr + cols))
  list(rho = rho, beta = beta, rho_V = rho_V, beta_V = beta_V,
       phi = if (n_ar > 0) rest[, pick(seq_len(n_ar))] else numeric(0),
       sr = rest[, pick(n_ar + seq_len(n_sr))],
       rho_all = NULL, beta_all = NULL, phi_all = NULL, sr_all = NULL,
       n_obs = length(y), valid_panels = length(Ys), hausman = NULL)
}


#' @keywords internal
.select_lag_order <- function(data, depvar, indepvars, lr, id, time,
                                pmax, qmax, criterion) {
  # Lag selection using BIC or AIC
  best_ic <- Inf
  best_p <- 1
  best_q <- rep(1, length(indepvars))
  
  panels <- unique(data[[id]])
  
  for (ip in 1:pmax) {
    for (iq in 1:qmax) {
      qlags <- rep(iq, length(indepvars))
      
      # Pool data and fit OLS
      total_rss <- 0
      total_n <- 0
      total_k <- 0
      
      for (panel_id in panels) {
        pdata <- data[data[[id]] == panel_id, ]
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
            total_k <- ncol(reg_data$X)
          }
        }
      }
      
      if (total_n > total_k) {
        if (criterion == "bic") {
          ic <- total_n * log(total_rss / total_n) + total_k * log(total_n)
        } else {
          ic <- total_n * log(total_rss / total_n) + 2 * total_k
        }
        
        if (ic < best_ic) {
          best_ic <- ic
          best_p <- ip
          best_q <- qlags
        }
      }
    }
  }
  
  list(p = best_p, q = best_q, ic = best_ic)
}
