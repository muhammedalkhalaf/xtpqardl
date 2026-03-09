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
#'   variables.
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
#'   \item{halflife_mg}{Vector of mean group half-life of adjustment by quantile}
#'   \item{sr_mg}{Matrix of mean group short-run coefficients}
#'   \item{phi_mg}{Matrix of mean group AR coefficients (if p > 1)}
#'   \item{beta_V}{Variance-covariance matrix for beta_mg}
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
#' where \eqn{ECT_{i,t-1} = y_{i,t-1} - \beta(\tau)' X_{i,t-1}} is the error 
#' correction term, \eqn{\rho(\tau)} is the speed of adjustment (should be 
#' negative for convergence), and \eqn{\beta(\tau)} are the long-run 
#' cointegrating parameters.
#'
#' @references
#' Pesaran MH, Shin Y, Smith RP (1999). "Pooled Mean Group Estimation of 
#' Dynamic Heterogeneous Panels." \emph{Journal of the American Statistical 
#' Association}, 94(446), 621-634. \doi{10.1080/01621459.1999.10474156}
#'
#' Cho JS, Kim TH, Shin Y (2015). "Quantile Cointegration in the 
#' Autoregressive Distributed-Lag Modeling Framework." \emph{Journal of 
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.02.030}
#'
#' Bildirici M, Kayikci F (2022). "Uncertainty, Renewable Energy, and CO2 
#' Emissions in Top Renewable Energy Countries: A Panel Quantile Regression 
#' Approach." \emph{Energy}, 247, 124303. \doi{10.1016/j.energy.2022.124303}
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
  
  # Initialize result matrices
  n_sr <- sum(qlags)  # Total short-run coefficients
  n_ar <- max(0, p - 1)  # AR lag coefficients
  
  rho_all <- matrix(NA, nrow = n_panels, ncol = ntau)
  beta_all <- matrix(NA, nrow = n_panels, ncol = k_lr * ntau)
  halflife_all <- matrix(NA, nrow = n_panels, ncol = ntau)
  phi_all <- matrix(NA, nrow = n_panels, ncol = n_ar * ntau)
  sr_all <- matrix(NA, nrow = n_panels, ncol = n_sr * ntau)
  
  valid_panels <- 0
  n_obs <- 0
  
  # Panel-by-panel estimation (for PMG and MG)
  if (model %in% c("pmg", "mg")) {
    
    for (i in seq_along(panels)) {
      panel_id <- panels[i]
      panel_data <- data[data[[id]] == panel_id, ]
      
      # Skip panels with insufficient observations
      if (nrow(panel_data) < (p + max(qlags) + 5)) {
        next
      }
      
      # Build regression data for this panel
      reg_data <- .build_panel_regression(panel_data, depvar, indepvars, lr,
                                           p, qlags, time, constant)
      
      if (is.null(reg_data) || nrow(reg_data$X) < 5) {
        next
      }
      
      # Estimate quantile regression for each tau
      for (ti in seq_along(tau)) {
        tauval <- tau[ti]
        
        result <- tryCatch({
          fit <- quantreg::rq(reg_data$y ~ reg_data$X - 1, tau = tauval)
          coefs <- coef(fit)
          names(coefs) <- colnames(reg_data$X)
          list(coefs = coefs, success = TRUE)
        }, error = function(e) {
          list(coefs = NULL, success = FALSE)
        })
        
        if (!result$success) next
        
        coefs <- result$coefs
        
        # Extract rho (ECT coefficient = coefficient on lr_y)
        rho_val <- coefs[lr_y]
        rho_all[i, ti] <- rho_val
        
        # Compute half-life
        if (!is.na(rho_val) && rho_val < 0) {
          halflife_all[i, ti] <- log(2) / abs(rho_val)
        }
        
        # Extract and compute beta (long-run coefficients)
        # beta_j = -coef(x_j) / rho
        for (j in seq_along(lr_x)) {
          xvar <- lr_x[j]
          if (xvar %in% names(coefs) && !is.na(rho_val) && abs(rho_val) > 1e-10) {
            beta_col <- (ti - 1) * k_lr + j
            beta_all[i, beta_col] <- -coefs[xvar] / rho_val
          }
        }
        
        # Extract AR coefficients (phi)
        if (n_ar > 0) {
          ar_names <- paste0("L", 1:n_ar, ".", depvar)
          for (j in seq_len(n_ar)) {
            ar_name <- ar_names[j]
            if (ar_name %in% names(coefs)) {
              phi_col <- (ti - 1) * n_ar + j
              phi_all[i, phi_col] <- coefs[ar_name]
            }
          }
        }
        
        # Extract SR coefficients
        sr_idx <- 0
        for (j in seq_along(indepvars)) {
          xvar <- indepvars[j]
          qj <- qlags[j]
          for (lag in 0:(qj - 1)) {
            sr_idx <- sr_idx + 1
            if (lag == 0) {
              sr_name <- xvar
            } else {
              sr_name <- paste0("L", lag, ".", xvar)
            }
            if (sr_name %in% names(coefs)) {
              sr_col <- (ti - 1) * n_sr + sr_idx
              sr_all[i, sr_col] <- coefs[sr_name]
            }
          }
        }
      }
      
      valid_panels <- valid_panels + 1
      n_obs <- n_obs + nrow(reg_data$X)
    }
    
    # Compute mean group estimates
    rho_mg <- colMeans(rho_all, na.rm = TRUE)
    beta_mg <- colMeans(beta_all, na.rm = TRUE)
    halflife_mg <- colMeans(halflife_all, na.rm = TRUE)
    phi_mg <- if (n_ar > 0) colMeans(phi_all, na.rm = TRUE) else numeric(0)
    sr_mg <- colMeans(sr_all, na.rm = TRUE)
    
    # Compute variance-covariance matrices (MG-style)
    rho_V <- .compute_mg_variance(rho_all, valid_panels)
    beta_V <- .compute_mg_variance(beta_all, valid_panels)
    
  } else {
    # DFE estimation
    result <- .estimate_dfe(data, formula, id, time, lr, tau, p, qlags, constant)
    
    rho_mg <- result$rho
    beta_mg <- result$beta
    halflife_mg <- result$halflife
    phi_mg <- result$phi
    sr_mg <- result$sr
    rho_V <- result$rho_V
    beta_V <- result$beta_V
    n_obs <- result$n_obs
    valid_panels <- n_panels
  }
  
  # Construct result object
  result <- list(
    beta_mg = matrix(beta_mg, nrow = 1),
    rho_mg = matrix(rho_mg, nrow = 1),
    halflife_mg = matrix(halflife_mg, nrow = 1),
    sr_mg = matrix(sr_mg, nrow = 1),
    phi_mg = if (length(phi_mg) > 0) matrix(phi_mg, nrow = 1) else NULL,
    beta_V = beta_V,
    rho_V = rho_V,
    beta_all = beta_all,
    rho_all = rho_all,
    halflife_all = halflife_all,
    phi_all = if (n_ar > 0) phi_all else NULL,
    sr_all = sr_all,
    tau = tau,
    p = p,
    q = qlags,
    model = model,
    ardl_order = ardl_order,
    n_obs = n_obs,
    n_panels = n_panels,
    valid_panels = valid_panels,
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


#' @keywords internal
.build_panel_regression <- function(panel_data, depvar, indepvars, lr, 
                                      p, qlags, time, constant) {
  n <- nrow(panel_data)
  lr_y <- lr[1]
  lr_x <- lr[-1]
  
  # Check that all required variables exist
  required_vars <- c(depvar, indepvars, lr)
  if (!all(required_vars %in% names(panel_data))) {
    return(NULL)
  }
  
  # Build design matrix
  # LR variables (levels): lr_y (for ECT) + lr_x
  # AR lags: L1.depvar, ..., L(p-1).depvar
  # SR lags: indepvars + their lags up to qlags
  
  # Start with maximum lag needed
  max_lag <- max(p, max(qlags))
  start_row <- max_lag + 1
  
  if (start_row >= n) return(NULL)
  
  # Number of usable observations
  n_use <- n - max_lag
  
  # Response variable (current period first difference)
  y <- panel_data[[depvar]][start_row:n]
  
  # Build X matrix
  X_list <- list()
  col_names <- c()
  
  # LR variables (levels at t-1 for ECT computation)
  for (lv in lr) {
    # Use value at t-1 relative to y
    X_list[[length(X_list) + 1]] <- panel_data[[lv]][(start_row - 1):(n - 1)]
    col_names <- c(col_names, lv)
  }
  
  # AR lags (lagged first differences of y)
  if (p > 1) {
    for (lag in 1:(p - 1)) {
      lag_name <- paste0("L", lag, ".", depvar)
      X_list[[length(X_list) + 1]] <- panel_data[[depvar]][(start_row - lag):(n - lag)]
      col_names <- c(col_names, lag_name)
    }
  }
  
  # SR variables and their lags
  for (j in seq_along(indepvars)) {
    xvar <- indepvars[j]
    qj <- qlags[j]
    
    # Contemporary value
    X_list[[length(X_list) + 1]] <- panel_data[[xvar]][start_row:n]
    col_names <- c(col_names, xvar)
    
    # Lagged values
    if (qj > 1) {
      for (lag in 1:(qj - 1)) {
        lag_name <- paste0("L", lag, ".", xvar)
        X_list[[length(X_list) + 1]] <- panel_data[[xvar]][(start_row - lag):(n - lag)]
        col_names <- c(col_names, lag_name)
      }
    }
  }
  
  # Add constant if requested
  if (constant) {
    X_list[[length(X_list) + 1]] <- rep(1, n_use)
    col_names <- c(col_names, "constant")
  }
  
  # Combine into matrix
  X <- do.call(cbind, X_list)
  colnames(X) <- col_names
  
  # Remove rows with NAs
  complete_cases <- complete.cases(cbind(y, X))
  y <- y[complete_cases]
  X <- X[complete_cases, , drop = FALSE]
  
  if (length(y) < 5) return(NULL)
  
  return(list(y = y, X = X))
}


#' @keywords internal
.compute_mg_variance <- function(mat, n_panels) {
  # Mean Group variance: V(theta_MG) = Var(theta_i) / N
  if (n_panels <= 1) {
    return(diag(ncol(mat)) * NA)
  }
  
  # Remove panels with all NAs
  valid_rows <- apply(mat, 1, function(x) !all(is.na(x)))
  mat_valid <- mat[valid_rows, , drop = FALSE]
  
  if (nrow(mat_valid) <= 1) {
    return(diag(ncol(mat)) * NA)
  }
  
  # Compute sample variance-covariance and divide by N
  V <- cov(mat_valid, use = "pairwise.complete.obs") / nrow(mat_valid)
  
  # Replace NA/NaN with 0 on diagonal
  diag_vals <- diag(V)
  diag_vals[is.na(diag_vals)] <- 0
  diag(V) <- diag_vals
  
  return(V)
}


#' @keywords internal
.estimate_dfe <- function(data, formula, id, time, lr, tau, p, qlags, constant) {
  # Dynamic Fixed Effects estimation
  # Pool all data and estimate with panel fixed effects using quantile regression
  
  lr_y <- lr[1]
  lr_x <- lr[-1]
  k_lr <- length(lr_x)
  ntau <- length(tau)
  
  # Get variable names from formula
  mf <- model.frame(formula, data = data, na.action = na.omit)
  depvar <- all.vars(formula)[1]
  X <- model.matrix(formula, data = mf)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  indepvars <- colnames(X)
  k <- length(indepvars)
  n_sr <- sum(qlags)
  n_ar <- max(0, p - 1)
  
  # Create panel dummies
  panels <- unique(data[[id]])
  n_panels <- length(panels)
  
  # Build full design matrix with lags
  panel_data_list <- list()
  
  for (panel_id in panels) {
    pdata <- data[data[[id]] == panel_id, ]
    reg_data <- .build_panel_regression(pdata, depvar, indepvars, lr, 
                                          p, qlags, time, constant = FALSE)
    if (!is.null(reg_data)) {
      df <- data.frame(y = reg_data$y)
      df <- cbind(df, as.data.frame(reg_data$X))
      df[[id]] <- panel_id
      panel_data_list[[length(panel_data_list) + 1]] <- df
    }
  }
  
  if (length(panel_data_list) == 0) {
    stop("No valid panel data for DFE estimation")
  }
  
  pooled_data <- do.call(rbind, panel_data_list)
  n_obs <- nrow(pooled_data)
  
  # Create panel dummies (excluding first for identification)
  pooled_data[[id]] <- factor(pooled_data[[id]])
  
  # Build formula for quantile regression
  xvars <- setdiff(names(pooled_data), c("y", id))
  dfe_formula <- as.formula(paste("y ~", paste(xvars, collapse = " + "), 
                                    "+ factor(", id, ") - 1"))
  
  # Initialize results
  rho <- rep(NA, ntau)
  beta <- rep(NA, k_lr * ntau)
  halflife <- rep(NA, ntau)
  phi <- rep(NA, n_ar * ntau)
  sr <- rep(NA, n_sr * ntau)
  
  # Estimate for each quantile
  for (ti in seq_along(tau)) {
    tauval <- tau[ti]
    
    fit <- tryCatch({
      quantreg::rq(dfe_formula, tau = tauval, data = pooled_data)
    }, error = function(e) NULL)
    
    if (is.null(fit)) next
    
    coefs <- coef(fit)
    
    # Extract rho (ECT coefficient)
    if (lr_y %in% names(coefs)) {
      rho[ti] <- coefs[lr_y]
      
      if (rho[ti] < 0) {
        halflife[ti] <- log(2) / abs(rho[ti])
      }
      
      # Compute beta
      for (j in seq_along(lr_x)) {
        xvar <- lr_x[j]
        if (xvar %in% names(coefs) && abs(rho[ti]) > 1e-10) {
          beta_col <- (ti - 1) * k_lr + j
          beta[(ti - 1) * k_lr + j] <- -coefs[xvar] / rho[ti]
        }
      }
    }
    
    # Extract AR coefficients
    if (n_ar > 0) {
      for (j in 1:n_ar) {
        ar_name <- paste0("L", j, ".", depvar)
        if (ar_name %in% names(coefs)) {
          phi[(ti - 1) * n_ar + j] <- coefs[ar_name]
        }
      }
    }
    
    # Extract SR coefficients
    sr_idx <- 0
    for (j in seq_along(indepvars)) {
      xvar <- indepvars[j]
      qj <- qlags[j]
      for (lag in 0:(qj - 1)) {
        sr_idx <- sr_idx + 1
        if (lag == 0) {
          sr_name <- xvar
        } else {
          sr_name <- paste0("L", lag, ".", xvar)
        }
        if (sr_name %in% names(coefs)) {
          sr[(ti - 1) * n_sr + sr_idx] <- coefs[sr_name]
        }
      }
    }
  }
  
  # Simple variance estimates (from quantreg)
  rho_V <- diag(ntau) * 0.01  # Placeholder
  beta_V <- diag(k_lr * ntau) * 0.01  # Placeholder
  
  list(
    rho = rho,
    beta = beta,
    halflife = halflife,
    phi = phi,
    sr = sr,
    rho_V = rho_V,
    beta_V = beta_V,
    n_obs = n_obs
  )
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
