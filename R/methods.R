#' Print Method for xtpqardl Objects
#'
#' @param x An object of class \code{"xtpqardl"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.xtpqardl <- function(x, ...) {
  cat("\nPanel Quantile ARDL Estimation\n")
  cat(rep("=", 50), "\n", sep = "")
  cat("Model:           ", toupper(x$model), "\n")
  cat("ARDL order:      ", x$ardl_order, "\n")

  cat("Quantiles:       ", paste(sprintf("%.2f", x$tau), collapse = ", "), "\n")
  cat("Panels:          ", x$n_panels, " (", x$valid_panels, " valid)\n", sep = "")
  cat("Observations:    ", x$n_obs, "\n")
  cat("Dependent var:   ", x$depvar, "\n")
  cat("LR variables:    ", paste(x$lr_x, collapse = ", "), "\n")
  cat("ECT variable:    ", x$lr_y, "\n")
  cat(rep("=", 50), "\n\n", sep = "")
  
  cat("Use summary() for detailed results\n")
  cat("Use wald_test() to test parameter equality across quantiles\n\n")
  
  invisible(x)
}


#' Summary Method for xtpqardl Objects
#'
#' Produces a detailed summary of Panel Quantile ARDL estimation results,
#' including long-run coefficients, ECT speed of adjustment, half-life
#' of adjustment, and short-run parameters by quantile.
#'
#' @param object An object of class \code{"xtpqardl"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.xtpqardl"} containing formatted
#'   tables of results.
#'
#' @examples
#' \donttest{
#' data(pqardl_sample)
#' fit <- xtpqardl(
#'   formula = d_y ~ d_x1 + d_x2,
#'   data = pqardl_sample,
#'   id = "country",
#'   time = "year",
#'   lr = c("L_y", "x1", "x2"),
#'   tau = c(0.25, 0.50, 0.75)
#' )
#' summary(fit)
#' }
#'
#' @export
summary.xtpqardl <- function(object, ...) {
  
  x <- object
  ntau <- length(x$tau)
  k_lr <- x$k_lr
  
  # Build long-run coefficient table
  lr_table <- data.frame(
    variable = character(),
    tau = numeric(),
    estimate = numeric(),
    std_error = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    signif = character(),
    stringsAsFactors = FALSE
  )
  
  for (ti in seq_along(x$tau)) {
    for (j in seq_along(x$lr_x)) {
      col_idx <- (ti - 1) * k_lr + j
      est <- x$beta_mg[1, col_idx]
      
      se <- NA
      z <- NA
      pval <- NA
      stars <- ""
      
      if (!is.null(x$beta_V) && col_idx <= ncol(x$beta_V)) {
        var_val <- x$beta_V[col_idx, col_idx]
        if (!is.na(var_val) && var_val > 0) {
          se <- sqrt(var_val)
          z <- est / se
          pval <- 2 * (1 - pnorm(abs(z)))
          
          if (!is.na(pval)) {
            if (pval < 0.01) stars <- "***"
            else if (pval < 0.05) stars <- "**"
            else if (pval < 0.10) stars <- "*"
          }
        }
      }
      
      lr_table <- rbind(lr_table, data.frame(
        variable = x$lr_x[j],
        tau = x$tau[ti],
        estimate = est,
        std_error = se,
        z_value = z,
        p_value = pval,
        signif = stars,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Build ECT table
  ect_table <- data.frame(
    tau = x$tau,
    rho = as.vector(x$rho_mg),
    std_error = NA,
    z_value = NA,
    p_value = NA,
    halflife = as.vector(x$halflife_mg),
    status = character(ntau),
    stringsAsFactors = FALSE
  )
  
  for (ti in seq_along(x$tau)) {
    rho_val <- x$rho_mg[1, ti]
    
    if (!is.null(x$rho_V) && ti <= ncol(x$rho_V)) {
      var_val <- x$rho_V[ti, ti]
      if (!is.na(var_val) && var_val > 0) {
        ect_table$std_error[ti] <- sqrt(var_val)
        ect_table$z_value[ti] <- rho_val / ect_table$std_error[ti]
        ect_table$p_value[ti] <- 2 * (1 - pnorm(abs(ect_table$z_value[ti])))
      }
    }
    
    if (!is.na(rho_val)) {
      if (rho_val < -0.5) {
        ect_table$status[ti] <- "Strong"
      } else if (rho_val < -0.1) {
        ect_table$status[ti] <- "Moderate"
      } else if (rho_val < 0) {
        ect_table$status[ti] <- "Weak"
      } else {
        ect_table$status[ti] <- "No conv."
      }
    } else {
      ect_table$status[ti] <- "N/A"
    }
  }
  
  result <- list(
    call = x$call,
    model = x$model,
    ardl_order = x$ardl_order,
    tau = x$tau,
    n_obs = x$n_obs,
    n_panels = x$n_panels,
    valid_panels = x$valid_panels,
    depvar = x$depvar,
    lr_x = x$lr_x,
    lr_y = x$lr_y,
    lr_table = lr_table,
    ect_table = ect_table,
    beta_mg = x$beta_mg,
    rho_mg = x$rho_mg
  )
  
  class(result) <- "summary.xtpqardl"
  return(result)
}


#' Print Method for summary.xtpqardl Objects
#'
#' @param x An object of class \code{"summary.xtpqardl"}.
#' @param digits Number of significant digits to display. Default is 4.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.xtpqardl <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  XTPQARDL - Panel Quantile Autoregressive Distributed Lag Model\n")
  cat(rep("=", 78), "\n", sep = "")
  
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nModel:           ", toupper(x$model), "\n")
  cat("ARDL order:      ", x$ardl_order, "\n")
  cat("Dependent var:   ", x$depvar, "\n")
  cat("LR variables:    ", paste(x$lr_x, collapse = ", "), "\n")
  cat("ECT variable:    ", x$lr_y, "\n")
  cat("Panels:          ", x$n_panels, " (", x$valid_panels, " valid)\n", sep = "")
  cat("Observations:    ", x$n_obs, "\n")
  cat("Quantiles:       ", paste(sprintf("%.2f", x$tau), collapse = ", "), "\n")
  
  # Table 1: Long-Run Coefficients
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  Table 1: Long-Run Cointegrating Parameters beta(tau)\n")
  cat("  beta_j(tau) = -coef(x_j) / rho(tau)\n")
  cat(rep("-", 78), "\n", sep = "")
  
  cat(sprintf("%14s %8s %12s %10s %10s %10s\n",
              "Variable", "Quantile", "Coef.", "Std.Err.", "z-stat", "P>|z|"))
  cat(rep("-", 78), "\n", sep = "")
  
  current_tau <- NA
  for (i in seq_len(nrow(x$lr_table))) {
    row <- x$lr_table[i, ]
    
    if (is.na(current_tau) || row$tau != current_tau) {
      cat(sprintf("  -- tau = %.2f %s\n", row$tau, paste(rep("-", 55), collapse = "")))
      current_tau <- row$tau
    }
    
    est_str <- if (is.na(row$estimate)) "n/a" else sprintf("%10.4f", row$estimate)
    se_str <- if (is.na(row$std_error)) "" else sprintf("%10.4f", row$std_error)
    z_str <- if (is.na(row$z_value)) "" else sprintf("%10.3f", row$z_value)
    p_str <- if (is.na(row$p_value)) "" else sprintf("%10.4f", row$p_value)
    
    cat(sprintf("%14s %8.2f %12s %10s %10s %10s %s\n",
                row$variable, row$tau, est_str, se_str, z_str, p_str, row$signif))
  }
  cat(rep("-", 78), "\n", sep = "")
  cat("Signif. codes: '***' 0.01, '**' 0.05, '*' 0.10\n")
  
  # Table 2: ECT Speed of Adjustment
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  Table 2: ECM Speed of Adjustment rho(tau)\n")
  cat("  rho(tau) should be negative for convergence to equilibrium\n")
  cat(rep("-", 78), "\n", sep = "")
  
  cat(sprintf("%10s %12s %10s %10s %10s %10s %12s\n",
              "Quantile", "rho(tau)", "Std.Err.", "z-stat", "P>|z|", "Half-Life", "Status"))
  cat(rep("-", 78), "\n", sep = "")
  
  for (i in seq_len(nrow(x$ect_table))) {
    row <- x$ect_table[i, ]
    
    rho_str <- if (is.na(row$rho)) "n/a" else sprintf("%10.4f", row$rho)
    se_str <- if (is.na(row$std_error)) "" else sprintf("%10.4f", row$std_error)
    z_str <- if (is.na(row$z_value)) "" else sprintf("%10.3f", row$z_value)
    p_str <- if (is.na(row$p_value)) "" else sprintf("%10.4f", row$p_value)
    hl_str <- if (is.na(row$halflife) || row$halflife <= 0) "Inf" else sprintf("%10.2f", row$halflife)
    
    cat(sprintf("%10.2f %12s %10s %10s %10s %10s %12s\n",
                row$tau, rho_str, se_str, z_str, p_str, hl_str, row$status))
  }
  cat(rep("-", 78), "\n", sep = "")
  cat("Half-life = ln(2)/|rho(tau)| - periods to close 50% of disequilibrium\n")
  
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  XTPQARDL v1.0.1                                    ", x$ardl_order, "\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("\n")
  
  invisible(x)
}


#' Coefficients Method for xtpqardl Objects
#'
#' Extract estimated coefficients from a Panel Quantile ARDL model.
#'
#' @param object An object of class \code{"xtpqardl"}.
#' @param type Character string specifying which coefficients to extract:
#'   \code{"beta"} for long-run coefficients (default), \code{"rho"} for
#'   ECT speed of adjustment, or \code{"all"} for both.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector or list of coefficients.
#'
#' @examples
#' \donttest{
#' data(pqardl_sample)
#' fit <- xtpqardl(
#'   formula = d_y ~ d_x1 + d_x2,
#'   data = pqardl_sample,
#'   id = "country",
#'   time = "year",
#'   lr = c("L_y", "x1", "x2"),
#'   tau = c(0.25, 0.50, 0.75)
#' )
#' coef(fit)
#' coef(fit, type = "rho")
#' }
#'
#' @export
coef.xtpqardl <- function(object, type = c("beta", "rho", "all"), ...) {
  type <- match.arg(type)
  
  # Build coefficient names
  beta_names <- paste0(rep(object$lr_x, length(object$tau)), "_tau",
                        rep(object$tau, each = object$k_lr))
  rho_names <- paste0("rho_tau", object$tau)
  
  beta <- as.vector(object$beta_mg)
  names(beta) <- beta_names
  
  rho <- as.vector(object$rho_mg)
  names(rho) <- rho_names
  
  if (type == "beta") {
    return(beta)
  } else if (type == "rho") {
    return(rho)
  } else {
    return(list(beta = beta, rho = rho))
  }
}


#' Variance-Covariance Matrix for xtpqardl Objects
#'
#' Extract the variance-covariance matrix of the estimated parameters.
#'
#' @param object An object of class \code{"xtpqardl"}.
#' @param type Character string specifying which covariance matrix to extract:
#'   \code{"beta"} for long-run coefficients (default), or \code{"rho"} for
#'   ECT coefficients.
#' @param ... Additional arguments (currently unused).
#'
#' @return A variance-covariance matrix.
#'
#' @export
vcov.xtpqardl <- function(object, type = c("beta", "rho"), ...) {
  type <- match.arg(type)
  
  if (type == "beta") {
    V <- object$beta_V
    if (!is.null(V)) {
      beta_names <- paste0(rep(object$lr_x, length(object$tau)), "_tau",
                            rep(object$tau, each = object$k_lr))
      rownames(V) <- colnames(V) <- beta_names
    }
    return(V)
  } else {
    V <- object$rho_V
    if (!is.null(V)) {
      rho_names <- paste0("rho_tau", object$tau)
      rownames(V) <- colnames(V) <- rho_names
    }
    return(V)
  }
}
