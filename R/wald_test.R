#' Wald Test for Parameter Equality Across Quantiles
#'
#' Performs Wald tests for the null hypothesis that parameters are equal
#' across different quantiles. Tests both long-run coefficients (beta) and
#' ECT speed of adjustment (rho).
#'
#' @param object An object of class \code{"xtpqardl"}.
#' @param joint Logical. If \code{TRUE} (default), performs joint test for
#'   all coefficients. If \code{FALSE}, performs individual tests for each
#'   variable.
#'
#' @return An object of class \code{"wald_test.xtpqardl"} containing:
#' \describe{
#'   \item{beta_test}{Wald test result for long-run coefficients}
#'   \item{rho_test}{Wald test result for ECT coefficients}
#'   \item{individual_beta}{Individual tests for each long-run variable (if joint = FALSE)}
#'   \item{tau}{Quantiles tested}
#' }
#'
#' @details
#' The Wald test statistic is computed as:
#' \deqn{W = (R\theta)'[R V R']^{-1}(R\theta)}
#'
#' where \eqn{\theta} is the vector of coefficients, \eqn{V} is the
#' variance-covariance matrix, and \eqn{R} is a restriction matrix testing
#' equality of coefficients across quantiles.
#'
#' Under the null hypothesis of equal coefficients, \eqn{W} follows a
#' chi-squared distribution with degrees of freedom equal to the number of
#' restrictions.
#'
#' @references
#' Koenker R, Bassett G (1982). "Tests of Linear Hypotheses and L1 Estimation."
#' \emph{Econometrica}, 50(6), 1577-1583. \doi{10.2307/1913398}
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
#' wald_test(fit)
#' }
#'
#' @export
wald_test <- function(object, joint = TRUE) {
  UseMethod("wald_test")
}


#' @rdname wald_test
#' @export
wald_test.xtpqardl <- function(object, joint = TRUE) {
  
  ntau <- length(object$tau)
  k_lr <- object$k_lr
  
  if (ntau < 2) {
    stop("At least 2 quantiles required for Wald test")
  }
  
  # Number of restrictions for equality across quantiles
  n_restrict <- ntau - 1
  
  # Test for beta equality
  beta_result <- .wald_test_equality(object$beta_mg, object$beta_V,
                                       object$tau, k_lr, joint)
  
  # Test for rho equality
  rho_result <- .wald_test_equality(object$rho_mg, object$rho_V,
                                      object$tau, 1, joint = TRUE)
  
  result <- list(
    beta_test = beta_result,
    rho_test = rho_result,
    tau = object$tau,
    k_lr = k_lr,
    lr_x = object$lr_x,
    joint = joint
  )
  
  class(result) <- "wald_test.xtpqardl"
  return(result)
}


#' @keywords internal
.wald_test_equality <- function(theta_mg, theta_V, tau, k, joint) {
  ntau <- length(tau)
  n_restrict <- ntau - 1
  
  if (is.null(theta_V) || any(is.na(theta_V))) {
    return(list(
      statistic = NA,
      df = NA,
      p_value = NA,
      message = "Variance matrix not available"
    ))
  }
  
  theta <- as.vector(theta_mg)
  n_coef <- length(theta)
  
  # For joint test: test all coefficients equal across quantiles
  # For individual: test each variable separately
  
  if (joint) {
    # Build restriction matrix R such that R*theta = 0 under H0
    # For each variable j and quantiles t1, t2, ..., tK:
    # theta_j(t1) = theta_j(t2) = ... = theta_j(tK)
    # This means: theta_j(t2) - theta_j(t1) = 0, theta_j(t3) - theta_j(t1) = 0, ...
    
    n_rows <- k * n_restrict
    R <- matrix(0, nrow = n_rows, ncol = n_coef)
    
    row_idx <- 0
    for (j in seq_len(k)) {
      for (ti in 2:ntau) {
        row_idx <- row_idx + 1
        # Coefficient at quantile ti
        col_ti <- (ti - 1) * k + j
        # Coefficient at quantile 1 (reference)
        col_t1 <- (1 - 1) * k + j
        
        if (col_ti <= n_coef && col_t1 <= n_coef) {
          R[row_idx, col_ti] <- 1
          R[row_idx, col_t1] <- -1
        }
      }
    }
    
    # Compute Wald statistic
    r_theta <- R %*% theta
    r_V_r <- R %*% theta_V %*% t(R)
    
    # Check for singularity
    if (any(is.na(r_V_r)) || det(r_V_r) < 1e-15) {
      return(list(
        statistic = NA,
        df = n_rows,
        p_value = NA,
        message = "Singular variance matrix"
      ))
    }
    
    tryCatch({
      r_V_r_inv <- solve(r_V_r)
      W <- as.numeric(t(r_theta) %*% r_V_r_inv %*% r_theta)
      p_value <- 1 - pchisq(W, df = n_rows)
      
      return(list(
        statistic = W,
        df = n_rows,
        p_value = p_value,
        message = NULL
      ))
    }, error = function(e) {
      return(list(
        statistic = NA,
        df = n_rows,
        p_value = NA,
        message = paste("Error:", e$message)
      ))
    })
    
  } else {
    # Individual tests for each variable
    results <- vector("list", k)
    
    for (j in seq_len(k)) {
      R <- matrix(0, nrow = n_restrict, ncol = n_coef)
      
      for (ti in 2:ntau) {
        row_idx <- ti - 1
        col_ti <- (ti - 1) * k + j
        col_t1 <- (1 - 1) * k + j
        
        if (col_ti <= n_coef && col_t1 <= n_coef) {
          R[row_idx, col_ti] <- 1
          R[row_idx, col_t1] <- -1
        }
      }
      
      r_theta <- R %*% theta
      r_V_r <- R %*% theta_V %*% t(R)
      
      if (any(is.na(r_V_r)) || det(r_V_r) < 1e-15) {
        results[[j]] <- list(
          statistic = NA,
          df = n_restrict,
          p_value = NA
        )
      } else {
        tryCatch({
          r_V_r_inv <- solve(r_V_r)
          W <- as.numeric(t(r_theta) %*% r_V_r_inv %*% r_theta)
          p_value <- 1 - pchisq(W, df = n_restrict)
          results[[j]] <- list(
            statistic = W,
            df = n_restrict,
            p_value = p_value
          )
        }, error = function(e) {
          results[[j]] <- list(
            statistic = NA,
            df = n_restrict,
            p_value = NA
          )
        })
      }
    }
    
    return(results)
  }
}


#' Print Method for wald_test.xtpqardl Objects
#'
#' @param x An object of class \code{"wald_test.xtpqardl"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.wald_test.xtpqardl <- function(x, ...) {
  
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  Wald Tests for Parameter Equality Across Quantiles\n")
  cat(rep("=", 78), "\n", sep = "")
  
  cat("\nQuantiles tested:", paste(sprintf("%.2f", x$tau), collapse = ", "), "\n")
  cat("H0: Parameters are equal across quantiles\n")
  cat("H1: At least one parameter differs\n\n")
  
  cat(rep("-", 78), "\n", sep = "")
  cat("  Test for Long-Run Coefficients beta(tau)\n")
  cat(rep("-", 78), "\n", sep = "")
  
  if (x$joint) {
    bt <- x$beta_test
    if (!is.na(bt$statistic)) {
      cat(sprintf("  Wald statistic:  %.4f\n", bt$statistic))
      cat(sprintf("  Degrees of freedom: %d\n", bt$df))
      cat(sprintf("  P-value:         %.4f", bt$p_value))
      if (bt$p_value < 0.01) {
        cat(" ***\n")
      } else if (bt$p_value < 0.05) {
        cat(" **\n")
      } else if (bt$p_value < 0.10) {
        cat(" *\n")
      } else {
        cat("\n")
      }
      
      if (bt$p_value < 0.05) {
        cat("  => Reject H0: Long-run effects vary significantly across quantiles\n")
      } else {
        cat("  => Fail to reject H0: No significant difference across quantiles\n")
      }
    } else {
      cat("  Test not available:", bt$message, "\n")
    }
  } else {
    # Individual tests
    cat(sprintf("%-15s %12s %8s %12s\n", "Variable", "Wald stat", "df", "P-value"))
    cat(rep("-", 50), "\n", sep = "")
    
    for (j in seq_along(x$beta_test)) {
      bt <- x$beta_test[[j]]
      vname <- if (j <= length(x$lr_x)) x$lr_x[j] else paste0("var", j)
      
      if (!is.na(bt$statistic)) {
        stars <- ""
        if (bt$p_value < 0.01) stars <- "***"
        else if (bt$p_value < 0.05) stars <- "**"
        else if (bt$p_value < 0.10) stars <- "*"
        
        cat(sprintf("%-15s %12.4f %8d %12.4f %s\n",
                    vname, bt$statistic, bt$df, bt$p_value, stars))
      } else {
        cat(sprintf("%-15s %12s %8s %12s\n", vname, "n/a", "", ""))
      }
    }
  }
  
  cat("\n")
  cat(rep("-", 78), "\n", sep = "")
  cat("  Test for ECT Speed of Adjustment rho(tau)\n")
  cat(rep("-", 78), "\n", sep = "")
  
  rt <- x$rho_test
  if (!is.na(rt$statistic)) {
    cat(sprintf("  Wald statistic:  %.4f\n", rt$statistic))
    cat(sprintf("  Degrees of freedom: %d\n", rt$df))
    cat(sprintf("  P-value:         %.4f", rt$p_value))
    if (rt$p_value < 0.01) {
      cat(" ***\n")
    } else if (rt$p_value < 0.05) {
      cat(" **\n")
    } else if (rt$p_value < 0.10) {
      cat(" *\n")
    } else {
      cat("\n")
    }
    
    if (rt$p_value < 0.05) {
      cat("  => Reject H0: Adjustment speed varies significantly across quantiles\n")
    } else {
      cat("  => Fail to reject H0: No significant difference in adjustment speed\n")
    }
  } else {
    cat("  Test not available:", rt$message, "\n")
  }
  
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("Signif. codes: '***' 0.01, '**' 0.05, '*' 0.10\n")
  cat(rep("=", 78), "\n\n", sep = "")
  
  invisible(x)
}


#' Compute Impulse Response Function
#'
#' Computes the impulse response function (IRF) for Panel Quantile ARDL models,
#' showing the response to a one-unit shock via the error correction mechanism.
#'
#' @param object An object of class \code{"xtpqardl"}.
#' @param horizon Integer specifying the number of periods for the IRF.
#'   Default is 20.
#'
#' @return A matrix with rows representing time periods and columns
#'   representing quantiles. Each entry shows the response at that period
#'   for that quantile.
#'
#' @details
#' The IRF for the error correction model is computed as:
#' \deqn{IRF_t(\tau) = (1 + \rho(\tau))^t}
#'
#' which shows the decay of a unit shock over time through the error
#' correction mechanism. Values approach zero as \eqn{t \to \infty} when
#' \eqn{-1 < \rho(\tau) < 0}.
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
#' irf <- compute_irf(fit, horizon = 15)
#' print(irf)
#' }
#'
#' @export
compute_irf <- function(object, horizon = 20) {
  
  if (!inherits(object, "xtpqardl")) {
    stop("'object' must be of class 'xtpqardl'")
  }
  
  ntau <- length(object$tau)
  rho <- as.vector(object$rho_mg)
  
  irf_matrix <- matrix(NA, nrow = horizon + 1, ncol = ntau)
  colnames(irf_matrix) <- paste0("tau_", object$tau)
  rownames(irf_matrix) <- paste0("t_", 0:horizon)
  
  for (ti in seq_along(object$tau)) {
    rho_val <- rho[ti]
    if (!is.na(rho_val) && rho_val < 0 && rho_val > -2) {
      for (t in 0:horizon) {
        irf_matrix[t + 1, ti] <- (1 + rho_val)^t
      }
    }
  }
  
  class(irf_matrix) <- c("irf.xtpqardl", "matrix")
  attr(irf_matrix, "tau") <- object$tau
  attr(irf_matrix, "rho") <- rho
  
  return(irf_matrix)
}


#' Print Method for irf.xtpqardl Objects
#'
#' @param x An object of class \code{"irf.xtpqardl"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.irf.xtpqardl <- function(x, ...) {
  cat("\nImpulse Response Function by Quantile\n")
  cat("Response to a 1-unit shock via ECM mechanism\n")
  cat(rep("-", 50), "\n", sep = "")
  
  tau <- attr(x, "tau")
  
  cat(sprintf("%8s", "Period"))
  for (ti in seq_along(tau)) {
    cat(sprintf(" %12s", paste0("tau=", tau[ti])))
  }
  cat("\n")
  cat(rep("-", 50), "\n", sep = "")
  
  for (t in seq_len(nrow(x))) {
    cat(sprintf("%8d", t - 1))
    for (ti in seq_len(ncol(x))) {
      val <- x[t, ti]
      if (is.na(val)) {
        cat(sprintf(" %12s", "div."))
      } else {
        cat(sprintf(" %12.4f", val))
      }
    }
    cat("\n")
  }
  
  cat(rep("-", 50), "\n", sep = "")
  
  invisible(x)
}
