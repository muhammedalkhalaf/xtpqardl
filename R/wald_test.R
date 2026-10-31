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
#' @param all Logical, for \code{model = "dfe"} only. If \code{FALSE}
#'   (default) the Koenker-Bassett tests in \code{kb_test} restrict only the
#'   economic slopes of the pooled regression (speed of adjustment, long-run
#'   levels, short run). If \code{TRUE} they restrict every non-intercept
#'   coefficient, including the \eqn{N - 1} panel-intercept contrasts, which
#'   is what \code{quantreg::anova(..., joint = TRUE)} does on the equivalent
#'   regression \code{y ~ slopes + factor(id)}.
#'
#' @return An object of class \code{"wald_test.xtpqardl"} containing:
#' \describe{
#'   \item{beta_test}{Wald test result for long-run coefficients}
#'   \item{rho_test}{Wald test result for ECT coefficients}
#'   \item{individual_beta}{Individual tests for each long-run variable (if joint = FALSE)}
#'   \item{kb_test}{For \code{model = "dfe"}: Koenker and Bassett (1982)
#'     equality tests across quantiles on the pooled quantile-regression
#'     coefficients, computed from the joint covariance \code{dfe_coef_V}
#'     (estimator \code{dfe_cov} of the fit): \code{joint} (one row: F
#'     statistic, numerator and denominator degrees of freedom, p-value) and
#'     \code{individual} (one row per restricted coefficient), the names of
#'     the restricted coefficients in \code{restricted}, the estimator in
#'     \code{se} and whether the panel contrasts were included in
#'     \code{all}. \code{NULL} for MG and PMG.}
#'   \item{tau}{Quantiles tested}
#' }
#'
#' @details
#' The Wald test statistic is computed as:
#' \deqn{W = (R\theta)'[R V R']^{-1}(R\theta)}
#'
#' where \eqn{\theta} is the vector of coefficients, \eqn{V} is the
#' variance-covariance matrix, and \eqn{R} is a restriction matrix testing
#' equality of coefficients across quantiles (each quantile against the
#' first). Under the null hypothesis of equal coefficients, \eqn{W} follows a
#' chi-squared distribution with degrees of freedom equal to the number of
#' restrictions. For MG the covariance is the Pesaran-Smith between-panel
#' covariance of the stacked panel estimates; for DFE it is the joint
#' covariance described in \code{\link{xtpqardl}} (one estimator for every
#' block, named in \code{dfe_cov}); for PMG the cross-quantile blocks are not
#' estimated and \code{beta_test} is not available.
#'
#' \strong{Koenker-Bassett tests (DFE).} \code{kb_test} applies the same
#' restriction scheme to the pooled quantile-regression coefficients
#' themselves (rho, theta, short run) with the joint covariance
#' \eqn{W = [(\min(\tau_a,\tau_b) - \tau_a\tau_b) H_a^{-1} X'X H_b^{-1}]_{ab}}
#' of \code{\link{xtpqardl}} (Koenker and Bassett 1982, Theorem 3.1, with
#' the estimator chosen by \code{dfe_cov}). As in \code{quantreg::anova} the
#' statistic is reported as \eqn{T_n/ndf} against an F distribution with
#' \eqn{ndf} restrictions and \eqn{ddf = nm - ndf} denominator degrees of
#' freedom (\eqn{n} observations, \eqn{m} quantiles); with
#' \code{dfe_cov = "nid"} and \code{all = TRUE} the joint and individual
#' rows reproduce \code{anova(rq_1, ..., rq_m, se = "nid")} on the
#' equivalent regression with \code{factor(id)}. With the default
#' \code{all = FALSE} the joint test has \eqn{ndf = k (m - 1)} where \eqn{k}
#' is the number of economic slopes, so the hypothesis is the one stated in
#' the printout, and the \eqn{N - 1} contrasts (whose equality across
#' quantiles is of no economic interest and dilutes the test) are not
#' restricted.
#'
#' \strong{Rank of the theta covariance.} Cho, Kim and Shin (2015, Theorem
#' 1(ii)) show that the asymptotic covariance of the long-run level
#' coefficients \eqn{\gamma(\tau)} (theta) has rank one under their
#' normalisation. In their words, "has a rank equal to unity by Theorem 1(ii).
#' This implies that, if R has a rank greater than unity, the corresponding
#' Wald test has a degenerate asymptotic null distribution." (the page numbers
#' of that passage are not given here; the library copy has no page markers).
#' The individual rows of \code{kb_test} for each theta coefficient restrict
#' one coefficient across quantiles; the joint test with more than one
#' long-run regressor restricts several theta coefficients at once and relies
#' on the finite-sample full rank of their sandwich covariance.
#'
#' The audit's Monte Carlo is consistent with this concern, but does not
#' isolate it. In the \code{kb_joint_all} column of its size grid (nid
#' estimator, all non-intercept coefficients restricted, iid N(0,1) errors,
#' N = 5, T = 20, R = 500) the rejection rate at nominal 5 percent was 16.2
#' percent with one long-run regressor and 23.2 percent with two. The two
#' designs also differ in the number of restrictions and the nid density
#' bias documented in the audit is a confounder, so the increase should be
#' read as consistent with, not as evidence for, the rank argument. The
#' long-run beta = -theta/rho is tested in \code{beta_test}, with the delta
#' method.
#'
#' @references
#' Cho JS, Kim TH, Shin Y (2015). "Quantile Cointegration in the
#' Autoregressive Distributed-Lag Modeling Framework." \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.05.003}
#'
#' Koenker R, Bassett G (1982). "Robust Tests for Heteroscedasticity Based on
#' Regression Quantiles." \emph{Econometrica}, 50(1), 43-61.
#' \doi{10.2307/1912528}
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
wald_test <- function(object, joint = TRUE, all = FALSE) {
  UseMethod("wald_test")
}


#' @rdname wald_test
#' @export
wald_test.xtpqardl <- function(object, joint = TRUE, all = FALSE) {

  ntau <- length(object$tau)
  k_lr <- object$k_lr

  if (ntau < 2) {
    stop("At least 2 quantiles required for Wald test")
  }

  # Test for beta equality
  beta_result <- .wald_test_equality(object$beta_mg, object$beta_V,
                                       object$tau, k_lr, joint)

  # Test for rho equality
  rho_result <- .wald_test_equality(object$rho_mg, object$rho_V,
                                      object$tau, 1, joint = TRUE)

  ## Koenker-Bassett (1982) tests on the underlying quantile regression
  ## coefficients, available when the fits share one design (model = "dfe").
  kb <- if (!is.null(object$dfe_coef_V)) .kb_tests(object, all = isTRUE(all)) else NULL

  result <- list(
    beta_test = beta_result,
    rho_test = rho_result,
    kb_test = kb,
    tau = object$tau,
    k_lr = k_lr,
    lr_x = object$lr_x,
    lr_y = object$lr_y,
    model = object$model,
    joint = joint
  )

  class(result) <- "wald_test.xtpqardl"
  return(result)
}


#' Koenker-Bassett equality tests across quantiles from the joint covariance
#'
#' Restrictions: consecutive differences across quantiles, as in
#' \code{quantreg::anova.rqlist}; statistic Tn/ndf against F(ndf, n m - ndf).
#' @keywords internal
#' @noRd
.kb_tests <- function(object, all = FALSE) {
  B <- object$dfe_coef
  V <- object$dfe_coef_V
  tau <- object$tau
  m <- length(tau); pn <- nrow(B); n <- object$n_obs
  cn <- rownames(B)
  sel <- if (all) setdiff(cn, "(Intercept)") else object$dfe_slopes
  b <- as.numeric(B)
  Dm <- diff(diag(m))
  one_test <- function(cols) {
    E <- matrix(0, length(cols), pn)
    E[cbind(seq_along(cols), match(cols, cn))] <- 1
    D <- kronecker(Dm, E)
    ndf <- nrow(D)
    ddf <- n * m - ndf
    DW <- D %*% V %*% t(D)
    if (anyNA(DW) || rcond(DW) < 1e-12) {
      return(data.frame(ndf = ndf, ddf = ddf, Tn = NA_real_, pvalue = NA_real_))
    }
    Db <- D %*% b
    Tn <- as.numeric(t(Db) %*% solve(DW, Db)) / ndf
    data.frame(ndf = ndf, ddf = ddf, Tn = Tn, pvalue = stats::pf(Tn, ndf, ddf, lower.tail = FALSE))
  }
  jt <- one_test(sel)
  ind <- do.call(rbind, lapply(sel, function(cc) one_test(cc)))
  rownames(ind) <- sel
  rownames(jt) <- "joint"
  msg <- if (anyNA(V[rep(match(sel, cn), m) + rep((seq_len(m) - 1) * pn, each = length(sel)),
                      rep(match(sel, cn), m) + rep((seq_len(m) - 1) * pn, each = length(sel))])) {
    "joint covariance not available for some quantiles (see panel_cov_status)"
  } else NULL
  list(joint = jt, individual = ind, restricted = sel, all = all,
       se = object$dfe_cov, coef_names = cn, message = msg,
       note = paste0("Tests on the pooled DFE quantile-regression coefficients",
                     if (all) " including the panel-intercept contrasts" else
                       " (rho, theta, short run; panel intercepts not restricted)",
                     ", joint covariance estimator dfe_cov = \"", object$dfe_cov,
                     "\"; the long-run beta = -theta/rho is tested in beta_test."))
}


#' @keywords internal
.wald_test_equality <- function(theta_mg, theta_V, tau, k, joint) {
  ntau <- length(tau)
  n_restrict <- ntau - 1

  if (is.null(theta_V) || any(is.na(theta_V))) {
    na_one <- list(statistic = NA, df = NA, p_value = NA,
                   message = "Variance matrix not available")
    ## joint = FALSE must return one element per variable (1.0.3 returned a
    ## single flat list here, which made print() fail with
    ## "$ operator is invalid for atomic vectors")
    if (!joint) return(rep(list(na_one), k))
    return(na_one)
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

    # Check for singularity (scale-free: det() of a matrix of small
    # variances is tiny even when the matrix is well conditioned)
    if (any(is.na(r_V_r)) || rcond(r_V_r) < 1e-12) {
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

      if (any(is.na(r_V_r)) || rcond(r_V_r) < 1e-12) {
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
        cat("  Reject H0 at 5 percent: long-run effects differ across quantiles\n")
      } else {
        cat("  Do not reject H0 at 5 percent: no evidence that long-run effects differ\n")
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

      if (is.list(bt) && !is.null(bt$statistic) && !is.na(bt$statistic)) {
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
      cat("  Reject H0 at 5 percent: adjustment speed differs across quantiles\n")
    } else {
      cat("  Do not reject H0 at 5 percent: no evidence that adjustment speed differs\n")
    }
  } else {
    cat("  Test not available:", rt$message, "\n")
  }

  if (!is.null(x$kb_test)) {
    kb <- x$kb_test
    cat("\n")
    cat(rep("-", 78), "\n", sep = "")
    cat("  Koenker-Bassett tests on the pooled quantile-regression coefficients\n")
    cat(sprintf("  (DFE design shared across quantiles; joint covariance dfe_cov = \"%s\")\n", kb$se))
    cat(rep("-", 78), "\n", sep = "")
    if (!is.null(kb$joint)) {
      jt <- kb$joint
      lab <- if (isTRUE(kb$all)) {
        sprintf("Joint equality of all %d non-intercept coefficients (slopes and panel contrasts):",
                length(kb$restricted))
      } else {
        sprintf("Joint equality of the %d economic slopes (rho, theta, short run):",
                length(kb$restricted))
      }
      cat("  ", lab, "\n", sep = "")
      if (is.na(jt[1, "Tn"])) {
        cat("    not available:", if (!is.null(kb$message)) kb$message else "singular covariance", "\n")
      } else {
        cat(sprintf("    F = %.4f  df = (%d, %d)  p = %.4f\n",
                    jt[1, "Tn"], as.integer(jt[1, "ndf"]), as.integer(jt[1, "ddf"]), jt[1, "pvalue"]))
      }
    }
    if (!is.null(kb$individual)) {
      tab <- kb$individual
      cat(sprintf("  %-20s %10s %6s %6s %10s\n", "Coefficient", "F", "ndf", "ddf", "p-value"))
      for (i in seq_len(nrow(tab))) {
        if (is.na(tab[i, "Tn"])) {
          cat(sprintf("  %-20s %10s %6d %6d %10s\n", rownames(tab)[i], "n/a",
                      as.integer(tab[i, "ndf"]), as.integer(tab[i, "ddf"]), ""))
        } else {
          cat(sprintf("  %-20s %10.4f %6d %6d %10.4f\n", rownames(tab)[i], tab[i, "Tn"],
                      as.integer(tab[i, "ndf"]), as.integer(tab[i, "ddf"]), tab[i, "pvalue"]))
        }
      }
    }
    if (!isTRUE(kb$all)) {
      cat("  Panel-intercept contrasts are not restricted (wald_test(..., all = TRUE) restricts them).\n")
    }
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
#' correction mechanism. Values approach zero as \eqn{t} grows when
#' \eqn{-2 < \rho(\tau) < 0}; for \eqn{-2 < \rho(\tau) \le -1} the response
#' oscillates in sign (the half-life of \code{\link{xtpqardl}} is defined
#' for \eqn{-1 < \rho < 0} only). Outside \eqn{(-2, 0)} the entries are
#' \code{NA} (printed as "div.").
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
