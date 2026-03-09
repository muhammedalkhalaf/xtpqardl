#' Simulated Panel Data for PQARDL Estimation
#'
#' A simulated panel dataset for demonstrating Panel Quantile ARDL estimation.
#' Contains 10 countries observed over 30 years with variables suitable for
#' error correction modeling.
#'
#' @format A data frame with 300 rows and 9 variables:
#' \describe{
#'   \item{country}{Factor indicating the panel unit (10 countries)}
#'   \item{year}{Integer year variable (1990-2019)}
#'   \item{y}{Dependent variable in levels (e.g., GDP per capita)}
#'   \item{x1}{First explanatory variable in levels (e.g., investment)}
#'   \item{x2}{Second explanatory variable in levels (e.g., trade openness)}
#'   \item{L_y}{Lagged dependent variable (y at t-1)}
#'   \item{d_y}{First difference of y}
#'   \item{d_x1}{First difference of x1}
#'   \item{d_x2}{First difference of x2}
#' }
#'
#' @details
#' The data are simulated from a panel error correction model with
#' heterogeneous adjustment speeds across countries. The true long-run
#' relationship is:
#'
#' \deqn{y_{it} = \beta_1 x_{1,it} + \beta_2 x_{2,it} + \mu_i + \varepsilon_{it}}
#'
#' with error correction dynamics:
#'
#' \deqn{\Delta y_{it} = \rho_i (y_{i,t-1} - \beta_1 x_{1,i,t-1} - \beta_2 x_{2,i,t-1}) 
#' + \gamma_1 \Delta x_{1,it} + \gamma_2 \Delta x_{2,it} + u_{it}}
#'
#' where \eqn{\rho_i \sim U(-0.6, -0.2)} varies by panel.
#'
#' @source Simulated data for package demonstration.
#'
#' @examples
#' data(pqardl_sample)
#' head(pqardl_sample)
#' 
#' # Check panel structure
#' table(pqardl_sample$country)
#'
#' @keywords datasets
"pqardl_sample"
