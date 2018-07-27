#' @title Theil-Sen slope of rolling change trend
#'
#' @description Calculates the Theil-Sen slope and p-value for the rolling change trend using the \code{mblm}
#' and \code{Kendall} packages.
#'
#' @details The two-sided p-value is calculated using Mann-Kendall trend test from
#' the \href{https://cran.r-project.org/web/packages/Kendall/Kendall.pdf}{Kendall} package. This is used in preference
#' to the confidence interval calculation provided by the \code{mblm} package, which suffers from inaccuracies
#' due to autocorrelation of the time series data.
#'
#' @param change_trend A rolling change trend plot output from the \code{rolling_change_trend} function.
#'
#' @return Data frame with two columns: the Theil-Sen slope and the two-sided p-value
#'
#' @export



theilsen_slope <- function(change_trend){

  theilsen <- mblm::mblm(trend_difference~date, change_trend$data)

  theilsen.slope <- theilsen$coefficients[2]
  signif <- Kendall::MannKendall(change_trend$data$trend_difference)$sl

  out <- data.frame("slope" = theilsen.slope, "p-value" = signif[1])

  return(out)

}



