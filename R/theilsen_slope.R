#' @title Theil-Sen slope of rolling change trend
#'
#' @description Calculates the Theil-Sen slope and 95% confidence intervals for the rolling change trend using the
#' \code{openair::TheilSen} function.
#'
#' @param change.trend A rolling change trend plot output from the \code{rolling_change_trend} function.
#'
#' @param avg.time Resolution for which to calculate Theil Sen slope (e.g. yearly, monthly, daily average data).
#'
#' @return Theil-Sen estimator plot and data (slope, 95% confidence intervals, significance level).
#'
#' @export



theilsen_slope <- function(change.trend, avg.time = "year"){

  data <- change.trend$data #%>%
    #dplyr::mutate(date = if(avg.time == "year"){
    #  lubridate::year(date)
    #} else if(avg.time == "month"){
    #  format(lubridate::date(date), "%Y-%m")
    #} else if(avg.time == "day"){
    #  lubridate::date(date)
    #})

  #theilsen <- mblm::mblm(trend_difference~date, data)

  #theilsen.slope <- theilsen$coefficients[2]
  #signif <- Kendall::MannKendall(change_trend$data$trend_difference)$sl

  theilsen <- openair::TheilSen(data, pollutant = "trend_difference", avg.time = avg.time, plot = FALSE)

  theilsen.data <- theilsen$data$res2 %>%
    dplyr::select(slope, lower, upper, slope.percent, lower.percent, upper.percent, p.stars) %>%
    dplyr::filter(!(is.nan(slope)))

  #out <- data.frame("slope" = theilsen.slope, "p-value" = signif[1])
  out <- list("plot" = theilsen$plot, "data" = theilsen.data)

  return(out)

}



