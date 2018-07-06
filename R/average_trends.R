#' @title Function to compare average trend with long term sites' trend (V2: 04-07-18)
#'
#' @description Produces plots of annual average trends in the pollutant of interest aggregated over all sites, and
#' over the long term sites only (i.e. sites open constantly over the period of interest).
#' Changes to v2:
#' \enumerate{
#'     \item Better organisation: generalised averaging and plotting functions (in helper functions) - both for averaging for each
#'     site and over all sites.
#' }
#'
#' @param obs  A data frame of ambient pollutant concentration data. Must contain the columns: site_code, date,
#' value. If 'pollutant' is a pollutant ratio, the data frames of the corresponding pollutants
#' must be supplied as a list of data frames in the order they are given in the ratio. E.g. for
#' \code{pollutant = "no2/nox"}, \code{obs = list(obs.no2, obs.nox)}.
#'
#' @param pollutant The pollutant of interest (character string). To calculate rolling change trend for a pollutant ratio,
#' separate the two pollutants with a forward slash e.g. \code{pollutant = "no2/nox"}.
#'
#' @param stat The metric used to average the ambient concentration data by year. Options: "median", "mean".
#'
#' @param start.date The starting date of the period of interest. Data capture to filter long term sites will be
#' conducted using this date as the starting date.
#'
#' @param end.date The end date of the period of interest. Data capture to filter long term sites will be
#' conducted using this date as the ending date.
#'
#' @param data.capture The data capture threshold used to define 'long term sites' (as a %). To qualify as a
#' long term site, the monitoring site must have data capture at least equal to this threshold over the period
#' defined by the 'start.date' and 'end.date' arguments.
#'
#' @param smooth.method The smoothing method to use in the \code{geom_smooth} function when plotting trends.
#'
#' @examples
#' \dontrun{
#' average_trends(obs.nox, pollutant = "nox",
#' start.date = "2000-01-01", end.date = "2017-12-31",
#' data.capture = 90)
#' }
#'
#' @import dplyr ggplot2
#'
#' @return A list of plots:
#' \itemize{
#'     \item \code{average.trend} - The average trend over all sites measuring during the period of interest.
#'     \item \code{longterm.trend} - The average trend over all long term sites measuring constantly throughout the period.
#'     \item \code{trends.comparison} - Comparison of the \code{average.trend} and the \code{longterm.trend}.
#'     \item \code{average.allsites} - The average trend for each individual sites measuring during the period of interest.
#'     \item \code{longterm.allsites} - The average trend for each individual long term site measuring constantly throughout the period.
#'  }
#'  Note: if there are no long term sites available for the period of interest, only the \code{average.trend} and
#'  \code{average.allsites} will be returned.
#'
#' @export

average_trends <- function(obs,
                         pollutant,
                         stat = "median",
                         start.date = "2000-01-01",
                         end.date = "2017-12-31",
                         data.capture = 90,
                         smooth.method = "loess"){

  ## Check arguments
  check_arguments(obs = obs,
                  pollutant = pollutant,
                  stat = stat,
                  start.date = start.date,
                  end.date = end.date,
                  data.capture = data.capture,
                  smooth.method = smooth.method)


  ## Data averaging (conditional on whether to calculate for pollutant or pollutant ratio)

  if(!(stringr::str_detect(pollutant, "/"))){

    allsites.sites <- average_data(obs, avg.time = "year", statistic = stat, sites = TRUE) # annual average for each site
    allsites.av <- average_data(obs, avg.time = "year", statistic = stat) # annual average over all sites

    longterm <- data_capture(obs, threshold = data.capture, start.date = start.date, end.date = end.date)

    longterm.sites <- longterm %>% average_data(avg.time = "year", statistic = stat, sites = TRUE) # annual average for each long term sites
    longterm.av <- longterm %>% average_data(longterm.sites, avg.time = "year", statistic = stat) # annual average over all long term sites

  } else{

    allsites.sites <- calculate_pollutant_ratio(obs, avg.time = "year", statistic = stat, sites = TRUE)
    allsites.av <- calculate_pollutant_ratio(obs, avg.time = "year", statistic = stat)

    longterm.sites <- pollutant_ratio_data_capture(obs, avg.time = "year", start.date, end.date, data.capture,
                                                   statistic = stat, sites = TRUE)
    longterm.av <- pollutant_ratio_data_capture(obs, avg.time = "year", start.date, end.date, data.capture, statistic = stat)

  }



  ## Plots

  allsites.av.plot <- trends_plots_helper(allsites.av, sites = FALSE, pollutant, stat, smooth.method, start.date, end.date)
  allsites.sites.plot <- trends_plots_helper(allsites.sites, sites = TRUE, pollutant, stat, smooth.method, start.date, end.date)

  # If long term sites...
  if(nrow(longterm.av) > 0){

    longterm.av.plot <- trends_plots_helper(longterm.av, sites = FALSE, pollutant, stat, smooth.method, start.date, end.date)
    longterm.sites.plot <- trends_plots_helper(longterm.sites, sites = TRUE, pollutant, stat, smooth.method, start.date, end.date)

    # Compare NOx time series from all monitoring sites and long term sites
    trends <- cowplot::plot_grid(allsites.av.plot,
                                 longterm.av.plot,
                                 labels = c("(a) all sites", "(b) long term sites"))

    out <- list("average.trend" = allsites.av.plot,
                "longterm.trend" = longterm.av.plot,
                "trends.comparison" = trends,
                "average.allsites" = allsites.sites.plot,
                "longterm.allsites" = longterm.sites.plot)

  } else{ # if no long term sites....
    print(paste0("No long term sites over the period ", start.date, " - ", end.date, ". Returning average trend only."))

    out <- list("average.trend" = allsites.av.plot,
                "average.allsites" = allsites.sites.plot)
  }


  return(out)
}

