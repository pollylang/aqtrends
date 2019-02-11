#' @title Compare average trend over all sites with average trend over long term sites
#'
#' @description Produces plots of annual average trends in the pollutant of interest aggregated over all sites, and
#' over the long term sites only (i.e. sites open constantly over the period of interest).
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
#' @param avg.time Resolution to average data to. Options: "year", "month", "day".
#'
#' @examples
#' \dontrun{
#' average_trends(london_nox_data, pollutant = "nox", start.date = "2000-01-01", end.date = "2017-12-31")
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
                         avg.time = "year"){

  ## Check arguments
  check_arguments(obs = obs,
                  pollutant = pollutant,
                  stat = stat,
                  start.date = start.date,
                  end.date = end.date)


  ## Data averaging (conditional on whether to calculate for pollutant or pollutant ratio)

  start <- lubridate::year(lubridate::ymd(start.date))
  end <- lubridate::year(lubridate::ymd(end.date))

  worker <- function(id = "all", include.sites = FALSE){

    if(!(stringr::str_detect(pollutant, "/"))){

      obs <- obs %>%
        dplyr::filter(date >= lubridate::ymd(start.date), date <= lubridate::ymd(end.date))

      if(id == "longterm"){
        longterm.sitecodes <- get_longterm_sites(obs, pollutant, start, end, resolution = avg.time)
        obs <- obs[obs$site_code %in% longterm.sitecodes, ]
      }

      if(nrow(obs) >= 1){            # if there are long term sites...
        obs <- average_data(obs, avg.time = avg.time, statistic = stat, sites = include.sites) # annual average for each site
      } else{                        # if there aren't any long term sites...
        return(NULL)
      }


    } else{

      # Get long term sites
      obs <- purrr::map(obs, function(x) dplyr::filter(x, date >= lubridate::ymd(start.date), date <= lubridate::ymd(end.date)))

      if(id == "longterm"){
        pollutants <- purrr::map_chr(obs, function(x) unique(x$variable))
        longterm.sitecodes <- purrr::map2(obs, pollutants, get_longterm_sites,
                                         start = start, end = end,
                                         resolution = avg.time) %>%
                                         {Reduce(intersect, .)}

        obs <- purrr::map(obs, function(x) x[x$site_code %in% longterm.sitecodes, ])
      }

      if(all(sapply(obs, function(x) nrow(x) > 1))){        # if there are long term sites...
        obs <- calculate_pollutant_ratio(obs, avg.time = avg.time, statistic = stat, sites = include.sites)
      } else {                                              # if there aren't any long term sites...
        return(NULL)
      }

    }

    ## Plots
    plot <- trends_plots_helper(obs, sites = include.sites, pollutant, stat, start.date, end.date)

    return(plot)

  }

  # Create plots for (i) all sites/overall trend, (ii) all sites/individual site trends, (iii) long term sites/overall trend,
  # (iv) long term sites/individual site trends
  args <- data.frame("id" = c(rep("all", 2), rep("longterm", 2)),
                     "include_sites" = rep(c(FALSE, TRUE)))

  plots <- purrr::map2(args$id, args$include_sites, worker) %>% plyr::compact()

  if(length(plots) == 4){
    # Compare NOx time series from all monitoring sites and long term sites
    trends <- cowplot::plot_grid(plots[[1]], plots[[3]],
                                 labels = c("All sites", "Long term sites"))

    out <- list("average.trend" = plots[[1]],
                "longterm.trend" = plots[[3]],
                "trends.comparison" = trends,
                "average.allsites" = plots[[2]],
                "longterm.allsites" = plots[[4]])

  } else if(length(plots) == 2){ # if no long term sites....

    print(paste0("No long term sites over the period ", start.date, " - ", end.date, ". Returning average trend only."))

    out <- list("average.trend" = plots[[1]],
                "longterm.trend" = plots[[2]])

  }
  return(out)
}

