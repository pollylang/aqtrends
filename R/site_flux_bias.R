#' @title The difference in concentration between opening and closing sites as a function of year
#'
#' @description Detects a difference in average concentration of the pollutant of interest between sites opening
#' and closing in a given year (i.e. a bias towards more/less polluted locations) and any change in this bias over
#' time.
#'
#' @param obs A data frame of ambient pollutant concentration data. Must contain the columns: site_code, date,
#' value. If 'pollutant' is a pollutant ratio, the data frames of the corresponding pollutants
#' must be supplied as a list of data frames in the order they are given in the ratio. E.g. for
#' \code{pollutant = "no2/nox"}, \code{obs = list(obs.no2, obs.nox)}.
#'
#' @param pollutant The pollutant of interest (character string). To calculate rolling change trend for a pollutant ratio,
#' separate the two pollutants with a forward slash e.g. \code{pollutant = "no2/nox"}.
#'
#' @param stat The metric (character string) used to average the ambient concentration data by year. Options: "median", "mean".
#'
#' @return A list of plots:
#' \itemize{
#'     \item \code{difference} - a plot of the difference between the average concentration of opening sites and closing
#'     sites as a function of year.
#'     \item \code{cumulative_difference} - the cumulative sum of differences in average concentration of opening sites and
#'     closing sites as a function of year.
#'     \item \code{weighted_cumulative_difference} - the cumulative sum of differences in average concentration of opening and
#'     closing sites weighted by the number of opening/closing sites as a function of year.
#'     \item \code{compare} - a side-by-side comparison of the \code{difference}, \code{cumulative_difference} and
#'     \code{weighted_cumulative_difference} plots.
#' }
#'
#' @import dplyr
#'
#'@examples
#'\dontrun{
#'site_flux_bias(london_nox_data, pollutant = "nox", stat = "median")
#'}
#'
#' @export


site_flux_bias <- function(obs,
                           pollutant,
                           stat = "median"){

  ## Check arguments
  check_arguments(obs = obs,
                  pollutant = pollutant,
                  stat = stat)

  ## Create data frame of moving sites as a function of year (i.e. which sites open/close in each year)

  if(!(stringr::str_detect(pollutant, "/"))){
    annual <- average_data(obs, avg.time = "year", statistic = stat, sites = TRUE)
  } else{
    annual <- calculate_pollutant_ratio(obs, avg.time = "year", statistic = stat, sites = TRUE)
  }

  start.end <- annual %>%
    dplyr::group_by(site_code) %>%
    dplyr::summarise(start_year = min(date), end_year = max(date))

  # Join observations data to moving sites data frame
  moving.sites <- start.end %>%
    tidyr::gather("change", "date", 2:3) %>%
    dplyr::mutate(change = ifelse(change == "start_year", "opening", "closing")) %>%
    dplyr::filter(!is.na(date)) %>%
    dplyr::left_join(annual, by = c("site_code" = "site_code", "date" = "date")) %>%
    dplyr::select(c(site_code, date, change, av_value))


  ## Data

  difference.counts <- plyr::count(moving.sites, c("date", "change")) %>%
    tidyr::spread(key = change, value = freq) %>%
    dplyr::rename(count.close = closing, count.open = opening)

  if((difference.counts %>%
      dplyr::filter(!(is.na(count.close)) & !(is.na(count.open))) %>%
    nrow()) == 0){
    stop("The input data frame contains no years in which sites both open AND close, therefore the difference between opening and closing sites cannot be calculated.")
  }

  difference <- moving.sites %>%
    tidyr::spread(change, av_value) %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(opening = ifelse(stat == "median", median(opening, na.rm = TRUE), mean(opening, na.rm = TRUE)),
                     closing = ifelse(stat == "median", median(closing, na.rm = TRUE), mean(closing, na.rm = TRUE))) %>%
    dplyr::mutate(difference = opening - closing) %>%
    dplyr::mutate(difference0 = ifelse(is.na(difference), 0, difference)) %>%
    dplyr::mutate(cumulative_difference = cumsum(difference0)) %>%
    dplyr::mutate(cumulative_difference = ifelse(cumulative_difference == 0, NA, cumulative_difference)) %>%
    dplyr::left_join(difference.counts, by = c("date" = "date")) %>%
    dplyr::mutate(weighted.opening = opening * count.open,
                  weighted.closing = closing * count.close) %>%
    dplyr::mutate(weighted.difference = weighted.opening - weighted.closing) %>%
    dplyr::mutate(weighted.difference0 = ifelse(is.na(weighted.difference), 0, weighted.difference)) %>%
    dplyr::mutate(w.cumulative_difference = cumsum(weighted.difference0)) %>%
    dplyr::mutate(w.cumulative_difference = ifelse(w.cumulative_difference == 0, NA, w.cumulative_difference))

  start.yr <- min(difference$date, na.rm = TRUE)
  end.yr <- max(difference$date, na.rm = TRUE)

  ## Plots

  # Difference between opening/closing sites
  plots <- function(y.var){

    if(y.var == "difference"){
      y.title <- openair::quickText(paste0("Difference in ", pollutant, " concentration (ugm-3)"))
    } else if(y.var == "cumulative_difference"){
      y.title <- openair::quickText(paste0("Cumulative difference in ", pollutant, " concentration (ugm-3)"))
    } else{
      y.title <- openair::quickText(paste0("Cumulative weighted difference in ", pollutant, " concentration (ugm-3)"))
    }

    if(nrow(difference) <= 5){
      smooth.method <- "gam"
    } else{
      smooth.method <- "loess"
    }

    plot.colour <- viridis::inferno(1, begin = 0.3, end = 0.8)

    plot <- ggplot2::ggplot(difference, ggplot2::aes_string(x = "date", y = y.var)) +
      ggplot2::geom_point(color = plot.colour) +
      ggplot2::geom_smooth(method = smooth.method,
                           color =  plot.colour,
                           fill =  plot.colour,
                           alpha=0.25) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      ggplot2::labs(y = y.title, x = "Year") +
      ggplot2::theme_minimal()

    return(plot)

  }

  difference.plot <- plots("difference")
  cum.difference <- plots("cumulative_difference")
  weighted.cum.diff <- plots("w.cumulative_difference")

  compare <- cowplot::plot_grid(difference.plot, cum.difference, weighted.cum.diff, ncol = 3)

  out <- list("difference" = difference.plot,
              "cumulative_difference" = cum.difference,
              "weighted_cumulative_difference" = weighted.cum.diff,
              "all" = compare)

  return(out)

}

