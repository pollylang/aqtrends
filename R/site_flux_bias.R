#' @title The difference between opening and closing sites as a function of year (v2: 05-07-2018)
#'
#' @description Detects a difference in average concentration of the pollutant of interest between sites opening
#' and closing in a given year (i.e. a bias towards more/less polluted locations) and any change in this bias over
#' time. Changes to v2:
#' \enumerate{
#'     \item Improved organisation: averaging and pollutant ratio functions wrapped in helper functions (generalised across
#'     all functions in the package to reduce redundancy).
#'     \item Plot layouts and aesthetics modified to match those of the other functions in the package.
#'     \item Generalised plotting function and applied over different parameters to produce the three plots (more efficient)
#'     \item Removed extraneous code and redundancies
#'     \item Included the option to specify the smoothing method when plotting - useful for data frames with few (moving) sites
#'     to be able to use GAM rather than loess.
#'     \item Slightly altered the output for greater consistency and flexibility of use.
#'  }
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
#' @param smooth.method The smoothing method to use in the \code{geom_smooth} function when plotting trends.
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
#' @import dplyr ggplot2
#'
#'@examples
#'\dontrun{
#'site_flux_bias(obs.nox, pollutant = "nox", stat = "median", smooth.method = "gam")
#'}
#'
#' @export


site_flux_bias <- function(obs,
                              pollutant,
                              stat = "median",
                              smooth.method = "loess"){

  ## Check arguments
  check_arguments(obs = obs,
                  pollutant = pollutant,
                  stat = stat,
                  smooth.method = smooth.method)


  ## Create data frame of moving sites as a function of year (i.e. which sites open/close in each year)

  if(!(stringr::str_detect(pollutant, "/"))){
    annual <- average_data(obs, avg.time = "year", statistic = stat, sites = TRUE)
  } else{
    annual <- calculate_pollutant_ratio(obs, avg.time = "year", statistic = stat, sites = TRUE)
  }

  start.end <- annual %>%
    group_by(site_code) %>%
    summarise(start_year = min(date), end_year = max(date))

  # Join observations data to moving sites data frame
  moving.sites <- start.end %>%
    tidyr::gather("change", "date", 2:3) %>%
    mutate(change = ifelse(change == "start_year", "opening", "closing")) %>%
    filter(!is.na(date)) %>%
    left_join(annual, by = c("site_code" = "site_code", "date" = "date")) %>%
    dplyr::select(c(site_code, date, change, av_value))


  ## Data

  difference.counts <- plyr::count(moving.sites, c("date", "change")) %>%
    tidyr::spread(key = change, value = freq) %>%
    rename(count.close = closing, count.open = opening)

  if((difference.counts %>%
    filter(!(is.na(count.close)) & !(is.na(count.open))) %>%
    nrow()) == 0){
    stop("The input data frame contains no years in which sites both open AND close, therefore the difference between opening and closing sites cannot be calculated.")
  }

  difference <- moving.sites %>%
    tidyr::spread(change, av_value) %>%
    group_by(date) %>%
    summarise(opening = ifelse(stat == "median", median(opening, na.rm = TRUE), mean(opening, na.rm = TRUE)),
              closing = ifelse(stat == "median", median(closing, na.rm = TRUE), mean(closing, na.rm = TRUE))) %>%
    mutate(difference = opening - closing) %>%
    mutate(difference0 = ifelse(is.na(difference), 0, difference),
           cumulative_difference = cumsum(difference0),
           cumulative_difference = ifelse(cumulative_difference == 0, NA, cumulative_difference)) %>%
    left_join(difference.counts, by = c("date" = "date")) %>%
    mutate(weighted.opening = opening * count.open,
           weighted.closing = closing * count.close) %>%
    mutate(weighted.difference = weighted.opening - weighted.closing) %>%
    mutate(weighted.difference0 = ifelse(is.na(weighted.difference), 0, weighted.difference),
           w.cumulative_difference = cumsum(weighted.difference0),
           w.cumulative_difference = ifelse(w.cumulative_difference == 0, NA, w.cumulative_difference))

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

    plot <- ggplot(difference, aes_string(x = "date", y = y.var)) +
      geom_point() +
      geom_smooth(method = smooth.method) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = round(seq(start.yr, end.yr, by = 2),1)) +
      ylab(y.title) + xlab("Year") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 10),
            axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
            axis.text.y = element_text(size = 10))

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

