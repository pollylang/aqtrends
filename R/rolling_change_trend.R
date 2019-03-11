#' @title Rolling change trends (extract true trend)
#'
#' @description Removal of distorting effect of site movement to reveal underlying trend by using changes in
#' concentration as a function of year as a proxy for the average trend. This method retains information about the
#' shape of the trend, while ignoring differences in magnitude, thus removing the leveraging effect of opening and
#' closing sites with extreme magnitudes.
#'
#' @param obs A data frame of ambient pollutant concentration data. Must contain the columns: site_code, date,
#' value. If 'pollutant' is a pollutant ratio, the data frames of the corresponding pollutants
#' must be supplied as a list of data frames in the order they are given in the ratio. E.g. for
#' \code{pollutant = "no2/nox"}, \code{obs = list(obs.no2, obs.nox)}.
#'
#' @param pollutant The pollutant of interest (character string). To calculate rolling change trend for a pollutant ratio,
#' separate the two pollutants with a forward slash e.g. \code{pollutant = "no2/nox"}.
#'
#' @param window.width The width of the moving window, n, over which the change in concentration is calculated (in years).
#'
#' @param avg.ts The resolution to which to average each time series, upon which the rolling regression is carried out.
#' For example, setting \code{avg.ts = "day"} means the rolling regression will be carried out on the daily average concentrations
#' from each time series (monitoring site). Options are: "year", "month", "week", and "day".
#'
#' @param stat The metric (character string) used to average the ambient concentration data by year. Options: "median", "mean".
#'
#' @param start.date,end.date The starting and ending dates (character string) of the period of interest over which
#' to calculate and plot the change trend.
#'
#' @param parallel Logical indicating whether the rolling changes should be computed in parallel. If \code{TRUE}, the
#' parallelisation will be implemented using the \code{foreach} function. The number of cores used will be the total
#' number of cores - 1.
#'
#' @param verbose Logical indicating whether to print the date range of the rolling window over which the calculation is
#' being applied.
#'
#' @return A plot of the rolling change trend.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' rolling_change_trend(london_nox_data,
#' pollutant = "nox",
#' window.width = 3,
#' avg.ts = "year",
#' stat = "median",
#' start.date = "2000-01-01", end.date = "2017-12-31")
#' }
#'
#' @export


rolling_change_trend <- function(obs,
                                 pollutant,
                                 window.width,
                                 avg.ts = "year",
                                 stat = "median",
                                 start.date = "2000-01-01",
                                 end.date = "2017-12-31",
                                 parallel = FALSE,
                                 verbose = FALSE){


  ## Check arguments
  check_arguments(obs = obs,
                  pollutant = pollutant,
                  window.width = window.width,
                  stat = stat,
                  start.date = start.date,
                  end.date = end.date,
                  parallel = parallel,
                  avg.ts = avg.ts)


  # Function to calculate data frame with row for difference between initial conc and trend for each moving window
  trend_difference <- function(rolling.x){

    # Initialise delta y_1 (concentration change value in the first year = raw annual average concentration in first year)
    first.row <- rolling.x[1, ] %>% dplyr::mutate(trend = av_value) # copy the first row to the data frame (initialise delta y1)

    rolling.x <- rolling.x %>%
      tibble::add_row(date = first.row$date, av_value = first.row$av_value, n = first.row$n,
              moving_window = paste0(first.row$moving_window, ".1"), window_width = first.row$window_width,
              trend = first.row$trend, .before = 1) %>% # add first row (initialise - delta y1)
      dplyr::mutate(moving_window = ifelse(moving_window == first.row$moving_window,
                                           paste0(first.row$moving_window, ".2"),
                                           moving_window)) %>%
      dplyr::group_by(moving_window) %>%
      dplyr::summarise(trend = unique(trend),
                       n = unique(n),
                       date = (min(date) + floor((max(date)-min(date))/2)),
                       window_width = unique(window_width)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(date) %>%
      dplyr::mutate(trend_difference = ifelse(moving_window == paste0(first.row$moving_window, ".1"),
                                              trend, NA))

    # Calculate concentration change: delta y_i = delta y_{i-1} + beta_i
    for(i in 2:nrow(rolling.x)){

      rolling.x[i, "trend_difference"] <- rolling.x[i-1, "trend_difference"] + rolling.x[i, "trend"]

    }

    return(rolling.x)
  }


  # Plot rolling change trend (concentration change as a function of the year)
  period_plot <- function(df){

    plot.colour <- viridis::inferno(1, begin=0.3, end=0.8)

    p <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = date, y = trend_difference)) +
      ggplot2::geom_point(color = plot.colour) +
      ggplot2::geom_smooth(method = "loess",
                  color =  plot.colour,
                  fill =  plot.colour,
                  alpha=0.25) +
      ggplot2::xlab("Year") +
      ggplot2::ylab(openair::quickText(paste0(pollutant, " concentration (ug m-3)"))) +
      ggplot2::geom_text(ggplot2::aes(label = n), vjust = 1.3, size = 3) +
      ggplot2::theme_minimal()

    return(p)

  }

  # Define the start and end dates
  window1 <- define_moving_windows(window.width, obs, start.date, end.date, avg.ts)

  # For each moving window, compute rolling regression
  if(parallel == TRUE){
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)

    rolling.df <- foreach::foreach(d1 = window1, .combine = "rbind", .packages = c("dplyr", "ggplot2", "lubridate"),
                          .export = c("average_data", "sites_open_throughout_window", "calculate_pollutant_ratio", "rolling_worker")) %dopar%
      try(rolling_worker(d1, obs = obs, pollutant = pollutant, window.width = window.width, stat = stat, avg.ts = avg.ts,
                         verbose = verbose))

    doParallel::stopImplicitCluster()
  } else{
    rolling.df <- tryCatch({purrr::map_dfr(window1, rolling_worker,
                                           obs = obs,
                                           pollutant = pollutant,
                                           window.width = window.width,
                                           stat = stat, avg.ts = avg.ts,
                                           verbose = verbose)
    }, error = function(e){
      print("Error message: ")
      print(e)
      data.frame("date" = as.POSIXct(character()),
                 "av_value" = numeric(),
                 "n" = numeric())
    })
  }


  if(nrow(rolling.df) < 1){
    print("Error in rolling trend calculation. Possibly insufficient data.")
    return(NULL)
  } else{

    # Calculate concentration change for each year (moving window)
    rolling.out <- trend_difference(rolling.df)

    # Correct the 'date' variable to always be in the form YYYY-01-01 (where YYYY is the middle year of the moving window)
    rolling.out <- rolling.out %>%
      mutate(year = case_when(lubridate::month(date) == 1 ~ lubridate::year(date),
                              lubridate::month(date) == 12 ~ lubridate::year(date)+1,
                              TRUE ~ lubridate::year(date))) %>%
      mutate(date = lubridate::ymd(paste0(year, "-01-01")))

    # Plot rolling change trend
    plots <- period_plot(rolling.out)

    return(plots)
  }

}

