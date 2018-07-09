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
#' @param stat The metric (character string) used to average the ambient concentration data by year. Options: "median", "mean".
#'
#' @param start.date,end.date The starting and ending dates (character string) of the period of interest over which
#' to calculate and plot the change trend.
#'
#' @param data.capture The data capture threshold used to filter data within each moving window. During calculation of the
#' rolling regression, only data from sites with data capture >= \code{data.capture} over the period encapsulated by the
#' moving window will be included in the rolling regression.
#'
#' @param smooth.method The smoothing method to use in the \code{geom_smooth} function when plotting trends.
#'
#' @param parallel Logical indicating whether the rolling changes should be computed in parallel. If \code{TRUE}, the
#' parallelisation will be implemented using the \code{foreach} function. The number of cores used will be the total
#' number of cores - 1.
#'
#' @return A plot of the proxy trend (change in concentration over an N year period as a function of time).
#'
#' @import dplyr ggplot2 parallel doParallel foreach
#'
#' @examples
#' \dontrun{
#' rolling_change_trend(obs.nox,
#' pollutant = "nox",
#' window.width = 3,
#' stat = "median",
#' start.date = "2000-01-01", end.date = "2017-12-31")
#' }
#'
#' @export


rolling_change_trend <- function(obs,
                             pollutant,
                             window.width,
                             stat = "median",
                             start.date = "2000-01-01",
                             end.date = "2017-12-31",
                             data.capture = 90,
                             smooth.method = "loess",
                             parallel = FALSE){

  ## Check arguments
  check_arguments(obs = obs,
                  pollutant = pollutant,
                  window.width = window.width,
                  stat = stat,
                  start.date = start.date,
                  end.date = end.date,
                  data.capture = data.capture,
                  smooth.method = smooth.method,
                  parallel = parallel)


  # Worker function for calculating rolling trends over n year periods
  rolling_worker <- function(d1, pollutant, window.width, stat){

    # Define start, and and range of the moving window
    time <- window.width - 1

    start <- d1
    end <- d1

    lubridate::year(end) <- lubridate::year(end) + time # create 'end' object to specify end of moving window where end = i + (n-1)
    end <- end %>% lubridate::year() %>% paste0("-12-31") %>% lubridate::ymd()

    range <- paste0(as.character(start), " - ", as.character(end)) # date range of window

    print(paste0("Calculating rolling trend for the window: ", range))

    # Conditional branching - different averaging functions for pollutants and pollutant ratios
    if(!(stringr::str_detect(pollutant, "/"))){

      dat <- obs %>% filter(date >= start, date <= end) # filter to include dates within moving window

      longterm.dat <- data_capture(dat, threshold = data.capture, start.date = as.character(start),
                                   end.date = as.character(end)) %>% # filter to include only sites with sufficient data capture over window
        average_data(avg.time = "year", statistic = stat) %>% # calculate annual average concentration
        mutate(moving_window = paste0(range, " (", n, ")"), # moving window variable (date range and number of sites)
               window_width = window.width) # window width variable (n)

    } else{

      longterm.dat <- pollutant_ratio_data_capture(obs, avg.time = "year", as.character(start),
                                                   as.character(end), data.capture, statistic = stat) %>% # helper function filters by moving window range and data capture, then averages data
        mutate(moving_window = paste0(range, " (", n, ")"), # moving window variable (date range and number of sites)
               window_width = window.width) # window width variable (n)

    }


    # Rolling regression
    if(nrow(longterm.dat) > 1){
      lm.trend <- lm(av_value ~ date, data = longterm.dat) # fit linear model (rolling regression)
      trend <- as.numeric(lm.trend$coefficients[2]) # extract linear model coefficient (beta_i)
    } else{
      trend <- NA # if there are no long term sites, will return NA in the trend column
    }

    longterm.dat <- longterm.dat %>%
      mutate(trend = rep(trend, times = nrow(longterm.dat))) # append linear coefficient to data frame

    return(longterm.dat)

  }

  # Function to calculate data frame with row for difference between initial conc and trend for each moving window
  trend_difference <- function(rolling.x){

    # Initialise delta y_1 (concentration change value in the first year = raw annual average concentration in first year)
    first.row <- rolling.x[1, ] %>% mutate(trend = av_value) # copy the first row to the data frame (initialise delta y1)

    rolling.x <- rolling.x %>%
      add_row(date = first.row$date, av_value = first.row$av_value, n = first.row$n,
              moving_window = paste0(first.row$moving_window, ".1"), window_width = first.row$window_width,
              trend = first.row$trend, .before = 1) %>% # add first row (initialise - delta y1)
      group_by(moving_window) %>%
      summarise(trend = unique(trend),
                n = unique(n),
                date = median(date, na.rm = TRUE),
                window_width = unique(window_width)) %>%
      ungroup() %>%
      arrange(date) %>%
      mutate(trend_difference = ifelse(moving_window == paste0(first.row$moving_window, ".1"),
                                       trend, NA))

    # Calculate concentration change: delta y_i = delta y_{i-1} + beta_i
    for(i in 2:nrow(rolling.x)){

      rolling.x[i, "trend_difference"] <- rolling.x[i-1, "trend_difference"] + rolling.x[i, "trend"]

    }

    return(rolling.x)
  }


  # Plot rolling change trend (concentration change as a function of the year)
  period_plot <- function(df){

    p <- df %>%
      ggplot(aes(x = date, y = trend_difference)) +
      geom_point() +
      geom_smooth(method = smooth.method) +
      xlab("Year") +
      ylab(openair::quickText(paste("Change in ", pollutant, " concentration (ug m-3)"))) +
      geom_text(aes(label = n), vjust = 1.3, size = 3) +
      scale_x_continuous(breaks = round(seq(as.numeric(lubridate::year(start.date)),
                                            as.numeric(lubridate::year(end.date)), by = 2),1)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 12),
            #axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
            axis.text = element_text(size = 12))

    return(p)

  }


  # Define the n moving windows
  start <- lubridate::ymd(start.date)
  end <- lubridate::ymd(end.date)

  window1 <- seq.Date(start, end, by = "1 year") %>%
    .[-length(.)] # remove last element (will be the end of the final moving window)

  end.window <- window1[length(window1)]
  lubridate::year(end.window) <- lubridate::year(end) - (window.width-1) # set the start point of the final moving window as the time period end - (moving window width - 1)

  window1 <- window1[window1 <= end.window]


  # For each moving window, compute rolling regression
  if(parallel == TRUE){
    no_cores <- detectCores() - 1
    cl<-makeCluster(no_cores)
    registerDoParallel(cl)

    rolling.df <- foreach(d1 = window1, .combine = "rbind", .packages = c("dplyr", "ggplot2"),
                          .export = c("data_capture", "average_data", "pollutant_ratio_data_capture", "calculate_pollutant_ratio")) %dopar%
      try(rolling_worker(d1, pollutant = pollutant, window.width = window.width, stat = stat))

    stopImplicitCluster()
  } else{
    rolling.df <- tryCatch({purrr::map_dfr(window1, rolling_worker,
                                 pollutant = pollutant,
                                 window.width = window.width,
                                 stat = stat)
    }, error = function(e){
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

    # Plot rolling change trend
    plots <- period_plot(rolling.out)

    return(plots)
  }


}
