#' @title Rolling trend plots
#'
#' @description Demonstrates a distortion of the overall trend as a consequence of the leveraging effect of opening and
#' closing of monitoring sites with different magnitudes over the period examined using rolling trend plots.
#'
#' @param obs A data frame of ambient pollutant concentration data. Must contain the columns: site_code, date,
#' value. If 'pollutant' is a pollutant ratio, the data frames of the corresponding pollutants
#' must be supplied as a list of data frames in the order they are given in the ratio. E.g. for
#' \code{pollutant = "no2/nox"}, \code{obs = list(obs.no2, obs.nox)}.
#'
#' @param pollutant The pollutant of interest (character string). To calculate rolling change trend for a pollutant ratio,
#' separate the two pollutants with a forward slash e.g. \code{pollutant = "no2/nox"}.
#'
#' @param window.width The width of the moving window, n, over which the change in concentration is calculated (in years). This
#' can be a vector (to compare rolling and average trends over a range of window widths) or a numeric (to return the rolling
#' trend for a single window width).
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
#' @return A plot (or, if \code{window.width} is a vector, a list of plots labelled by the value of the window width) of the
#' rolling trends (left) and the average trend across the same data as was used in the rolling trends (i.e. filtered by the
#' \code{window.width}) (right).
#'
#' @import dplyr ggplot2 parallel doParallel foreach
#'
#' @examples
#' \dontrun{
#' rolling_trends(obs.nox, pollutant = "nox",
#' window.width = c(2, 5, 7, 10), stat = "median",
#' start.date = "2000-01-01", end.date = "2017-12-31",
#' data.capture = 90, smooth.method = "gam",
#' parallel = FALSE)
#' }
#'
#' @export

rolling_trends <- function(obs,
                           pollutant,
                           window.width = c(2, 3, 5, 7, 10, 12, 15),
                           stat = "median",
                           start.date = "2000-01-01",
                           end.date = "2017-12-31",
                           data.capture = 90,
                           smooth.method = "loess",
                           parallel = TRUE){


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


  # Function to plot rolling trend plots (compared to overall trend)
  rollingav_plot <- function(df, pollutant, smooth.method){

    rolling.yr <- df %>%
      ggplot(aes(x = date, y = av_value, group = as.factor(moving_window), color=trend)) +
      geom_point(show.legend = FALSE) +
      geom_smooth(show.legend = FALSE, se = FALSE, method = "gam") +
      scale_color_gradient2(low = "darkblue", mid = "darkviolet", high = "darkred", midpoint = 0) +
      xlab("Year") +
      ylab(openair::quickText(paste(pollutant, " concentration (ug m-3)"))) +
      scale_x_continuous(breaks = round(seq(lubridate::year(lubridate::ymd(start.date)), lubridate::year(lubridate::ymd(end.date)), by = 2),1)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 12),
           axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
           axis.text.y = element_text(size = 12))

    rolling.all <- df %>%
      ggplot(aes(x = date, y = av_value, color=trend)) +
      geom_point() +
      geom_smooth(method = smooth.method, show.legend = FALSE) +
      scale_color_gradient2(midpoint=0, low="darkblue", mid="darkviolet",
                            high="darkred", space ="Lab" ) +
      labs(x = "Year", colour = "Trend", fill = "Trend") +
      ylab(openair::quickText(paste(pollutant, " concentration (ug m-3)"))) +
      scale_x_continuous(breaks = round(seq(lubridate::year(lubridate::ymd(start.date)), lubridate::year(lubridate::ymd(end.date)), by = 2),1)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
            axis.text.y = element_text(size = 12))

    rolling.comp <- cowplot::plot_grid(rolling.yr, rolling.all, ncol = 2, rel_widths = c(1.65, 2))

    return(rolling.comp)

  }

  # Define the n moving windows
  define_moving_windows <- function( window.width, start.date, end.date){

    start <- lubridate::ymd(start.date)
    end <- lubridate::ymd(end.date)

    window1 <- seq.Date(start, end, by = "1 year") %>%
      .[-length(.)] # remove last element (will be the end of the final moving window)

    end.window <- window1[length(window1)]
    lubridate::year(end.window) <- lubridate::year(end) - (window.width-1) # set the start point of the final moving window as the time period end - (moving window width - 1)

    window1 <- window1[window1 <= end.window]

  }

  windows <- purrr::map(window.width, define_moving_windows, start.date = start.date, end.date = end.date)


  # Compute rolling trends for each moving window for each window width
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)

  rolling <- list()

  for(i in 1:length(window.width)){

      w <- windows[[i]]
      w.width <- window.width[i]

      if(parallel == TRUE){
        no_cores <- detectCores() - 1
        cl<-makeCluster(no_cores)
        registerDoParallel(cl)

        rolling.df <- foreach(d1 = w, .combine = "rbind", .packages = c("dplyr", "ggplot2"),
                              .export = c("data_capture", "average_data", "pollutant_ratio_data_capture", "calculate_pollutant_ratio")) %dopar%
          try(rolling_worker(d1, pollutant = pollutant, window.width = w.width, stat = stat))

        stopImplicitCluster()
      } else{
        rolling.df <- tryCatch({purrr::map_dfr(w, rolling_worker,
                                               pollutant = pollutant,
                                               window.width = w.width,
                                               stat = stat)
        }, error = function(e){
          data.frame("date" = as.POSIXct(character()),
                     "av_value" = numeric(),
                     "n" = numeric())
        })
      }

      rolling[[i]] <- rolling.df

    }


  rolling.names <- paste0("moving_window_width=", window.width)
  names(rolling) <- rolling.names

  # Plots
  rolling.plots <- purrr::map(rolling, rollingav_plot, pollutant = pollutant, smooth.method = smooth.method)

  return(rolling.plots)
}



