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
#' @return A plot (or, if \code{window.width} is a vector, a list of plots labelled by the value of the window width) of the
#' rolling trends (left) and the average trend across the same data as was used in the rolling trends (i.e. filtered by the
#' \code{window.width}) (right).
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' rolling_trends(london_nox_data, pollutant = "nox",
#' window.width = c(2, 5, 7, 10), stat = "median",
#' start.date = "2000-01-01", end.date = "2017-12-31",
#' parallel = FALSE)
#' }
#'
#' @export

rolling_trends <- function(obs,
                           pollutant,
                           window.width = c(2, 3, 5, 7, 10, 12, 15),
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
                  parallel = parallel)


  # Function to plot rolling trend plots (compared to overall trend)
  rollingav_plot <- function(df, pollutant){

    rolling.yr <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = date, y = av_value, group = as.factor(moving_window), color=trend)) +
      ggplot2::geom_point(show.legend = FALSE) +
      ggplot2::geom_smooth(show.legend = FALSE, se = FALSE, method = "gam") +
      #ggplot2::scale_color_gradient2(low = "darkblue", mid = "darkviolet", high = "darkred", midpoint = 0) +
      viridis::scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.8) +
      ggplot2::labs(x = "Year", y = openair::quickText(paste0(pollutant, " concentration (ug m-3)"))) +
      ggplot2::theme_minimal()

    if(length(unique(df$date)) <= 5){
      smooth.method <- "gam"
    } else{
      smooth.method <- "loess"
    }

    rolling.all <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = date, y = av_value, color=trend)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = smooth.method,
                           show.legend = FALSE,
                           color = "black") +
      #scale_color_gradient2(midpoint=0, low="darkblue", mid="darkviolet", high="darkred", space ="Lab" ) +
      viridis::scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.8) +
      ggplot2::labs(x = "Year", y = openair::quickText(paste0(pollutant, " concentration (ug m-3)")),
                    color = "Slope") +
      ggplot2::theme_minimal()

    rolling.comp <- cowplot::plot_grid(rolling.yr, rolling.all, ncol = 2, rel_widths = c(1.65, 2))

    return(rolling.comp)

  }

  windows <- purrr::map(window.width, define_moving_windows, obs = obs,
                        start.date = start.date, end.date = end.date, avg.ts = avg.ts)


  # Compute rolling trends for each moving window for each window width
  rolling <- list()

  for(i in 1:length(window.width)){

      w <- windows[[i]]
      w.width <- window.width[i]

      if(parallel == TRUE){
        no_cores <- parallel::detectCores() - 1
        cl<-parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)

        rolling.df <- foreach::foreach(d1 = w, .combine = "rbind", .packages = c("dplyr", "ggplot2"),
                              .export = c("average_data", "calculate_pollutant_ratio")) %dopar%
          try(rolling_worker(d1, obs = obs, pollutant = pollutant, window.width = w.width, stat = stat, avg.ts = avg.ts,
                             verbose = verbose))

        doParallel::stopImplicitCluster()
      } else{
        rolling.df <- tryCatch({purrr::map_dfr(w, rolling_worker,
                                               obs = obs,
                                               pollutant = pollutant,
                                               window.width = w.width,
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

      rolling[[i]] <- rolling.df

    }

  rolling.names <- paste0("moving_window_width=", window.width)
  names(rolling) <- rolling.names

  # Plots
  rolling.plots <- purrr::map(rolling, rollingav_plot, pollutant = pollutant)

  return(rolling.plots)
}



