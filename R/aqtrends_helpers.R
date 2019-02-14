#' Average raw hourly data by year or month-year
#'
#' @param df Data frame containing the variables \code{date}, \code{site_code}, \code{value}.
#'
#' @param avg.time Character vector specifying the time period over which to average the data. Option are 'year' or
#' 'month' (to return the annual/monthly average)
#'
#' @param statistic The statistic to use for averaging. Options are 'mean' or 'median'.
#'
#' @param sites Logical indicating whether averaging should be grouped by site (TRUE) or data should be averaged over all sites (FALSE).
#'
#' @return Data frame of annual/monthly average concentration data

average_data <- function(df, avg.time = "year", statistic = "mean", sites = FALSE, date.format = FALSE){


  if(nrow(df) > 1){

    df.out <- df %>%
      dplyr::mutate(date = if(avg.time == "year"){
        lubridate::year(date)
      } else if(avg.time == "month"){
        format(lubridate::date(date), "%Y-%m")
      } else if(avg.time == "day"){
        lubridate::date(date)
      }) %>%
      dplyr::group_by(date) %>%
      dplyr::mutate(n = length(unique(site_code))) %>%
      dplyr::ungroup()

    if(sites == FALSE){
      df.out <- df.out %>%
        dplyr::group_by(date) %>%
        dplyr::summarise(av_value = if(statistic == "median"){
          median(value, na.rm = TRUE)}
          else if(statistic == "mean"){
            av_value = mean(value, na.rm = TRUE)},
          n = unique(n)) %>%
        dplyr::ungroup()
    } else{
      df.out <- df.out %>%
        dplyr::group_by(date, site_code) %>%
        dplyr::summarise(av_value = if(statistic == "median"){
          median(value, na.rm = TRUE)}
          else if(statistic == "mean"){
            av_value = mean(value, na.rm = TRUE)},
          n = unique(n)) %>%
        dplyr::ungroup()
    }

    if(date.format == TRUE){
      df.out <- df.out %>%
        dplyr::mutate(date = if(avg.time == "year"){
          paste0(date, "-01-01")
        } else if(avg.time == "month"){
          paste0(date, "-01")
        } else if(avg.time == "day"){
          date
        }) %>%
        dplyr::mutate(date = lubridate::ymd(date))
    }

    return(df.out)

  } else{

    df.out <- data.frame("date" = as.POSIXct(character()),
                         "av_value" = numeric(),
                         "n" = numeric())

    return(df.out)

  }

}


#' Calculate average pollutant ratios using linear model
#'
#' @param df.list List of data frames of observation data for each pollutant in the ratio. First element of list should be
#' data frame for the pollutant which is the numerator in the ratio (second element should be for the denominator pollutant).
#' Both data frames should contain the variables \code{date}, \code{site_code}, \code{value}.
#'
#' @param avg.time Character vector specifying the time period over which to average the data. Option are 'year' or
#' 'month' (to return the annual/monthly average)
#'
#' @param statistic The statistic to use for averaging. Options are 'mean' or 'median'.
#'
#' @param sites Logical indicating whether averaging should be grouped by site (TRUE) or data should be averaged over all sites (FALSE).
#'
#' @return Data frame of annual/monthly average concentration data for the pollutant ratio


calculate_pollutant_ratio <- function(df.list, avg.time, statistic, sites = FALSE){

  obs.num <- df.list[[1]]
  obs.den <- df.list[[2]]

  obs.num$variable <- "var1" # var1 = variable on numerator of ratio
  obs.den$variable <- "var2" # var2 = variable on denominator of ratio

  sites.int <- intersect(unique(obs.num$site_code), unique(obs.den$site_code)) # sites measuring both var1 and var2

  # Create var1/var2 variable (coefficient of the linear model of var1 hourly concentration as a function of
  # var2 hourly concentration)
  obs <- obs.num %>%
    dplyr::bind_rows(obs.den) %>% # row-bind observations data frames of var1 and var2
    dplyr::distinct(site_code, date, variable, .keep_all = TRUE) %>%
    dplyr::filter(site_code %in% sites.int) %>% # filter to only include sites with both var1 and var2 measurements
    tidyr::spread(variable, value) %>%
    dplyr::mutate(date = if(avg.time == "year"){ # convert date to year/year-month
        lubridate::year(date)
      } else if(avg.time == "month"){
        format(lubridate::date(date), "%Y-%m")
      } else if(avg.time == "week"){
        format(paste0(lubridate::year(date), "-", lubridate::week(date)))
      } else if(avg.time == "day"){
        lubridate::date(date)
      }) %>%
    dplyr::group_by(site_code) %>%
    dplyr::filter(!(all(is.na(var1))) & !(all(is.na(var2)))) %>% # remove NA values
    dplyr::group_by(date) %>%
    dplyr::mutate(n = length(unique(site_code))) %>% # number of sites in each avg.time
    dplyr::ungroup() %>%
    dplyr::group_by(date,n, site_code) %>% # calculate linear model of var1 as function of var2 for each year
    do(model = try(lm(var1 ~ var2, data = .), silent = TRUE)) %>%
    dplyr::mutate(av_value = tryCatch({summary(model)$coeff[2]
    }, error = function(e){
      NA}),
    n = unique(n)) %>% # extract gradient of linear model as var1/var2
    dplyr::select(-model)

  if(sites == FALSE){
    obs <- obs %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(av_value = if(statistic == "median"){
        median(av_value, na.rm = TRUE)}
        else if(statistic == "mean"){
          av_value = mean(av_value, na.rm = TRUE)},
        n = unique(n))
  }

  obs <- obs %>%
    dplyr::mutate(date = if(avg.time == "year"){
      paste0(date, "-01-01")
    } else if(avg.time == "month"){
      paste0(date, "-01")
    } else if(avg.time == "day"){
      date
    }) %>%
    dplyr::mutate(date = lubridate::ymd(date))

  return(obs)

}



#' Check that a character string is in correct date format
#'
#' @param date.string The character string of the date to be tested for whether the format is correct (\%Y-\%m-\%d)
#'
#' @return Logical value indicating whether the string is in the form YYYY-mm-dd (TRUE) or not (FALSE)


check_date_format <- function(date.string){

  # Checks that the input character string is in the format %Y-%m-%d
  out <- stringr::str_detect(date.string, "[0-2][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]")

  return(out)
}


#' Plot average trends in pollutant concentration
#'
#' @param df Data frame of average pollutant concentrations containing the variables \code{date}, \code{site_code}, \code{av_value}.
#'
#' @param sites Logical indicating whether to plot trends by site (TRUE) or average trend over all sites (FALSE).
#'
#' @param poll The pollutant for which to plot the trends.
#'
#' @param stat The statistic which was used for averaging (for the y axis label)
#'
#' @param smooth.method The smoothing method to use in the \code{geom_smooth} function when plotting trends.
#'
#' @param start.date,end.date A string in the format "YYYY-MM-DD" specifying the starting date and the ending date for the trends.
#'
#' @return A plot of the average trend in \code{poll} concentration over the period \code{start.date} - code{end.date} either
#' for each monitoring site or averaged over all monitoring sites


trends_plots_helper <- function(df, sites, poll, stat, start.date, end.date){

  if(nrow(df) <= 5){
    smooth.method <- "gam"
  } else{
    smooth.method <- "loess"
  }

  plot.colour <- viridis::inferno(1, begin=0.3, end=0.8)

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = av_value)) +
    ggplot2::geom_point() +
    ggplot2::xlab("Year") +
    ggplot2::ylab(openair::quickText(paste0(stringr::str_to_title(stat), " ", poll, " concentration (ugm-3)"))) +
    ggplot2::theme_minimal()

  if(sites == TRUE){
    plot <- plot +
      ggplot2::geom_smooth(size = 2,
                           method = smooth.method,
                           color = "black") +
      ggplot2::geom_smooth(ggplot2::aes(color = site_code), se = FALSE, show.legend = FALSE,
                           method = smooth.method,
                           span = 1, size = 0.5) +
      viridis::scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.8, discrete = TRUE)
  } else{
    plot <- plot +
      ggplot2::geom_smooth(method = "loess",
                           color =  plot.colour,
                           fill =  plot.colour,
                           alpha=0.25) +
      ggplot2::geom_text(ggplot2::aes(label = n), vjust = 1.3, size = 3)
  }

  return(plot)

}





#' Check function arguments
#'
#' @param obs,pollutant,window.width,stat,start.date,end.date,data.capture,smooth.method,parallel Pass the arguments to the
#' main function directly to the \code{check_arguments} function for checking.
#'
#' @return If one of the arguments violates the argument requirements, will return an informative error. Otherwise, the function
#' does nothing.

check_arguments <- function(obs = NULL,
                            pollutant = NULL,
                            window.width = NULL,
                            stat = NULL,
                            start.date = NULL,
                            end.date = NULL,
                            parallel = NULL,
                            avg.ts = NULL){

  # 1. obs
  if(!(is.null(obs))){
    if(!(stringr::str_detect(pollutant, "/"))){
      if(!(all(c("date", "site_code", "value") %in% colnames(obs)))){
        stop("Required columns missing from data frame. Please ensure 'obs' input has the following columns:'site_code', 'date', 'value'.")
      }
    } else{
      if(!(is.list(obs))){
        stop("For this pollutant, the input data frame 'obs' must be a list of data frames of observations data for each pollutant in the ratio (in the order of ratio notation).")
      }
      if(!(all(c("date", "site_code", "value") %in% colnames(obs[[1]]))) & all(c("date", "site_code", "value") %in% colnames(obs[[2]]))){
        stop("The input data frames in the 'obs' list must have the following columns: 'site_code', 'date', 'value'.")
      }
    }
  }


  # 2. pollutant
  if(!(is.null(pollutant))){
    if(!(is.character(pollutant))){
      stop("The 'pollutant' argument must be a character string specifying the name of the pollutant of interest.")
    }
  }


  # 4. stat
  if(!(is.null(stat))){
    if(!(stat %in% c("median", "mean"))){
      stop("The 'stat' argument must be one of the following: 'median', 'mean'.")
    }
  }


  # 5. start/end.date
  if(!(is.null(start.date)) & !(is.null(end.date))){
    if(!(check_date_format(start.date) & check_date_format(end.date))){ # pattern matching to check input dates are in correct format
      stop("The 'start.date' and 'end.date' arguments must be character strings in the format '%Y-%m-%d'.")
    }
  }


  # 8. parallel
  if(!(is.null(parallel))){
    if(!(is.logical(parallel))){
      stop("The 'parallel' argument must be a logical.")
    }
  }


  # 9. avg.ts
  if(!(is.null(avg.ts))){
    if(!(avg.ts %in% c("year", "month"))){
      stop("The 'avg.ts' argument must be one of the following: 'year', 'month'.")
    }
  }


}



#' Extract site codes of monitoring sites open throughout the duration of the moving window
#'
#' @param df Data frame of observations data from which to extract site codes
#'
#' @param resolution The time resolution for which to apply filtering. I.e. if \code{resolution} is "month", the function
#' will select and return monitoring sites measuring over a certain threshold of months during the moving window period.
#' Options: "month", "day"
#'
#' @param data.capture The data capture threshold to apply to filtering the sites. A value of 1 means that the site must
#' have measurements for every month/day within the moving window, while a value of 0 means no filtering takes place at all.
#' A value of 0.9 means that the site must have measurements for at least 90% of the months/days spanned by the moving window.
#'
#' @return A vector of site codes for monitoring sites open for the entire duration of the moving window (as specified by
#' the filtering parameters in the arguments).


sites_open_throughout_window <- function(df, resolution = "month", data.capture = 1){

  if(resolution == "month"){
    dates.window <- unique(format(as.Date(df$date), "%Y-%m"))
  } else if(resolution == "day"){
    dates.window <- unique(lubridate::date(df$date))
  }


  # Vector of sites open over entire duration of moving window (i.e. measuring in each month of window)
  open.sites <- df %>%
    dplyr::mutate(date = if(resolution == "month"){
      format(as.Date(date), "%Y-%m")
      } else if(resolution == "day"){
        lubridate::date(date)
      }) %>%
    dplyr::select(date, site_code) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site_code) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n >= (length(dates.window))*data.capture) %>%
    dplyr::select(site_code) %>%
    dplyr::pull()

  return(open.sites)

}



#' Return the site codes of the 'long term sites' (recorded measurements during every year of the period specified in the
#' function call with >= \code{data.capture} % data capture).
#'
#' @param df Data frame of observation data from which to get long term sites
#'
#' @param pollutant Pollutant of interest
#'
#' @param start,end Starting and ending year of period of interest (e.g. 2000, 2017)
#'
#' @param resolution Resolution at which to determine whether a site is long term. E.g. for a resolution of 'year', to qualify
#' as 'long term', a monitoring site must have measurements of the pollutant of interest in every year between \code{start}
#' and \code{end}. For a resolution of 'month', the site must have data measured during **each month** between \code{start}
#' and \code{end}.
#'
#' @param data.capture Data capture threshold to apply to 'long term sites' data. To qualify as a long term site,
#' the site must (in addition to having measurements in every year/month/day of the period) have >= \code{data.capture}
#' % data capture over the period of analysis (\code{start} - \code{end}). Default is 90% data capture.
#'
#' @return Vector of site codes of sites meeting the 'long term' criterion (measurements during every year/month/day of the
#' period of interest, and sufficient data capture over period of interest).

get_longterm_sites <- function(df, pollutant, start, end, resolution = "year", data.capture = 90){

  start.date <- paste0(start, "-01-01")
  end.date <- paste0(end, "-12-31")

  if(resolution == "year"){
    dates <- start:end
  } else if(resolution == "month"){
    dates <- seq.Date(as.Date(start.date), as.Date(end.date), by = "1 month") %>% format("%Y-%m")
  } else if(resolution == "day"){
    dates <- seq.Date(as.Date(start.date), as.Date(end.date), by = "1 day")
  } else if(resolution == "hour"){
    dates <- seq.POSIXt(lubridate::ymd(start.date, tz = "GMT"), lubridate::ymd(end.date, tz = "GMT"), by = "1 hour")
  } else{
    stop("The input time resolution is not supported. Supported resolutions are: 'day', 'month', 'year'.")
  }

  # Vector of sites open over entire duration of moving window (i.e. measuring in each month of window)

  longterm.sites <- df %>%
    dplyr::filter(date >= lubridate::ymd(start.date) & date <= lubridate::ymd(end.date),
                  variable == pollutant) %>%
    dplyr::mutate(date = if(resolution == "month"){
      format(as.Date(date), "%Y-%m")
    } else if(resolution == "day"){
      lubridate::date(date)
    } else if(resolution == "year"){
      lubridate::year(date)
    }) %>%
    dplyr::select(date, site_code) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site_code) %>%
    dplyr::mutate(longterm = dplyr::case_when(all(dates %in% date) ~ "yes",
                                              TRUE ~ "no")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(longterm == "yes") %>%
    dplyr::select(site_code) %>%
    dplyr::pull() %>%
    unique()

  longterm.data <- df %>%
    filter(site_code %in% longterm.sites) %>%
    data_capture(threshold = data.capture, start.date = start.date, end.date = end.date)

  longterm.sites <- longterm.sites %>% .[. %in% unique(longterm.data$site_code)]

  return(longterm.sites)

}



#' Define the n moving windows of the rolling regressions
#'
#' @param window.width Width of the moving window
#'
#' @param start.date,end.date Starting and ending dates of the entire period of analysis (e.g. "2000-01-01", "2017-12-31").
#'
#' @return Vector of starting dates of the n moving windows (in Date data type).

define_moving_windows <- function(window.width, obs, start.date, end.date, avg.ts){

  start <- lubridate::ymd(start.date)
  end <- lubridate::ymd(end.date)

  if(is.list(obs)){
    obs.all <- dplyr::bind_rows(obs)
  } else{
    obs.all <- obs
  }

  measurements <- obs.all %>%
    dplyr::filter(!(is.na(value))) %>%
    dplyr::mutate(year = lubridate::year(date)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(n = nrow(.)) %>%
    dplyr::ungroup() %>%
    dplyr::select(year, n) %>%
    dplyr::distinct() %>%
    dplyr::filter(n > 0) %>%
    dplyr::arrange(year)

  start <- lubridate::ymd(paste0(max(c(lubridate::year(start), min(measurements$year)), na.rm = T), "-01-01"))
  end <- lubridate::ymd(paste0(min(c(lubridate::year(end), max(measurements$year)), na.rm = T), "-12-31"))


  # Define the n moving windows
  if(avg.ts == "year"){
    interval <- "1 year"
  } else if(avg.ts == "month"){
    interval <- "1 month"
  }

  window1 <- seq.Date(start, end, by = interval)
  window1 <- window1[1:(length(window1) - (window.width-1))]  # make sure that start date of final moving window is the time period end - (moving window width - 1)

  return(window1)
}


#' Worker function for calculating rolling trend over a single window (i.e. one of the moving windows)
#'
#' @param d1 Date of the starting year of the moving window for which to calculate the rolling trend
#'
#' @param pollutant Pollutant for which to calculate rolling trend
#'
#' @param window.width Width of the moving window
#'
#' @param stat Statistic used to average data (mean, median)
#'
#' @import lubridate
#'
#' @return Dataframe of rolling trend for the moving window.

rolling_worker <- function(d1, obs, pollutant, window.width, stat, avg.ts, verbose){

  # Define start, and and range of the moving window
  time <- window.width - 1

  start <- d1
  end <- d1

  if(avg.ts == "year"){
    lubridate::year(end) <- lubridate::year(end) + time # create 'end' object to specify end of moving window where end = i + (n-1)
    end <- end %>% lubridate::year() %>% paste0("-12-31") %>% lubridate::ymd()
  } else if(avg.ts == "month"){
    end <- end %m+% months(time)
    end <- end %>%
      format("%Y-%m") %>%
      paste("-", lubridate::days_in_month(end), sep="") %>%
      lubridate::ymd()
  }

  range <- paste0(as.character(start), " - ", as.character(end)) # date range of window

  if(verbose == TRUE) print(paste0("Calculating rolling trend for the window: ", range))

  empty.df <- data.frame(date = numeric(),
                         av_value = numeric(),
                         trend_difference = numeric(),
                         n = integer(),
                         moving_window = character(),
                         window_width = numeric(),
                         trend = numeric())

  # Conditional branching - different averaging functions for pollutants and pollutant ratios
  if(!(stringr::str_detect(pollutant, "/"))){

    longterm.dat <- obs %>%
      dplyr::filter(date >= start, date <= end) # filter to include dates within moving window

    if(nrow(longterm.dat) > 0){

      open.sites <- sites_open_throughout_window(longterm.dat,
                                                 resolution = ifelse(avg.ts == "year", "month", "day"),
                                                 data.capture = 1)

      longterm.dat <- longterm.dat %>%
        dplyr::filter(site_code %in% open.sites) # filter to only include sites open for entire duration of window

      if(nrow(longterm.dat) == 0){

        return(empty.df)

      } else{
        longterm.dat <- longterm.dat %>%
          average_data(avg.time = avg.ts, statistic = stat, sites = FALSE, date.format = TRUE) # calculate annual average concentration
      }
    } else{

      return(empty.df)

    }


  } else{

    longterm.dat <- purrr::map(obs, function(x) x %>% dplyr::filter(date >= start, date <= end)) # filter to include dates within moving window)

    if(all(sapply(longterm.dat, nrow) > 0)){         # check that neither data frame in list is empty

      open.sites <- purrr::map(longterm.dat, sites_open_throughout_window,
                               resolution = ifelse(avg.ts == "year", "month", "day"),
                               data.capture = 1) %>%
                               {Reduce(intersect, .)}

      longterm.dat <- purrr::map(longterm.dat, function(x) x[x$site_code %in% open.sites, ])

      if(all(sapply(longterm.dat, nrow) == 0)){

        return(empty.df)

      } else{

        longterm.dat <- longterm.dat %>%
          calculate_pollutant_ratio(avg.time = avg.ts, statistic = stat)

      }

    } else{

      return(empty.df)

    }

  }

  longterm.dat <- longterm.dat %>%
    dplyr::mutate(moving_window = range, # moving window variable (date range and number of sites)
                  window_width = window.width)  # window width variable (n)

  # Rolling regression

  if(nrow(longterm.dat) > 1){
    # Convert dates to different scale (relative to start of window - i.e. start from 1 = start)
    if(avg.ts == "year"){
      start.yr <- lubridate::year(start)
      longterm.dat <- longterm.dat %>%
        dplyr::mutate(x = lubridate::year(date) - start.yr)
    } else if(avg.ts == "month"){
      seq.months <- seq.Date(start, end, "1 month")
      key <- data.frame("x" = 1:length(seq.months), "date" = seq.months)
      longterm.dat <- dplyr::left_join(longterm.dat, key, by = "date")
    }

    # Linear regression
    lm.trend <- lm(av_value ~ x, data = longterm.dat) # fit linear model (rolling regression)
    trend <- as.numeric(lm.trend$coefficients[2]) # extract linear model coefficient (beta_i)

    # Add regression coefficient to output data frame
    longterm.dat <- longterm.dat %>%
      dplyr::mutate(trend = rep(trend, times = nrow(longterm.dat))) %>% # append linear coefficient to data frame
      dplyr::mutate(date = as.Date(date))
  } else{
    trend <- NA  # if there are no long term sites, will return NA in the trend column
    longterm.dat <- longterm.dat %>%
      dplyr::mutate(trend = rep(trend, times = nrow(longterm.dat))) %>%
      dplyr::mutate(moving_window = as.character(moving_window),
                    trend = as.numeric(trend)) %>% # mutate columns into same classes as those with data (for purrr::map_dfr rowbind operation with other moving window data)
      dplyr::mutate(date = as.Date(date))
  }



  return(longterm.dat)

}



#' Filter observations data to only include sites with sufficient data capture
#'
#' @param data A data frame of observation data with the columns: site_code, date, value (i.e. the same as
#' the column names in the observations table of the air_quality_data database).
#'
#' @param threshold A value between 0-100 specifying the data capture threshold by which to filter the data
#'
#' @param start.date,end.date A string in the format "YYYY-MM-DD" specifying the starting date and the ending date between
#' which to measure data capture.
#'
#' @return Data frame of observation data filtered by data capture over the defined period
#'
#' @import dplyr


data_capture <- function(data, threshold, start.date = NULL, end.date = NULL) {

  # Colnames case
  colnames(data) <- tolower(colnames(data))

  # Check function arguments
  args <- c("site_code", "date", "value")

  if(!(all(args %in% colnames(data)))){
    stop("The input data frame must have the columns 'site code', 'date', 'value'.")
  }

  if(!is.numeric(threshold)){
    stop("The 'threshold' argument must be numeric.")
  }
  if(!(threshold >= 0 & threshold <= 100)) {
    stop("The 'threshold' argument must be between 0 and 100.")
  }

  # Check start.date and end.date
  if(!(is.null(start.date))){
    if(!is.character(start.date)){
      stop("start.date (if specified) must be a character string.")
    }
  }

  if(!(is.null(end.date))){
    if(!is.character(end.date)){
      stop("end.date (if specified) must be a character string.")
    }
  }

  # Convert start.date and end.date arguments to datetime format
  if(!(is.null(start.date))) start.date <- paste0(start.date, "00:00:00") %>% lubridate::ymd_hms()
  if(!(is.null(end.date))) end.date <- paste0(end.date, "00:00:00") %>% lubridate::ymd_hms()


  # Define find.time.interval functions from openair
  find.time.interval <- function(dates) {

    ## could have several sites, dates may be unordered
    ## find the most common time gap in all the data
    dates <- unique(dates) ## make sure they are unique

    # work out the most common time gap of unique, ordered dates
    id <- which.max(table(diff(as.numeric(unique(dates[order(dates)])))))
    seconds <- as.numeric(names(id))

    if ("POSIXt" %in% class(dates)) seconds <- paste(seconds, "sec")

    if (class(dates)[1] == "Date") {
      seconds <- seconds * 3600 * 24
      seconds <- paste(seconds, "sec")
    }

    seconds
  }

  # Convert data capture to a number between 0 and 1
  threshold <- threshold / 100

  # Split data by site code
  dat <- split(data, data$site_code)

  fun <- function(data){
    # if one line, just return
    if (nrow(data) < 2) return(data)

    ## time zone of data
    TZ <- attr(data$date, "tzone")
    if (is.null(TZ)) TZ <- "GMT" ## as it is on Windows for BST
    if(TZ == "") TZ <- "GMT"

    ## function to fill missing data gaps
    ## assume no missing data to begin with

    ## If start.date and end.date ARE NOT defined in function call, set them as the max and min dates
    ## of the date range in the data set (if they are defined, leave unchanged and pad data between the
    ## specified dates).

    ## pad out missing data

    if(is.null(start.date)) start.date <- min(data$date, na.rm = TRUE)
    if(is.null(end.date)) end.date <- max(data$date, na.rm = TRUE)


    ## interval in seconds
    interval <- find.time.interval(data$date)

    ## equivalent number of days, used to refine interval for month/year
    days <- as.numeric(strsplit(interval, split = " ")[[1]][1]) /
      24 / 3600

    ## find time interval of data
    if (class(data$date)[1] == "Date") {
      interval <- paste(days, "day")
    } else {
      ## this will be in seconds
      interval <- find.time.interval(data$date)
    }

    ## better interval, most common interval in a year
    if (days == 31) interval <- "month"
    if (days %in% c(365, 366)) interval <- "year"

    ## only pad if there are missing data
    if (length(unique(diff(data$date))) != 1L) {
      all.dates <- data.frame(date = seq(start.date, end.date, by = interval))
      data <- data %>% dplyr::full_join(all.dates, by = "date")

    }

    # add missing site codes
    data$site_code <- rep(unique(data$site_code[!is.na(data$site_code)]), times = nrow(data))

    # if the data capture for this site code is < the specified data capture, delete it
    if((sum(is.na(data$value))/length(data$value)) <= 1 - threshold) {data.out <- data
    } else data.out <- data.frame()

    return(data.out)

  }

  # Apply function to all data frames of different site codes and recombine
  data.out <- list()
  data.out <- purrr::map(dat, fun)
  data <- dplyr::bind_rows(data.out)

  # Remove NA values from 'values' column
  if(nrow(data) > 0) data <- data[!is.na(data$value), ]

  return(data)
}

