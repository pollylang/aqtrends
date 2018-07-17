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

average_data <- function(df, avg.time = "year", statistic = "mean", sites = FALSE){


  if(nrow(df) > 1){

    df.out <- df %>%
      mutate(date = if(avg.time == "year"){
        lubridate::year(date)
      } else if(avg.time == "month"){
        format(as.Date(date), "%Y-%m")
      }) %>%
      group_by(date) %>%
      mutate(n = length(unique(site_code))) %>%
      ungroup()

    if(sites == FALSE){
      df.out <- df.out %>%
        group_by(date) %>%
        summarise(av_value = if(statistic == "median"){
          median(value, na.rm = TRUE)}
          else if(statistic == "mean"){
            av_value = mean(value, na.rm = TRUE)},
          n = unique(n))
    } else{
      df.out <- df.out %>%
        group_by(date, site_code) %>%
        summarise(av_value = if(statistic == "median"){
          median(value, na.rm = TRUE)}
          else if(statistic == "mean"){
            av_value = mean(value, na.rm = TRUE)},
          n = unique(n))
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
    bind_rows(obs.den) %>% # row-bind observations data frames of var1 and var2
    distinct(site_code, date, variable, .keep_all = TRUE) %>%
    filter(site_code %in% sites.int) %>% # filter to only include sites with both var1 and var2 measurements
    tidyr::spread(variable, value) %>%
    mutate(date = if(avg.time == "year"){ # convert date to year/year-month
        lubridate::year(date)
      } else if(avg.time == "month"){
        format(as.Date(date), "%Y-%m")
      }) %>%
    group_by(site_code) %>%
    filter(!(all(is.na(var1))) & !(all(is.na(var2)))) %>% # remove NA values
    group_by(date) %>%
    mutate(n = length(unique(site_code))) %>% # number of sites in each avg.time
    ungroup() %>%
    group_by(date,n, site_code) %>% # calculate linear model of var1 as function of var2 for each year
    do(model = try(lm(var1 ~ var2, data = .), silent = TRUE)) %>%
    mutate(av_value = tryCatch({summary(model)$coeff[2]
    }, error = function(e){
      NA}),
    n = unique(n)) %>% # extract gradient of linear model as var1/var2
    dplyr::select(-model)

  if(sites == FALSE){
    obs <- obs %>%
      group_by(date) %>%
      summarise(av_value = if(statistic == "median"){
        median(av_value, na.rm = TRUE)}
        else if(statistic == "mean"){
          av_value = mean(av_value, na.rm = TRUE)},
        n = unique(n))
  }

  return(obs)

}


#' Calculate average pollutant ratios (with data capture constraint)
#'
#' @param df.list List of data frames of observation data for each pollutant in the ratio. First element of list should be
#' data frame for the pollutant which is the numerator in the ratio (second element should be for the denominator pollutant).
#' Both data frames should contain the variables \code{date}, \code{site_code}, \code{value}.
#'
#' @param avg.time Character vector specifying the time period over which to average the data. Option are 'year' or
#' 'month' (to return the annual/monthly average)
#'
#' @param start,end A string in the format "YYYY-MM-DD" specifying the starting date and the ending date between
#' which to measure data capture.
#'
#' @param data.capture A value between 0-100 specifying the data capture threshold by which to filter the data
#'
#' @param statistic The statistic to use for averaging. Options are 'mean' or 'median'.
#'
#' @param sites Logical indicating whether averaging should be grouped by site (TRUE) or data should be averaged over all sites (FALSE).
#'
#' @return Data frame of annual/monthly average concentration data for the (filtered) pollutant ratio


pollutant_ratio_data_capture <- function(df.list, avg.time, start, end, data.capture, statistic = "mean", sites = FALSE){

  # Filter to only include dates within start-end date range, then filter by data capture over that period
  data.cap <- function(df, begin, finish, dc){
    out <- df %>%
      filter(date >= lubridate::ymd(begin), date <= lubridate::ymd(finish)) %>%
      data_capture(threshold = dc, start.date = begin, end.date = finish)

    return(out)
  }

  # Map data capture function over both obs data frames (numerator pollutant and denominator pollutant)
  dat <- purrr::map(df.list, data.cap, begin = start, finish = end, dc = data.capture)

  if(nrow(dat[[1]]) > 1 & nrow(dat[[2]]) > 1){

    # Calculate pollutant ratio for the filtered data
    out <- calculate_pollutant_ratio(dat, avg.time = avg.time, statistic = statistic, sites = sites)

  } else{

    out <- data.frame("date" = as.POSIXct(character()),
                         "av_value" = numeric(),
                         "n" = numeric())

  }


  return(out)

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

  plot <- ggplot(df, aes(x = date, y = av_value)) +
    geom_point() +
    xlab("Year") +
    ylab(openair::quickText(paste0(stringr::str_to_title(stat), " ", poll, " concentration (ugm-3)"))) +
    scale_x_continuous(breaks = round(seq(as.numeric(lubridate::year(start.date)), as.numeric(lubridate::year(end.date)), by = 2),1)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  if(sites == TRUE){
    plot <- plot +
      geom_smooth(size = 2, method = smooth.method) +
      geom_smooth(aes(color = site_code), se = FALSE, show.legend = FALSE, method = smooth.method, span = 1, size = 0.5)
  } else{
    plot <- plot +
      geom_smooth(method = smooth.method) +
      geom_text(aes(label = n), vjust = 1.3, size = 3)
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
                            data.capture = NULL,
                            smooth.method = NULL,
                            parallel = NULL){

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


  # 3. window.width
  if(!(is.null(window.width))){
    if(!(window.width >= 2 & window.width <= (lubridate::year(end.date) - lubridate::year(start.date)))){
      stop("Incorrect 'window.width' argument. The argument must be an integer between 2 and the range of the input dates in years.")
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


  # 6. data.capture
  if(!(is.null(data.capture))){
    if(!(is.numeric(data.capture) & data.capture >= 0 & data.capture <= 100)){
      stop("The 'data.capture' argument must be a numeric between 0 and 100 (%).")
    }
  }


  # 7. smooth.method
  if(!(is.null(smooth.method))){
    if(!(smooth.method %in% c("lm", "glm", "gam", "loess", "rlm"))){
      stop("The 'smooth.method' argument must be a character string with one of the following values: 'lm', 'loess', 'glm', 'gam', 'rlm'.")
    }
  }


  # 8. parallel
  if(!(is.null(parallel))){
    if(!(is.logical(parallel))){
      stop("The 'parallel' argument must be a logical.")
    }
  }


}

