% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_flux_bias.R
\name{site_flux_bias}
\alias{site_flux_bias}
\title{The difference in concentration between opening and closing sites as a function of year}
\usage{
site_flux_bias(obs, pollutant, stat = "median")
}
\arguments{
\item{obs}{A data frame of ambient pollutant concentration data. Must contain the columns: site_code, date,
value. If 'pollutant' is a pollutant ratio, the data frames of the corresponding pollutants
must be supplied as a list of data frames in the order they are given in the ratio. E.g. for
\code{pollutant = "no2/nox"}, \code{obs = list(obs.no2, obs.nox)}.}

\item{pollutant}{The pollutant of interest (character string). To calculate rolling change trend for a pollutant ratio,
separate the two pollutants with a forward slash e.g. \code{pollutant = "no2/nox"}.}

\item{stat}{The metric (character string) used to average the ambient concentration data by year. Options: "median", "mean".}
}
\value{
A list of plots:
\itemize{
    \item \code{difference} - a plot of the difference between the average concentration of opening sites and closing
    sites as a function of year.
    \item \code{cumulative_difference} - the cumulative sum of differences in average concentration of opening sites and
    closing sites as a function of year.
    \item \code{weighted_cumulative_difference} - the cumulative sum of differences in average concentration of opening and
    closing sites weighted by the number of opening/closing sites as a function of year.
    \item \code{compare} - a side-by-side comparison of the \code{difference}, \code{cumulative_difference} and
    \code{weighted_cumulative_difference} plots.
}
}
\description{
Detects a difference in average concentration of the pollutant of interest between sites opening
and closing in a given year (i.e. a bias towards more/less polluted locations) and any change in this bias over
time.
}
\examples{
\dontrun{
site_flux_bias(london_nox_data, pollutant = "nox", stat = "median")
}

}
