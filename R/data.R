#' Employed people in Italy by age class
#'
#' Quarterly time series of the number of employed people by age class in
#' Italy, over the period 1992Q4-2024Q1.
#'
#' @docType data
#'
#' @usage employed_IT
#'
#' @format A multivariate time series (class \code{ts}) with 126 rows and 10
#'   columns:
#' \describe{
#'   \item{Y15.24}{Thousands of employed in the age class 15-24.}
#'   \item{Y15.64}{Thousands of employed in the age class 15-64.}
#'   \item{Y20.64}{Thousands of employed in the age class 20-64.}
#'   \item{Y25.54}{Thousands of employed in the age class 25-54.}
#'   \item{Y55.64}{Thousands of employed in the age class 55-64.}
#'   \item{Y15.29}{Thousands of employed in the age class 15-29.}
#'   \item{Y15.74}{Thousands of employed in the age class 15-74.}
#'   \item{Y25.29}{Thousands of employed in the age class 25-29.}
#'   \item{Y29.54}{Thousands of employed in the age class 29-54.}
#'   \item{Y29.74}{Thousands of employed in the age class 29-74.}
#' }
#'
#' @keywords datasets
#'
#' @source Eurostat.
"employed_IT"

#' Hourly electricity load of Northern Italy and temperature of Milan
#'
#' Hourly observations from 2023-01-01 01:00:00 to 2025-09-30 23:00:00. The
#' electricity load for Northern Italy was downloaded from the FTP site of the
#' GME (mercatoelettrico.org). The temperature was downloaded from the
#' Open-Meteo archive (latitude 45.4642, longitude 9.1900).
#'
#' @docType data
#'
#' @usage data(load_temp)
#'
#' @format A \code{data.frame} with 24092 rows and three columns:
#' \describe{
#'   \item{Datetime}{A \code{POSIXct} object with day and hour.}
#'   \item{Load}{The electricity load of Northern Italy, in MWh.}
#'   \item{Temperature}{The temperature of Milan, in degrees Celsius.}
#' }
#'
#' @keywords datasets
#'
#' @source \href{https://www.mercatoelettrico.org/}{GME},
#'   \href{https://open-meteo.com/}{Open-Meteo}.
"load_temp"

#' Electricity demand and temperature: the analysis sample
#'
#' The homogeneous sub-sample of \code{\link{load_temp}} used in the electricity
#' demand versus temperature application. It keeps only Tuesdays to Thursdays,
#' the non-vacation months (February-March, May-July, September-November) and
#' the central hours of the working day (10-15), and drops loads below
#' 15000 MWh and temperatures of 35 degrees Celsius or above. These filters
#' remove the weekday, holiday and time-of-day effects that would otherwise
#' confound the relationship between demand and temperature.
#'
#' The script that builds this data set from \code{load_temp} is in
#' \code{data-raw/demand_temp.R} of the package sources;
#' \code{inst/scripts/load_temperature.R} shows how \code{load_temp} itself is
#' assembled from the raw GME and Open-Meteo sources.
#'
#' @docType data
#'
#' @usage data(demand_temp)
#'
#' @format A \code{data.frame} with 1691 rows and two columns:
#' \describe{
#'   \item{Temperature}{The temperature of Milan, in degrees Celsius.}
#'   \item{Load}{The electricity load of Northern Italy, in MWh.}
#' }
#'
#' @keywords datasets
#'
#' @seealso \code{\link{load_temp}} for the full hourly series.
#'
#' @source \href{https://www.mercatoelettrico.org/}{GME},
#'   \href{https://open-meteo.com/}{Open-Meteo}.
"demand_temp"
