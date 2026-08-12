#' Dataset: employed_IT
#'
#' A dataset containing the quarterly time series of employed
#' people by age class in Italy in the time span 1992Q4-2024Q1.
#'
#' @docType data
#' 
#' @usage employed_IT
#' 
#' @format A multivariate time series with 126 rows and 10 columns:
#' \describe{
#'   \item{Y15.24}{Thousand of employed in the age class 15-24}
#'   \item{Y15.64}{Thousand of employed in the age class 15-64}
#'   \item{Y20.64}{Thousand of employed in the age class 20-64}
#'   \item{Y25.54}{Thousand of employed in the age class 25-54}
#'   \item{Y55.64}{Thousand of employed in the age class 55-64}
#'   \item{Y15.29}{Thousand of employed in the age class 15-29}
#'   \item{Y15.74}{Thousand of employed in the age class 15-74}
#'   \item{Y25.29}{Thousand of employed in the age class 25-29}
#'   \item{Y29.54}{Thousand of employed in the age class 29-54}
#'   \item{Y29.74}{Thousand of employed in the age class 29-74}
#' }
#' @source Eurostat
"employed_IT"

#' Hourly electricity load of Northern Italy and temperature of Milan
#'
#' The dataset start on 2023-01-01 01:00:00 and ends at 2025-09-30 23:00:00.
#' The hourly electricity load for Northern Italy has been downloaded from
#' the FTP of the GME (mercatoelettrico.org).
#' The temperature has been dowloaded form the Open-Meteo site
#' (archive-api.open-meteo.com) with latitude 45.4642 and longitude 9.1900.
#' The load is in MWh and the temperature is in Celsius degrees.
#'
#' @docType data
#'
#' @usage data(load_temp)
#'
#' @format An object of class \code{data.frame} with three columns:
#' #' \describe{
#'   \item{Datetime}{A "POSIXct" object with day and hour.}
#'   \item{Load}{The electricity load of Northern Italy in MWh.}
#'   \item{Temperature}{The temperature of Milan in Celsius degrees.}
#' }
#'
#' @keywords datasets
#'
#' @source \href{https://www.mercatoelettrico.org/}{GME},
#' \href{https://archive-api.open-meteo.com}{Open-Meteo}
#'
"load_temp"