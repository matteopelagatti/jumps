# Replication script for the electricity-demand vs. temperature example.
#
# It documents how the two data sets shipped with the package were assembled
# from their original sources:
#
#   data(load_temp)     # full hourly series (load + temperature)
#   data(demand_temp)   # the filtered analysis sample used below
#
# You do NOT need to run this script to reproduce the analysis: both data sets
# are already available with data().  Section 4 of the introduction vignette
# reproduces the figures below with base graphics only.
#
# NOTE: the raw day-ahead market file MGPquant_20230101_20260407.csv is *not*
# shipped with the package (it is ~11 MB).  Download it from the GME FTP site
# (https://www.mercatoelettrico.org/) and place it in ./inst/scripts/ if you
# want to re-run the assembly from scratch.
#
# This script requires the packages httr, jsonlite, tidyverse, cowplot and
# scales, none of which is a dependency of `jumps`.

library(httr)
library(jsonlite)
library(tidyverse)


## ==========================================
## 1. LOAD TEMPERATURES (Open-Meteo)
## ==========================================

# Define the coordinates for the city of Milan
lat <- 45.4642
lon <- 9.1900
start_date <- "2023-01-01"
end_date   <- "2025-09-30"

# Build the URL for the Open-Meteo API
url_meteo <- sprintf(
  "https://archive-api.open-meteo.com/v1/archive?latitude=%s&longitude=%s&start_date=%s&end_date=%s&hourly=temperature_2m",
  lat, lon, start_date, end_date
)

# Call GET and parse the JSON
response_meteo <- GET(url_meteo)
dati_grezzi_meteo <- fromJSON(rawToChar(response_meteo$content))

# Create a clean dataframe
df_temperature <- data.frame(
  datetime = as.POSIXct(dati_grezzi_meteo$hourly$time,
                        format="%Y-%m-%dT%H:%M", tz="UTC"),
  temperatura_C = dati_grezzi_meteo$hourly$temperature_2m
)

head(df_temperature)

## ==========================================
## 2. LOAD DAY-AHEAD MARKET DATA FOR THE LOAD
## ==========================================


df_load <- read.csv("./inst/scripts/MGPquant_20230101_20260407.csv")

df_load$Data <- as.Date(as.character(df_load$Data), format = "%Y%m%d")
df_load <-
  df_load |>
  mutate(datetime = lubridate::ymd_h(paste(Data, Ora)))

df_load |>
  select(datetime, MWh = NORD_ACQUISTI) |> # NORTHEN ITALY
  inner_join(df_temperature) ->
  df

df <- df |> rename(Celsius = temperatura_C)

# FINAL SELECTION
df |>
  filter(wday(datetime, week_start = 1) %in% (2:4), # TUE, WED, THU
         month(datetime) %in% c(2:3, 5:7, 9:11),    # NO VACATION MONTHS
         hour(datetime) %in% c(10:15),              # CENTRAL HOURS
         MWh > 15000, Celsius < 35) |>              # LIMIT EXTREMES
  select(Celsius, MWh) ->
  dt


# PASS SMOOTHING SPLINES WITH JUMPS
library(jumps)
set.seed(202609)

result_mle <- auto_ssj_mle(y = dt$MWh, x = dt$Celsius,
                           ic = "ebic", ebic_xi = 0.3)

result_atw <- auto_ssj_atw_cv(y = dt$MWh, x = dt$Celsius,
                              lambda_grid = ssj_lambda_grid(dt$Celsius,
                                                            n_lambda = 15L,
                                                            from = 1e-3,
                                                            to = 1e+1),
                              ic = "ebic", ebic_xi = 0.3)


sw <- sqrt(mean(tapply(dt$MWh, dt$Celsius,
                       function(z) if (length(z) > 1) var(z) else NA), na.rm = TRUE))
result_cssd <- auto_ssj_cssd(y = dt$MWh, x = dt$Celsius, delta = sw)

p_mle <- plot(result_mle, main = "SSJ-MLE", xlab = "Temperature C°", ylab = "Load (MWh)")
ggplot2::ggsave("demand_vs_temperature_mle.pdf", p_mle, width = 4.5, height = 4.5, dpi = 300, bg = "white")

p_atw <- plot(result_atw, main = "SSJ-ATW (CV for lambda)", xlab = "Temperature C°", ylab = "Load (MWh)")
ggplot2::ggsave("demand_vs_temperature_atw.pdf", p_atw, width = 4.5, height = 4.5, dpi = 300, bg = "white")

plot(result_cssd, main = "CSSD", xlab = "Temperature C°", ylab = "Load (MWh)")

# Sensitivity to the choice of xi
vxi <- seq(0.10, 0.5, 0.05)
res_mle <- vector("list", length(vxi))
res_atw <- vector("list", length(vxi))
for (i in seq_along(vxi)) {
  set.seed(202609)
  res_mle[[i]] <- auto_ssj_mle(y = dt$MWh, x = dt$Celsius,
                             ic = "ebic", ebic_xi = vxi[i])
  
  res_atw[[i]] <- auto_ssj_atw_cv(y = dt$MWh, x = dt$Celsius,
                                lambda_grid = ssj_lambda_grid(dt$Celsius,
                                                              n_lambda = 15L,
                                                              from = 1e-3,
                                                              to = 1e+1),
                                ic = "ebic", ebic_xi = vxi[i])
}

names(res_mle) <- vxi
names(res_atw) <- vxi

sapply(res_mle, function(x) ssj_breaks(x) |> nrow())
sapply(res_atw, function(x) ssj_breaks(x) |> nrow())

df_warps <- rbind(
  cbind(method = "MLE", purrr::map_dfr(res_mle, ssj_breaks, .id = "xi")),
  cbind(method = "ATW", purrr::map_dfr(res_atw, ssj_breaks, .id = "xi"))
)

df_warps |>
  ggplot(aes(y = x, x = xi, size = warp, fill = method)) +
  geom_point(shape = 21, color = "black", stroke = 1) +
  facet_grid(~method) +
  xlab(expression(xi)) +
  ylab("Temperature (°C)") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(expression("Warp place and size as a function of" ~ xi))
ggsave("warp_place_size.pdf", width = 8, height = 4)

