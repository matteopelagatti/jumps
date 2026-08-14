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
ss1 <- auto_ssj_mle(y = dt$MWh, x = dt$Celsius,
                    ic = "ebic", ebic_xi = 0.1)

# BASIC PLOT
plot(ss1$x, ss1$smoothed_level, type = "l",
     lwd = 2, main = "Electricity demand vs. temperature",
     xlab = "Temperature (C°)",
     ylab = "Demand (MWh)")
points(x = dt$Celsius, y = dt$MWh,
       pch = 16, col = "blue", cex = 0.5)
abline(v  = ss1$x[-1][ss1$gamma > .2],
       lty = 2, lwd = 2)
grid()


## ==========================================
## 4. FINAL PLOT (ggplot2)
## ==========================================
# needs cowplot, only for plot_grid(), which aligns the two panels vertically

jump_thr <- 0.2   # threshold on gamma above which a jump is declared

# fitted level (ss1$x is sorted, ties included) with its +/- 2 s.e. band
fit_df <- tibble(
  Celsius = ss1$x,
  level   = ss1$smoothed_level,
  se      = sqrt(ss1$var_smoothed_level)
)

# gamma[i] refers to the interval ending at x[i + 1] -> align it with ss1$x[-1]
gamma_df <- tibble(
  Celsius = ss1$x[-1],
  gamma   = ss1$gamma,
  is_jump = gamma > jump_thr
)

jump_x <- gamma_df$Celsius[gamma_df$is_jump]

col_pts  <- "#4C78A8"
col_fit  <- "#D55E00"
col_jump <- "grey30"

# common horizontal scale, so the two panels line up exactly
x_lim   <- range(ss1$x) + c(-1, 1) * 0.02 * diff(range(ss1$x))
x_scale <- scale_x_continuous(limits = x_lim, expand = c(0, 0),
                              breaks = scales::breaks_width(5))

theme_paper <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey88"),
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "grey35", size = 10),
        plot.caption     = element_text(colour = "grey45", size = 8.5, hjust = 0),
        plot.margin      = margin(5.5, 10, 5.5, 10))

# ---- upper panel: data, fitted spline, jump locations --------------------
p_fit <- ggplot() +
  geom_vline(xintercept = jump_x, linetype = "dashed",
             linewidth = 0.4, colour = col_jump) +
  geom_point(data = dt, aes(Celsius, MWh),
             colour = col_pts, alpha = 0.3, size = 1.1) +
  geom_ribbon(data = fit_df,
              aes(Celsius, ymin = level - 2 * se, ymax = level + 2 * se),
              fill = col_fit, alpha = 0.3) +
  geom_path(data = fit_df, aes(Celsius, level),   # path, not line: keeps the
            colour = col_fit, linewidth = 1) +    # vertical steps at the jumps
  x_scale +
  scale_y_continuous(labels = scales::label_comma()) +
  # zoom in (coord_*, not scale limits: the few points below are just hidden)
  coord_cartesian(ylim = c(18000, NA)) +
  labs(title    = "Electricity demand vs. temperature",
       subtitle = "Northern Italy, Tue–Thu, hours 10–15",
       x = NULL, y = "Demand (MWh)") +
  theme_paper +
  theme(axis.text.x = element_blank())

# ---- lower panel: estimated jump strengths -------------------------------
p_gamma <- ggplot(subset(gamma_df, is_jump), aes(Celsius, gamma)) +
  geom_vline(xintercept = jump_x, linetype = "dashed",
             linewidth = 0.4, colour = col_jump) +
  geom_hline(yintercept = jump_thr, linetype = "dotted",
             linewidth = 0.4, colour = col_jump) +
  geom_segment(aes(xend = Celsius, yend = 0),
               colour = col_fit, linewidth = 0.7) +
  geom_point(colour = col_fit, size = 1.6) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Temperature (°C)",
       y = expression("Jump strength " * gamma)#,
       # caption = paste0("Dashed lines mark the ", length(jump_x),
       #                  " estimated jumps (jump strength above ", jump_thr,
       #                  "); the shaded band is the fitted level ± 2 standard errors.")
       ) +
  theme_paper

fig <- cowplot::plot_grid(p_fit, p_gamma, ncol = 1, align = "v", axis = "lr",
                          rel_heights = c(2, 1))
fig

ggsave("demand_vs_temperature.pdf", fig, width = 6, height = 6, dpi = 300, bg = "white")


