# Builds data/demand_temp.rda, the analysis sample used in the electricity
# demand vs. temperature application of inst/scripts/load_temperature.R.
#
# Source: the `load_temp` data set shipped with the package, which is the
# hourly Northern-Italy electricity load (GME day-ahead market, column
# NORD_ACQUISTI of MGPquant_20230101_20260407.csv) merged with the Milan
# temperature from the Open-Meteo archive.  See inst/scripts/load_temperature.R
# for how `load_temp` itself is assembled from the raw sources.
#
# The filters below isolate a homogeneous sample in which the load-temperature
# relationship is not confounded by weekday, holiday or time-of-day effects:
#   * Tuesday to Thursday only (no weekend or Monday/Friday effects);
#   * February-March, May-July and September-November (no vacation months);
#   * hours 10-15 (central hours of the working day);
#   * loads above 15000 MWh and temperatures below 35 C (drop extremes).

load("data/load_temp.rda")

wd <- as.integer(format(load_temp$Datetime, "%u"))   # 1 = Monday, 7 = Sunday
mo <- as.integer(format(load_temp$Datetime, "%m"))
hr <- as.integer(format(load_temp$Datetime, "%H"))

keep <- wd %in% 2:4 &                       # Tue, Wed, Thu
        mo %in% c(2:3, 5:7, 9:11) &         # no vacation months
        hr %in% 10:15 &                     # central hours
        load_temp$Load > 15000 &            # drop low-demand outliers
        load_temp$Temperature < 35          # drop extreme temperatures

demand_temp <- data.frame(
  Temperature = load_temp$Temperature[keep],
  Load        = load_temp$Load[keep]
)
row.names(demand_temp) <- NULL

stopifnot(nrow(demand_temp) > 0, !anyNA(demand_temp))

save(demand_temp, file = "data/demand_temp.rda",
     version = 3, compress = "xz")
