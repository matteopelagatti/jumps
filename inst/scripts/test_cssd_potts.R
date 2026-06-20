library(jumps)

# --- Simple step function: jump at x=0.5 ---
cat("=== Step function, jump at 0.5 ===\n")
x_step <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
y_step <- c(0, 0, 0, 1, 1, 1)
for (gam in c(0.1, 0.5, 1, 2)) {
  fit <- ssj_cssd_potts(y_step, x_step, p=0.99, gamma=gam, delta=0.01)
  cat("gamma=", gam, " n_disc=", fit[["n_disc"]])
  if (fit[["n_disc"]] > 0) cat("  locs:", round(fit[["disc_locs"]],3))
  cat("\n")
}

# --- g1 signal, sigma=0.1 (paper's Figure 3 setup) ---
cat("\n=== g1 signal, sigma=0.1, delta=sigma ===\n")
set.seed(2024)
x   <- sort(runif(100))
sig <- 0.10
y   <- besselJ(20*x, 1) + x*(x>=0.3 & x<=0.4) - x*(x>=0.6) + rnorm(100, sd=sig)
for (gam in c(1, 2, 4, 8, 16)) {
  fit <- ssj_cssd_potts(y, x, p=0.999, gamma=gam, delta=sig)
  cat("gamma=", gam, " n_disc=", fit[["n_disc"]])
  if (fit[["n_disc"]] > 0) cat("  locs:", round(fit[["disc_locs"]],3))
  cat("\n")
}

# --- auto_ssj_cssd ---
cat("\n=== auto_ssj_cssd with 5-fold CV ===\n")
fit_auto <- auto_ssj_cssd(y, x, delta=sig)
cat("CV-selected: p=", round(fit_auto[["cv_par"]]["p"], 4),
    " gamma=", round(fit_auto[["cv_par"]]["gamma"], 4), "\n")
cat("n_disc=", fit_auto[["n_disc"]])
if (fit_auto[["n_disc"]] > 0)
  cat("  locs:", round(fit_auto[["disc_locs"]], 3))
cat("\n")
