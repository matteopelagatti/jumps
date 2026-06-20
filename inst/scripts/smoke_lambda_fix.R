# Smoke test for the two lambda-fix changes in ssj_mle.R:
#   1. lambda_init_candidates now includes small values {1e-3, 1e-1, 1, 10, base}
#   2. auto_ssj_mle breaks the warm-start when obs_var = lambda*sigma^2 < 1e-4*var(y)
#
# Expected: for g1 at sigma=0.15, most reps should now find lambda in [0.1, 100]
# and detect ~3 discontinuities (was: lambda -> 0 in ~95% of reps, 0 disc detected).

library(jumps)

set.seed(123)
n_rep <- 20L
N     <- 100L
sigma <- 0.05

g1 <- function(x) besselJ(20 * x, nu = 1L) + x * (x >= 0.3 & x <= 0.4) - x * (x >= 0.6)
true_disc <- c(0.3, 0.4, 0.6)

active_gamma <- function(gamma, x) {
  h <- mean(diff(sort(x)))
  gamma > pmax(max(gamma) * 0.05, 0.30 * h)
}

lam_vals  <- numeric(n_rep)
n_disc    <- integer(n_rep)

cat("rep | lambda_mle       | n_disc\n")
cat("----+------------------+-------\n")
for (r in seq_len(n_rep)) {
  x   <- sort(runif(N))
  y   <- g1(x) + rnorm(N, sd = sigma)
  msg <- seq(0, sd(y) / mean(diff(x)), length.out = 21L)
  out <- tryCatch(auto_ssj_mle(y, x, grid = msg, ic = "bic"), error = function(e) NULL)
  if (!is.null(out)) {
    plot(x, g1(x), type = "l")
    lines(x, out$smoothed_level, col = "blue")
    abline(v = x[-1][active_gamma(out$gamma, x)])
    lam_vals[r] <- out$pars["lambda"]
    midpts      <- (x[-N] + x[-1L]) / 2L
    n_disc[r]   <- sum(active_gamma(out$gamma, x))
  }
  cat(sprintf("%3d | %16.4g | %d\n", r, lam_vals[r], n_disc[r]))
}

cat("\nSummary:\n")
cat("  lambda: median =", signif(median(lam_vals), 3),
    "  range = [", signif(min(lam_vals), 3), ",", signif(max(lam_vals), 3), "]\n")
cat("  n_disc distribution:\n")
print(table(n_disc))
cat("  Reps with exactly 3 disc:", sum(n_disc == 3L), "/", n_rep, "\n")
