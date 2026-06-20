library(jumps)
set.seed(42)
N <- 100L
x <- sort(runif(N, 0, 5))   # non-unit range to stress the scaling
g <- function(x) sin(4*x) + ifelse(x >= 2.5, 1, 0)
y <- g(x) + rnorm(N, sd = 0.2)

# ssj_kfs
k1 <- ssj_kfs(y, x, lambda = 1, maxsum = 2)
stopifnot(isTRUE(all.equal(k1$x, sort(x))))
stopifnot(abs(mean(k1$smoothed_level) - mean(g(x))) < 0.3)
cat("ssj_kfs OK: range(f) =", round(range(k1$smoothed_level), 3), "\n")

# auto_ssj_kfs
ak <- auto_ssj_kfs(y, x, lambda = 1)
cat("auto_ssj_kfs OK: n_disc =",
    sum(ak$gamma > 0.05 * diff(range(x)) / (N - 1L)), "\n")

# ssj_atw
a1 <- ssj_atw(y, x, lambda = 1, maxsum = 2)
stopifnot(isTRUE(all.equal(a1$x, sort(x))))
stopifnot(abs(mean(a1$smoothed_level) - mean(g(x))) < 0.3)
cat("ssj_atw OK: range(f) =", round(range(a1$smoothed_level), 3), "\n")

# auto_ssj_atw
aa <- auto_ssj_atw(y, x, lambda = 1)
cat("auto_ssj_atw OK: n_disc =",
    sum(aa$gamma > 0.05 * diff(range(x)) / (N - 1L)), "\n")

# auto_ssj_atw_gcv
ag <- auto_ssj_atw_gcv(y, x)
cat("auto_ssj_atw_gcv OK: lambda =", round(ag$lambda, 4),
    ", n_disc =", sum(ag$gamma > 0.05 * diff(range(x)) / (N - 1L)), "\n")

cat("All checks passed.\n")
