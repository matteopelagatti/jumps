% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpfj.R
\name{hpfj}
\alias{hpfj}
\title{HP filter with automatic jumps detection.}
\usage{
hpfj(y, maxsum = sd(y), edf = TRUE, parinit = NULL)
}
\arguments{
\item{y}{vector with the time series;}

\item{maxsum}{maximum sum of additional level standard deviations;}

\item{edf}{boolean if TRUE computes effective degrees of freedom otherwise computes
the number of degrees of freedom in the LASSO-regression way;}

\item{parinit}{either NULL or vector of 3+n parameters with starting values for the
optimizer; the order of the parameters is sd(slope disturbnce), sd(observatio noise),
square root of gamma, n additional std deviations for the slope.}
}
\value{
list with the following slots:
\itemize{
 \item opt: the output of the optimization function (nloptr)
 \item nobs: number of observations
 \item df: number of estimated parameters (model's degrees of freedom)
 \item loglik: value of the log-likelihood at maximum
 \item ic: vector of information criteria (aic, aicc, bic, hq)
 \item smoothed_level: vector with smoothed level with jumps (hp filter with jumps)
 \item var_smoothed_level: variance of the smoothed level
}
}
\description{
This is the lower-level function for the HP filter with jumps.
The user should use the \code{hpj} function instead, unless in need of more control
and speed. The function estimates the HP filter with jumps.
Jumps happen contextually in the level and in the slope: the standard deviation
of the slope disturbance is \eqn{\gamma} times the standard deviation of the
level disturbance at time \eqn{t}.
The HP smoothing parameter \eqn{\lambda} is estimated via MLE (assuming normally
distributed disturbances) as in Wahba (1978):
\eqn{\lambda = \sigma^2_\varepsilon / \sigma^2_\zeta}.
}
\examples{
set.seed(202311)
n <- 100
mu <- 100*cos(3*pi/n*(1:n)) - ((1:n) > 50)*n - c(rep(0, 50), 1:50)*10
y <- mu + rnorm(n, sd = 20)
plot(y, type = "l")
lines(mu, col = "blue")
hp <- hpfj(y, 60)
lines(hp$smoothed_level, col = "red")
}
\references{
Whaba (1978)
"Improper priors, spline smoothing and the problem of guarding against model errors in regression",
*Journal of the Royal Statistical Society. Series B*, Vol. 40(3), pp. 364-372.
DOI:10.1111/j.2517-6161.1978.tb01050.x
}
