#' Trigonometric seasonal variables
#' 
#' It produces a matrix with seasonal sinusoids to be used as
#' regressors. By default, for \eqn{t=1,\ldots,n} and
#' \eqn{j = 1, \ldots, \lfloor s/2\rfloor}, it computes
#' \deqn{\cos(2\pi j/s), \qquad \sin(2\pi j/s)}.
#' Notice that if $s$ is even the sine function at highest
#' frequency is omitted because it equals zero for all $t$.
#' The used can select a different set of harmonics.
#' 
#' @param n length of time series;
#' @param s seasonal period;
#' @param harmonics vector of harmonics to be used: the cosine and sine
#' functions are computed at frequencies \code{2*pi*harmonics/s}, if there
#' is an element of \code{harmonics} for which \code{2*pi*harmonics/s} is equal
#' to \eqn{\pi} the corresponding sine is omitted.
#' 
#' @returns It returns a matrix with n rows and so many columns
#' as the harmonics (by default these are \eqn{s-1}).
#' 
#' @examples
#' y <- log(AirPassengers)
#' X <- trigseas(length(y), 12)
#' X <- cbind(X, t = 1:length(y))
#' reg <- lm(y~X)
#' @export
trigseas <- function(n, s, harmonics = NULL) {
  if (is.null(harmonics)) harmonics <- 1:floor(s/2)
  fr_cos <- outer(seq(0, n-1), harmonics)*2*pi/s
  fr_sin <- fr_cos[, harmonics != s/2]
  matrix(cbind(cos(fr_cos), sin(fr_sin)), nrow = n,
         dimnames = list(NULL,
                         c(paste0("cos", harmonics),
                           paste0("sin", harmonics[harmonics != s/2])
                         )
         )
  )
}

#' Seasonal dummy variables
#' 
#' It produces a matrix with seasonal dummies summing to zero
#' to be used as regressors. If \eqn{s} is the seasonal period,
#' then the  \eqn{j}-th dummy equals 1 in season \eqn{j} and
#' \eqn{-1/(s-1)} in the other seasons, so that the variable
#' sums to 0 over one year period.
#' 
#' @param n length of time series;
#' @param s seasonal period;
#' @param remove to avoid collinearity with the constant
#' remove a column, by default it is column s; if NULL no
#' column is removed.
#' 
#' @returns It returns a matrix with \eqn{n} rows and \eqn{s-1}
#' column unless \code{remove = NULL}.
#' 
#' @examples
#' y <- log(AirPassengers)
#' X <- dummyseas(length(y), 12)
#' X <- cbind(X, t = 1:length(y))
#' reg <- lm(y~X)
#' @export
dummyseas <- function(n, s, remove = s) {
  vt <- seq(0, n-1)
  vs <- seq(0, s-1)
  L <- outer(vt, vs, function(x, y) floor(x %% s) == y )
  X <- matrix(0, n, s, dimnames = list(NULL, paste0("D", vs+1)))
  X[L] <- 1
  X[!L] <- -1/(s-1)
  if (is.null(remove)) X else X[, -remove]
}

#' Mean squared error
#' 
#' It computes the mean squared difference between the elements
#' of the vectors x and y.
#' 
#' @param y vector
#' @param x vector
#' 
#' @returns The scalar mean squared difference/error.
#' @export
mse <- function(y, x) mean((y-x)^2)
