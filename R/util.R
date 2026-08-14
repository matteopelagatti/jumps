#' Trigonometric seasonal variables
#'
#' Produces a matrix of seasonal sinusoids to be used as regressors. By
#' default, for \eqn{t = 1, \ldots, n} and \eqn{j = 1, \ldots, \lfloor s/2
#' \rfloor}, it computes \deqn{\cos(2\pi j t/s), \qquad \sin(2\pi j t/s).}
#' Notice that if \eqn{s} is even the sine function at the highest frequency
#' is omitted, because it equals zero for all \eqn{t}. The user can select a
#' different set of harmonics.
#'
#' @param n length of the time series.
#' @param s seasonal period.
#' @param harmonics vector of harmonics to be used: the cosine and sine
#'   functions are computed at frequencies \code{2*pi*harmonics/s}; if there is
#'   an element of \code{harmonics} for which \code{2*pi*harmonics/s} equals
#'   \eqn{\pi}, the corresponding sine is omitted.
#'
#' @returns A matrix with \code{n} rows and as many columns as the harmonics
#'   (by default \eqn{s-1}).
#'
#' @seealso \code{\link{dummyseas}} for the seasonal-dummy alternative.
#'
#' @examples
#' y <- log(AirPassengers)
#' X <- trigseas(length(y), 12)
#' X <- cbind(X, t = seq_along(y))
#' reg <- lm(y ~ X)
#'
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
#' Produces a matrix of seasonal dummies summing to zero, to be used as
#' regressors. If \eqn{s} is the seasonal period, then the \eqn{j}-th dummy
#' equals 1 in season \eqn{j} and \eqn{-1/(s-1)} in the other seasons, so that
#' the variable sums to 0 over one year.
#'
#' @param n length of the time series.
#' @param s seasonal period.
#' @param remove index of the column to remove in order to avoid collinearity
#'   with the constant; by default it is column \code{s}. If \code{NULL}, no
#'   column is removed.
#'
#' @returns A matrix with \eqn{n} rows and \eqn{s-1} columns, unless
#'   \code{remove = NULL}, in which case it has \eqn{s} columns.
#'
#' @seealso \code{\link{trigseas}} for the trigonometric alternative.
#'
#' @examples
#' y <- log(AirPassengers)
#' X <- dummyseas(length(y), 12)
#' X <- cbind(X, t = seq_along(y))
#' reg <- lm(y ~ X)
#'
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
#' Computes the mean squared difference between the elements of the vectors
#' \code{y} and \code{x}.
#'
#' @param y numeric vector.
#' @param x numeric vector of the same length as \code{y}.
#'
#' @returns The scalar mean squared difference (error).
#'
#' @examples
#' mse(1:3, 2:4)
#'
#' @export
mse <- function(y, x) mean((y-x)^2)
