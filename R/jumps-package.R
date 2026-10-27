#' jumps: Hodrick-Prescott Filter and Spline Smoothing with Jumps
#'
#' Extraction of a smooth trend that is allowed to have a small number of
#' automatically selected discontinuities (jumps, or breaks).
#'
#' The package provides three families of functions.
#' \itemize{
#'   \item \strong{Hodrick-Prescott filter with jumps.} The jumps are selected
#'     by a LASSO-like regularization, and the optimal regularization constant
#'     is chosen on a grid using information criteria. The main entry point is
#'     \code{\link{hpj}}; \code{\link{hpfj}}, \code{\link{hpfj_fix}} and
#'     \code{\link{hpfjx}} (with regressors) are the lower-level engines.
#'   \item \strong{Smoothing splines with jumps by time warping.} The physical
#'     spacing \eqn{\delta_i} between consecutive sampling sites is replaced by
#'     a latent spacing \eqn{\tau_i = \delta_i + \gamma_i}, and the warping
#'     variables \eqn{\gamma_i \ge 0} are estimated under a budget constraint.
#'     The method is that of Pelagatti, Maranzano and Gatto Ronchieri (2026).
#'     See \code{\link{ssj_mle}}, \code{\link{ssj_kfs}}, \code{\link{ssj_atw}}
#'     and \code{\link{ssj_fix}}, together with their \code{auto_*} wrappers.
#'   \item \strong{Cubic smoothing splines with discontinuities (CSSD).} An
#'     implementation of the dynamic-programming algorithm of Storath and
#'     Weinmann (2024): \code{\link{ssj_cssd_dp}},
#'     \code{\link{auto_ssj_cssd}} and \code{\link{predict_cssd_dp}}, plus
#'     the faster piecewise-constant approximation \code{\link{ssj_cssd}}.
#' }
#'
#' Deterministic seasonal regressors can be built with \code{\link{trigseas}}
#' and \code{\link{dummyseas}}.
#'
#' @references
#' Maranzano, P. and Pelagatti, M. (2025). A Hodrick-Prescott filter with
#' automatically selected breaks. \emph{Economic Modelling}, \strong{150},
#' 107132. \doi{10.1016/j.econmod.2025.107132}
#'
#' Pelagatti, M., Maranzano, P. and Gatto Ronchieri, G. C. (2026). Smoothing
#' splines with some discontinuities. \emph{arXiv preprint} arXiv:2610.33137.
#' \doi{10.48550/arXiv.2610.33137}
#'
#' Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
#' signals. \emph{Journal of Computational and Graphical Statistics},
#' \strong{33}(2), 651--664. \doi{10.1080/10618600.2023.2262000}
#'
#' @keywords internal
#' @useDynLib jumps, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom nloptr nloptr
#' @importFrom stats BIC logLik nobs optim qnorm sd time var
#' @importFrom graphics abline layout lines mtext par points polygon
#' @importFrom grDevices adjustcolor
"_PACKAGE"

# The ggplot2 branch of plot.ssj() refers to columns through the .data pronoun,
# which the tidy-evaluation data mask supplies at draw time.  R CMD check has no
# way of knowing that, so declare it here rather than take a NOTE.  Declaring
# the pronoun is enough: the column names themselves are reached through it.
utils::globalVariables(".data")
