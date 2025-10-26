\name{jumps-package}
\alias{jumps-package}
\alias{jumps}
\docType{package}
\title{
  Hodrick-Prescott Filter with Jumps
}
\description{
  Extraction of a smooth trend with automatically selected jumps.
}
\details{
  This package implements a version of the HP filter and of smoothing splines that allow discontinuities (jumps). The jumps are automatically selected using a LASSO-like regularization and the optimal regularization constant can be found providing a grid of alternatives and using information criteria. The user can also add regressors such as deterministic seasonal components.
}
\author{
Matteo Pelagatti, Paolo Maranzano.

Maintainer: Matteo Pelagatti <matteo.pelagatti@unimib.it>, Paolo Maranzano <paolo.maranzano@unimib.it> 
}
\references{
  Maranzano and Pelagatti (2024) A Hodrick-Prescott filter with automatically selected jumps.
}
\keyword{ package }