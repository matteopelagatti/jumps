% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpj.R
\name{BIC.hpj}
\alias{BIC.hpj}
\title{BIC method for the class hpj}
\usage{
\method{BIC}{hpj}(object, ...)
}
\arguments{
\item{object}{an object of class hpj;}

\item{...}{additional objects of class hpj.}
}
\value{
If just one object is provided,
a numeric value with the corresponding BIC.
If multiple objects are provided, a data.frame with rows corresponding to
the objects and columns representing the number of parameters in the
model (df) and the BIC.
}
\description{
BIC method for the class hpj
}
