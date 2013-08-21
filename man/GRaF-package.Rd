\name{GRaF-package}
\alias{GRaF-package}
\alias{GRaF}
\docType{package}
\title{
Species distribution modelling using latent Gaussian random fields
}
\description{
GRaF fits Gaussian random field species distribution models using Laplace approximations.
}

\details{
\tabular{ll}{
Package: \tab GRaF\cr
Type: \tab Package\cr
Version: \tab 0.1-11\cr
Date: \tab 2013-06-12\cr
License: \tab GPL (>= 2)\cr
}
GRaF provides functions to fit, visualise, and predict from Gaussian random field models for binary data, with an emphasis on species distribution models. The principle function is \code{\link{graf}} which returns a \code{graf} object, which can be manipulated using \code{plot}, \code{print} and \code{predict} methods among others.

GRaF is under active development and future releases will enable a range of additional methods, which may (or may not) include: models for count data, spatially and temporally explicit models, the ability to limit interactions between covariates, multi-species models and input and output of raster imagery. Please feel free to contact the maintainer with bug reports and suggestions!
}
\author{
Nick Golding <nick.golding@linacre.ox.ac.uk>
}
\references{
Golding & Purse (2013) GRaF: Fast and Flexible Bayesian Species Distribution Modelling Using Gaussian
Random Fields. Manuscript in preparation.
}
\keyword{ package }
\seealso{
\code{\link{graf}}
}

