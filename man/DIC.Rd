\name{DIC}
\alias{DIC}

\title{
Deviance information criterion for a \code{graf} model
}
\description{
Calculates the deviance information criterion (DIC) and effective parameters for a \code{graf} model.
}
\usage{
DIC(object)
}

\arguments{
  \item{object}{
A \code{graf} object.
}

}
\details{
Information criteria can be used to compare the goodness of fit of two models to the same dataset whilst accounting for model complexity. Commonly used methods such as AIC are difficult to implement for models such a GRaF where the number of parameters in the model is not obvious. Spiegelhalter et al. (2002) proposed the \emph{effective} number of parameters for these models and the associated DIC for model comparison. Here we calculate DIC by numerical integration. If weights are included in the model fitting they are taken into account when calculating the DIC. DIC should therefore not be used to compare between models with different observation weights.\cr
}
\value{
A vector containing:
\item{DIC}{
The deviance information criterion.
}
\item{pD}{
The estimated number of effective parameters in the model.
}
}
\references{
Spiegelhalter, D.J.. Best, N.G., Carlin, B.P., van der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society, Series B}, 64 (4): 583-639.
}

\seealso{
\code{\link{graf}}
}

