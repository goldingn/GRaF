\name{plot3d}
\alias{plot3d}

\title{
Surface plots of \code{graf} models
}
\description{
\code{plot3d} provides basic surface plots to visualise interactions betweens fitted terms in \code{graf} models. The predicted probability of presence along two covariates is plotted with other covariates are held either at their average value or near the peak (point of highest probability) of the latent random field. The prior function and credible intervals may also be plotted. 3d plots for factors are not implemented.\cr

}
\usage{
plot3d(object, dims = c(1, 2), resolution = 20, CI = 0.95, prior = TRUE, CIs = TRUE,
		theta = 30, phi = 30, peak = FALSE, ...)
}

\arguments{
  \item{object}{
A \code{graf} object.
}
  \item{dims}{
A vector of length 2 giving the column numbers of the covariates to plot.
}
  \item{resolution}{
The resolution of the prediction surface along each covariate.
}
  \item{CI}{
The level at which to calculate predictive credible intervals. The default value returns upper and lower 95\% credible intervals.
}
  \item{prior}{
Whether to additionally plot the value of the prior function.
}
  \item{CIs}{
Whether to additionally plot the credible intervals.
}
  \item{theta}{
Viewing angle for \code{persp}, theta gives the azimuthal direction.
}
  \item{phi}{
Viewing angle for \code{persp}, phi gives the colatitude.
}
  \item{peak}{
Whether to hold other covariates at the value of the training datapoint with the highest probability of presence. If \code{FALSE} covariates are held at their mean value in the training data set.
}
  \item{\dots}{
Additional arguments to pass to \code{persp}.
}
}

\value{
One or more surface plots of the marginal effect of a pair of predictors on the probability of species presence.
}

\seealso{
\code{\link{graf}}, \code{\link{plot.graf}}
}