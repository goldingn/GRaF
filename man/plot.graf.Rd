\name{plot.graf}
\alias{plot.graf}

\title{
Plot \code{graf} models
}
\description{
\code{plot.graf} provides basic plots to visualise fitted terms in \code{graf} models. The predicted probability of presence along each covariate is plotted, along with associate credible intervals, with other covariates are held either at their average value or near the peak (point of highest probability) of the latent random field.\cr
}
\usage{
\method{plot}{graf}(x, vars = NULL, resolution = 50, CI = 0.95, prior = FALSE, data = TRUE, jitter = 1,
		peak = FALSE, ...)
}

\arguments{
  \item{x}{
A \code{graf} object.
}
  \item{vars}{
A vector giving the column index of covariates for which to plot fitted terms. If \code{NULL} fitted terms are plotted for all covariates.
}
  \item{resolution}{
The resolution of the prediction line for each covariate.
}
  \item{CI}{
The level at which to calculate predictive credible intervals. The default value returns upper and lower 95\% credible intervals.
}
  \item{prior}{
Whether to additionally plot the value of the prior function.
}
  \item{data}{
Whether to add rug plots of the training data.
}
  \item{jitter}{
The amount to jitter the observed data points (only applies if \code{data = TRUE}). Set to 0 for no jittering.
}
  \item{peak}{
Whether to hold other covariates at the value of the training datapoint with the highest probability of presence. If \code{FALSE} covariates are held at their mean value in the training data set.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}

\value{
A series of plots (one per covariate) of the marginal effect of covariates (near the peak) on the probability of species presence.
}

\seealso{
\code{\link{graf}}, \code{\link{plot3d}}
}
