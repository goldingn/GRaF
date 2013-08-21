\name{predict.graf}
\alias{predict.graf}

\title{
Predictions for \code{graf} objects.
}
\description{
Predictions and associated credible intervals for \code{graf} objects, either for the data used for fitting or for a different dataset.
}
\usage{
\method{predict}{graf}(object, newdata = NULL, type = c("response", "latent"), CI = 0.95,
			maxn = NULL, ...)
}

\arguments{
  \item{object}{
A \code{graf} object.
}
  \item{newdata}{
An optional dataframe giving covariates to predict to. If \code{NULL} a prediction for the covariates used to fit the model is returned.
}
  \item{type}{
The level of the prediction. \code{"response"} (default) gives predictions on the probability scale and \code{"latent"} gives predictions on the scale of the latent Gaussian.
}
  \item{CI}{
The level at which to calculate predictive credible intervals.
The default value returns upper and lower 95\% credible intervals.
If credible intervals are not required this can be set to \code{NULL}.
If \code{type = 'latent'} setting \code{CI = 'std'} returns the mean and standard deviation of the latent field.
}
  \item{maxn}{
The maximum number of records to predict to in each batch. To avoid computationally expensive operations on large matrices, \code{predict.graf} splits the dataset for prediction into batches with maximum number of records \code{maxn}. This can be adjusted by the user to optimise computational efficiency on different machines. If \code{maxn = NULL} maxn is set at approximately one tenth the number rows used to fit the model.
}
  \item{\dots}{
Additional arguments for future versions.
}
}

\value{
A matrix of posterior modes and optionally credible intervals.
}

\seealso{
\code{\link{graf}}
}