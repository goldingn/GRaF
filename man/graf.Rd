\name{graf}
\alias{graf}

\title{
Fit a GRaF model to binary data
}
\description{
Sets up and fits a GRaF model, returning an object of class \code{graf} for which functions are available to visualise, compare and make predictions from models.
}
\usage{
graf(y, x, error = NULL, weights = NULL, prior = NULL, l = NULL, opt.l = FALSE,
		theta.prior.pars = c(log(10), 1), hessian = FALSE, opt.control = list(),
		verbose = FALSE, method = c('Laplace', 'EP'))
}
\details{
In addition to the dataframe of covariates \code{x} a matrix of corresponding standard deviations, describing the measurement error, may be specified via the \code{error} argument.\cr\cr
If prior knowledge about the species' ecology is available (e.g. environmental thresholds beyond which the species cannot persist), this may be incorporated into the GRaF model by providing an R function via the \code{prior} argument.\cr\cr
The smoothness of the fitted random field is controlled by a vector of lengthscales which is used to parameterise a squared exponential covariance function, with one lengthscale for each covariate. \code{graf} centres and standardises covariates before fitting and the lengthscales therefore relate to these scaled covariates.
Lengthscales can be specified by the user via the the \code{l} argument, roughly approximated by \code{graf} by setting \code{l = NULL} or optimised by setting \code{opt.l = TRUE}.
}
\arguments{
  \item{y}{
A vector denoting presence (1) and absence (0) records.
}
  \item{x}{
A dataframe of covariates which may contain factors as well as continuous variables. Factor covariates can still interact with other covariates, but they are not infuenced by either the mean function or the length-scale parameters and therefore give a maximum-likelihood estimate.
}
  \item{error}{
An optional matrix of standard deviations associated with x. If this is missing, covariates are assumed to be measured without error.
}
  \item{weights}{
An optional vector of weights to be used in the fitting process. Weights must be positive or zero. Observations with a weight of zero will be ignored.
}
  \item{prior}{
An optional R function providing an \emph{a priori} estimate of the probability of presence of the species given the covariates. The function must take a dataframe of environmental covariates (matching \code{x}) as input and return a corresponding vector of the probability of presence.
If \code{NULL} a flat prior is used which gives the species' prevalence as the probability of presence at all sites. 
}
  \item{l}{
An optional vector of lengthscales for each dimension of the Gaussian field. The lengthscale controls the smoothness of the fitted function, with a higher value giving a flatter curve
If \code{NULL} lengthscales are approximated as 8 times the ratio between the standard deviation of environmental variables at presence sites and the standard deviation of environmental variables at all sites. If \code{opt.l = TRUE} and \code{l} is specified it gives the starting position for the opimisation algorithm. Lengthscales specified for factors will be ignored.
}
  \item{opt.l}{
Whether to optimise the lengthscale parameters with respect to the model marginal likelihood using \code{\link{optim}} (with method \code{BFGS}) and return a model fitted with the maximum \emph{a posteriori} estimate of the lengthscales. With a large number of data points this can be quite time consuming.
}
  \item{theta.prior.pars}{
If \code{opt.l = TRUE}, allows the user to specify the mean and standard deviation of the normal hyperprior over the log-lengthscale hyperparameter theta. By default the mean is set at log(10) and the standard deviation at 1. 
}
  \item{hessian}{
If \code{opt.l = TRUE}, whether to calculate the hessian matrix whilst optimising lengthscales. Calculating the hessian is more time consuming and is not recommended for general use, but advanced users could use it to approximate the full posterior density of the lengthscale hyperparameters.
}
  \item{opt.control}{
If \code{opt.l = TRUE}, allows the user to pass a list of options directly to the \code{control} argument of the \code{optim} function used to  carry out the lengthscale optimisation.
For large models, adjusting control parameters such as the convergence tolerance (\code{reltol}) can reduce reduce computation time, with an associated loss of precision in the optimisation. See \code{control} in \code{\link{optim}} for more details.
}
  \item{verbose}{
Logical, whether to display (limited) progress information for the fitting algorithms.
}
  \item{method}{
The approximation method used to fit model. Currently only \code{'Laplace'} (Laplace approximation - the default) and \code{'EP'} (the expectation-propagation algorithm) are available. Whilst the EP algorithm is known to provide a more accurate approximation than Laplace under certain circumstances (see e.g. section 3.7.2 of Rasmussen & Williams (2006)), it has the disadvantages of much greater computational cost and no method for incorporating regression weights (at least in this implementation). In practice, models fitted using the two methods are likely to be very similar. The two approximations are described in detail by Rasmussen & Williams (2006).
}
}

\value{
A \code{graf} object which can be viewed and manipulated using \code{print}, \code{plot} and \code{predict} functions.
}
\references{
Golding & Purse (2013) GRaF: Fast and Flexible Bayesian Species Distribution Modelling Using Gaussian
Random Fields. Manuscript in preparation, preprint available here: \url{http://dx.doi.org/10.6084/m9.figshare.767289}

Rasmussen & Williams(2006) Gaussian Processes for Machine Learning. \url{http://www.gaussianprocess.org/gpml/chapters/}
}

\seealso{
\code{\link{print.graf}}, \code{\link{plot.graf}}, \code{\link{predict.graf}}, \code{\link{plot3d}}, \code{\link{DIC}}, \code{\link{optim}}.
}
\examples{
# load Anguilla data from the dismo package
data(Anguilla_train)

# get the first 100 presence-absence records
y <- Anguilla_train$Angaus[1:100]
# get covariates, removing LocSed (contains NAs) and make DSDam a factor
x <- Anguilla_train[1:100, 3:13]
x$DSDam <- factor(x$DSDam)

# fit a presence-absence model to the data
m1 <- graf(y, x)

#print a brief summary of the model
m1  # print(m1)

# plot the (conditional) effect of each term
par(mfrow = c(3, 4))
plot(m1)

# visualise the interaction between SegTSeas and SegLowFlow
par(mfrow = c(2, 2), mar = rep(2, 4))
plot3d(m1, dims = c(2, 3), theta = -110)

# fit another model with flatter responses
m2 <- graf(y, x, l = m1$ls * 2)

# compare the two models by DIC
DIC(m1)
DIC(m2)

# m2 has lower DIC so could be considered preferable

# predict back to dataset
pred <- predict(m2, x)
head(pred)
}