
# density and gradient of the prior over the log hyperparameters
theta.prior <- function(theta, pars) {
  
  if (any(is.na(pars))) {
    return(list(density = 0,
                gradient = rep(0, length(theta))))
  }
  
  density <- sum(dnorm(theta, pars[1], pars[2], log = TRUE))
  
  gradient <- (pars[1] - theta) / pars[2] ^ 2
  
  return(list(density = density,
              gradient = gradient))
  
}

# define the objective and gradient functions to optimise the hyperparameters
# of a GRaF model
objective <- function(theta, prior.pars, isfac, args, fun) {
  
  # unpack theta
  l <- ifelse(isfac, 0.01, NA)
  l[!isfac] <- exp(theta)
  
  # set the lengthscales
  args$l <- l
  
  # run the model
  m <- do.call(fun, args)
  
  # log likelihood and prior
  llik <- -m$mnll
  lpri <- theta.prior(theta, prior.pars)$density
  
  # log posterior
  lpost <- llik + lpri
  
  # and objective
  objective <- -lpost
  
  return (objective)
  
}


gradient <- function(theta, prior.pars, isfac, args, fun) {
  
  # unpack theta
  l <- ifelse(isfac, 0.01, NA)
  l[!isfac] <- exp(theta)
  
  # set the lengthscales
  args$l <- l
  
  # run the model
  m <- do.call(fun, args)
  
  # gradient of llik w.r.t. l
  dLdl <- m$l_grads[!isfac]
  
  # gradient of l w.r.t. theta
  dldtheta <- exp(theta)
  
  # gradient of lpri w.r.t. theta
  dpdtheta <- theta.prior(theta, prior.pars)$gradient
  
  # gradient of llik w.r.t. theta
  dLdtheta <- dLdl * dldtheta
  
  # gradient of lpost w.r.t. theta
  dPdtheta <- dLdtheta + dpdtheta
  
  # gradient of objective w.r.t. lpost
  dOdP <- -1
  
  # gradient of objective w.r.t. theta
  dOdtheta <- dOdP * dPdtheta 
  
  return (dOdtheta)
  
}

optimise.graf <- function(args) {
  # pass all the arguments of a call to graf, memoize and
  # optimise the model, and return afitted version
  
  # set opt.l to FALSE
  args[['opt.l']] <- FALSE
  
  # memoise graf
#   mgraf <- memoise(graf)
  mgraf <- graf
  
  # set up initial lengthscales
  k <- ncol(args$x)
  l <- rep(1, k)
  
  # find factors and drop them from theta
  notfacs <- 1:k
  
  facs <- which(unlist(lapply(args$x, is.factor)))
  
  if (length(facs) > 0) {
    notfacs <- notfacs[-facs]
    l[facs] <- 0.01
  }
  
  # log them
  theta <- log(l[notfacs])

  # logical vector of factors
  isfac <- 1:k %in% facs
  
  # optimisation arguments
  if (args$method == 'Laplace') {
    meth <- 'L-BFGS-B'
    grad <- gradient
  } else {
    meth <- 'BFGS'
    grad <- NULL
  }

  low <- -Inf 
  up <- Inf
  
  opt <- optim(theta,
               fn = objective,
               gr = grad,
               prior.pars = args$theta.prior.pars,
               isfac = isfac,
               args = args,
               fun = mgraf,
               hessian = args$hessian,
               lower = low,
               upper = up,
               method = meth,
               control = args$opt.control)
  
  # get the resultant lengthscales
  l[notfacs] <- exp(opt$par)
  
  args$l <- l
  
  # fit the final model and return
  m <- do.call(mgraf, args)

  # replace hessian with the hessian matrix or NULL
  if (args$hessian) {
    m$hessian <- opt$hessian
  } else {
    m$hessian <- NULL
  }
  
  # un-memoize graf
  forget(mgraf)
  
  return (m)
  
}