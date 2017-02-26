graf <-
  function (y, x, error = NULL, weights = NULL, prior = NULL, l = NULL, opt.l = FALSE,
            theta.prior.pars = c(log(10), 1), hessian = FALSE, opt.control = list(),
            verbose = FALSE, method = c('Laplace', 'EP')) {
    
    method <- match.arg(method)
    
    # optionally optimise graf (by recursively calling this function)
    if (opt.l) {
      
      # get all visible object as a list
      args <- capture.all()
      
      # get the expected objects
      expected_args <- names(formals(graf))
      
      # remove any unexpected arguments
      args <- args[names(args) %in% expected_args]
      
      # pass this to optimiser
      fit <- optimise.graf(args)
      
      # skip out of this function and return
      return (fit)
      
    }
    
    
    if (!is.data.frame(x)) stop ("x must be a dataframe")
    
    # convert any ints to numerics
    for(i in 1:ncol(x)) if (is.integer(x[, i])) x[, i] <- as.numeric(x[, i])
    
    obsx <- x
    k <- ncol(x)
    n <- length(y)
    
    if (is.null(weights)) {
      # if weights aren't provided
      weights <- rep(1, n)
    } else {
      # if they are, run some checks
      # throw an error if weights are specified with EP
      if (method == 'EP') {
        stop ('weights are not implemented for the EP algorithm (yet)')
      }
      # or if any are negative
      if (any(weights < 0)) {
        stop ('weights must be positive or zero')
      }
    }
    
    # find factors and convert them to numerics
    notfacs <- 1:k
    facs <- which(unlist(lapply(x, is.factor)))
    if (length(facs) > 0) notfacs <- notfacs[-facs]
    for (fac in facs) {
      x[, fac] <- as.numeric(x[, fac])
    }
    x <- as.matrix(x)
    
    # scale the matrix, retaining scaling
    scaling <- apply(as.matrix(x[, notfacs]), 2, function(x) c(mean(x), sd(x)))
    for (i in 1:length(notfacs)) {
      x[, notfacs[i]] <- (x[, notfacs[i]] - scaling[1, i]) / scaling[2, i]
    }
    
    # set up the default prior, if not specified
    exp.prev <- sum(weights[y == 1]) / sum(weights)
    if (is.null(prior))  mnfun <- function(x) rep(exp.prev, nrow(x))
    else mnfun <- prior
    
    # give an approximation to l, if not specified (or optimised)
    if (is.null(l)) {
      l <- rep(0.01, k)
      l[notfacs] <- apply(x[y == 1, notfacs, drop = FALSE], 2, sd) * 8
      l[l == 0] <- 1
    }
    
    # calculate mean (on unscaled data and probability scale)
    mn <- mnfun(obsx)
    
    # fit model
    if (method == 'Laplace') {
      # by Laplace approximation
      fit <- graf.fit.laplace(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = verbose)
    } else {
      # or using the expectation-propagation algorithm
      fit <- graf.fit.ep(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = FALSE)
    }
    
    fit$mnfun = mnfun
    fit$obsx <- obsx
    fit$facs <- facs
    fit$hessian <- hessian
    fit$scaling <- scaling
    fit$peak = obsx[which(fit$MAP == max(fit$MAP))[1], ]
    class(fit) <- "graf"
    fit
  }

capture.all <- function() {
  # capture all visible objects in the parent environment and pass to a list
  env <- parent.frame()
  object_names <- objects(env)
  objects <- lapply(object_names,
                    get,
                    envir = env)
  names(objects) <- object_names
  return (objects)
}