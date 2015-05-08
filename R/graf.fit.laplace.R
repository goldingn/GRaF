graf.fit.laplace <-
  function (y, x, mn, l, wt, e = NULL, tol  = 10 ^ -6, itmax = 50,
            verbose = FALSE) {
    
    if (is.vector(x)) x <- as.matrix(x)
    mn <- qnorm(mn)
    n <- length(y)
    
    # create the covariance matrix
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
    
    # an identity matrix for the calculations
    eye <- diag(n) 
    
    # initialise
    a <- rep(0, n)
    f <- mn
    obj.old <- Inf
    obj <- -sum(wt * d0(f, y))
    it <- 0
    
    # start newton iterations
    while (obj.old - obj > tol & it < itmax) {
      it <- it + 1
      obj.old <- obj
      W <- -(wt * d2(f, y))
      rW <- sqrt(W)
      cf <- f - mn
      mat1 <- rW %*% t(rW) * K + eye
      L <- tryCatch(chol(mat1),
                    error = function(x) return(NULL))
      b <- W * cf + wt * d1(f, y)
      mat2 <- rW * (K %*% b)
      adiff <- b - rW * backsolve(L, forwardsolve(t(L), mat2)) - a 
      dim(adiff) <- NULL
      
      # find optimum step size using Brent's method
      res <- optimise(psiline, c(0, 2), adiff, a, as.matrix(K), y, d0, mn, wt)
      a <- a + res$minimum * adiff
      f <- K %*% a + mn
      obj <- psi(a, f, mn, y, d0, wt)
      
    }
    
    # recompute key components
    cf <- f - mn
    W <- -(wt * d2(f, y))
    rW <- sqrt(W)
    mat1 <- rW %*% t(rW) * K + eye
    L <- tryCatch(chol(mat1),
                  error = function(x) return(NULL))
    
    # return marginal negative log-likelihood
    mnll <- (a %*% cf)[1, 1] / 2 + sum(log(diag(L)) - (wt * d0(f, y)))
    
    # get partial gradients of the objective wrt l
    
    # gradient components
    W12 <- matrix(rep(rW, n), n)
    R <- W12 * backsolve(L, forwardsolve(t(L), diag(rW)))
    C <- forwardsolve(t(L), (W12 * K))
    
    # partial gradients of the kernel
    dK <- cov.SE.d1(x, e, l)
    
    # rate of change of likelihood w.r.t. the mode
    s2 <- (diag(K) - colSums(C ^ 2)) / 2 * d3(f, y)
    
    # vector to store gradients
    l_grads <- rep(NA, length(l))
    
    for (i in 1:length(l)) {
      
      grad <- sum(R * dK[[i]]) / 2
      grad <- grad - (t(a) %*% dK[[i]] %*% a) / 2
      b <- dK[[i]]%*% d1(f, y)
      grad <- grad - t(s2) %*% (b - K %*% (R %*% b))
      l_grads[i] <- grad
      
    }
    
    if(verbose ) cat(paste("  ", it, "Laplace iterations\n"))
    if(it == itmax) print("timed out, don't trust the inference!")
    return(list(y = y, x = x, MAP = f, ls = l, a = a, W = W, L = L, K = K,
                e = e, obsx = x, obsy = y, mnll = mnll, wt = wt,
                l_grads = l_grads))
  }