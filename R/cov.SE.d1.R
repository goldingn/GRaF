cov.SE.d1 <- function (x, e = NULL, l) {
  # get gradients (matrices) of the kernel wrt. the parameters
  # CURRENTLY IGNORES e!!
  
  # number of parameters
  n <- length(l)
  
  # assign vector for gradients
  grads <- list()
  
  # get full covariance matrix
  K <- cov.SE(x1 = x, e1 = e, l = l)
  
  # transform parameters
  s <- 2 / l ^ 2
  
  # loop through them
  for (i in 1:n) {
    
    # squared distances
    d2_i <- as.matrix(dist(x[, i]) ^ 2)
    
    # gradient for each parameter
    grads[[i]] <- 0.5 * K * s[i] * d2_i
    
  }
  
  # return as a list
  return (grads)
  
}