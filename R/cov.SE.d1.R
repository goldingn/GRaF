cov.SE.d1 <- function (x, e = NULL, l) {
  # get gradients (matrices) of the kernel wrt. the parameters
  # CURRENTLY IGNORES e!!
  
  # number of parameters
  n <- length(l)
  
  # assign vector for gradients
  grads <- list()
  
  # get full covariance matrix
  K <- cov.SE(x1 = x, e1 = e, l = l)
  
  # loop through them
  for (i in 1:n) {
    
    # squared distances
    d2_i <- as.matrix(dist(x[, i]) ^ 2)
    
    # gradient for each parameter
    grads[[i]] <- K * (1 / l[i] ^ 2) * d2_i / 2
    
  }
  
  # return as a list
  return (grads)
  
}