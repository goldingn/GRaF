n <- 100
m <- 3
x <- matrix(rnorm(m * n), ncol = m)
b <- rnorm(m)
f <- x %*% b
y <- rbinom(n, 1, plogis(f))
l <- rep(1, m)
mn <- rep(mean(y), n)
wt <- rep(1, n)
x_df <- data.frame(x)

f <- function(l) {
  m <- graf(y, x_df, l = l)
  return (m$mnll)
}

g <- function(l) {
  m <- graf(y, x_df, l = l)
  return (m$l_grads)
}

o <- optim(par = l, fn = f, gr = g, method = 'BFGS')
o2 <- optim(par = l, fn = f, method = 'BFGS')

o
o2
