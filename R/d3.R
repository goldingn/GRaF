d3 <- function(z, y) {
  pr <- y > 0 & y < 1
  npr <- !pr
  ans <- vector('numeric', length(y))
  y[npr] <- 2 * y[npr] - 1
  ans[pr] <- 0
  n <- dnorm(z[npr])
  p <- pnorm(y[npr] * z[npr])
  f <- z[npr]
  a <- n / p ^ 3
  b <- f * (y[npr] ^ 2 + 2) * n * p
  c <- y[npr] * (f ^ 2 - 1) * p ^ 2 + 2 * y[npr] * n ^ 2
  ans[npr] <- a * (b + c)
  ans
}