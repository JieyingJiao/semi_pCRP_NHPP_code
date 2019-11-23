library(spatstat)


##' data generation
##' lamnbda0.mat: matrix of baseline intensity
##' beta is a two dimension vector for (x, y) coefficients
##' grid = c(nrow, ncol)
##' win1 (for x), win2 are two dimension vectors showing area to generate data
data_simu <- function(lambda0.mat, beta, grid, win1, win2, proc.no) {
  p <- length(beta)
  set.seed(1234)
  X_grid <- rnorm(prod(grid)*p, mean = 0, sd = 1)
  X_grid <- matrix(X_grid, nrow = prod(grid), ncol = p)
  
  lambda.mat <- as.vector(t(lambda0.mat)) * exp(X_grid %*% beta)
  
  lambda.mat <- matrix(lambda.mat, nrow = grid[1], byrow = TRUE)
  
  set.seed(proc.no)
  pp <- rpoispp(as.im(lambda.mat, W = owin(win1, win2)))
  
  
  
  return(list(x = pp$x, y = pp$y, X_grid = matrix(X_grid, ncol = p)))
}