rm(list=ls())
library(MHadaptive)
library(purrr)
library(reshape)
library(ggplot2)
library(matrixStats)
library(coda)

source("MCMC_Gibbs.R")
source("data_generation.R")



##' calculate Random index, given true index vector and full posterior samples
##' samples: full sample from pp_Gibbs function (including other pars)
RI <- function(z_true, samples) {
  n <- length(z_true)
  
  z0_pair <- combn(z_true, m = 2)
  
  ri <- rep(0, length(samples))
  for (i in 1:length(ri)) {
    z_pair <- combn(samples[[i]]$index_cluster, m = 2)
    a <- sum((z0_pair[1, ] == z0_pair[2, ]) * (z_pair[1, ] == z_pair[2, ]))
    b <- sum((z0_pair[1, ] != z0_pair[2, ]) * (z_pair[1, ] != z_pair[2, ]))
    ri[i] <- (a+b) / choose(n, 2)
  }
  return(ri)
}


##' given posterior samples after burn-in and thinning
dahl <- function(samples) {
  n <- length(samples[[1]]$index_cluster)
  niter <- length(samples)
  z_sample <- matrix(0, nrow = niter, ncol = n)
  for (i in 1:niter) {
    z_sample[i, ] <- samples[[i]]$index_cluster
  }
  BList <- map(1:niter, ~outer(z_sample[.x,], z_sample[.x,], "=="))
  BBar <- Reduce("+", BList) / niter
  SSE <- map_dbl(BList, ~sum((.x - BBar)^2))
  return(list(min.sse = min(SSE), cLS = which.min(SSE), 
              cluster = as.numeric(z_sample[which.min(SSE),])))
}

## calculate MSE over the grid
MSE_count <- function(simudata, samples, Dahl_index, grid = c(20, 20), A = 1) {
  n <- prod(grid)
  X_grid <- simudata$X_grid
  colindex <- ceiling(simudata$x)
  rowindex <- ceiling(simudata$y)
  index_grid <- (rowindex - 1) * grid[2] + colindex
  count_grid0 <- as.matrix(table(index_grid)) # may have missing grid
  count_grid <- rep(0, n)
  count_grid[as.numeric(rownames(count_grid0))] <- count_grid0
  

  index_cluster <- samples[[Dahl_index]]$index_cluster
  lambda0 <- samples[[Dahl_index]]$lambda0[index_cluster]
  
  beta_sample <- matrix(0, nrow = length(samples), ncol = length(samples[[1]]$beta))
  for (i in 1:length(samples)) {
    beta_sample[i, ] <- samples[[i]]$beta
  }
  betafit <- colMeans(beta_sample)
  
  lambda <- lambda0 * drop(exp(X_grid %*% betafit))
  
  return(mean((count_grid - A * lambda)^2))
}

## calculate modified BIC, given samples after burn-in and thinning
BIC <- function(samples, Dahl_index, X_grid, index_grid, A = 1) {
  N <- length(index_grid)
  p <- length(samples[[1]]$beta)
  
  beta_sample <- matrix(0, nrow = length(samples), ncol = p)
  for (i in 1:length(samples)) {
    beta_sample[i, ] <- samples[[i]]$beta
  }
  beta_fitted <- colMeans(beta_sample)
  
  lambda_data <- samples[[Dahl_index]]$lambda0[samples[[Dahl_index]]$index_cluster[index_grid]]
  lambda_data <- lambda_data * drop(exp(X_grid[index_grid, ] %*% beta_fitted))
  lambda_area <- samples[[Dahl_index]]$lambda0[samples[[Dahl_index]]$index_cluster]
  lambda_area <- lambda_area * drop(exp(X_grid %*% beta_fitted))
  
  loglikelihood <- sum(log(lambda_data)) - A * sum(lambda_area)
  
  n_par <- length(samples[[Dahl_index]]$lambda0)
  
  return(-2 * loglikelihood + n_par * log(N))
}



beta <- c(0.5, 0.5, 0, 0)
## generate data on 20 by 20 grid
grid <- c(20, 20)
win1 <- c(0, 20)
win2 <- c(0, 20)
A <- 1
proc.no <- 1

lambda0.mat <- matrix(0.2, nrow = grid[1], ncol = grid[2])
lambda0.mat[1:5, 1:5] <- 5
lambda0.mat[14:18, 2:6] <- 5
lambda0.mat[8:12, 12:18] <- 20
z_true <- rep(1, prod(grid))
z_true[c(1:5, 21:25, 41:45, 61:65, 81:85, 262:266, 282:286, 
         302:306, 322:326, 342:346)] <- 2
z_true[c(152:158, 172:178, 192:198, 212:218, 232:238)] <- 3



simudata <- data_simu(lambda0.mat = lambda0.mat, beta = beta, grid = grid, 
                      win1 = win1, win2 = win2, proc.no = proc.no)

n <- prod(grid)
colindex <- ceiling(simudata$x)
rowindex <- ceiling(simudata$y)
index_grid <- (rowindex - 1) * grid[2] + colindex
a <- 1
b <- 1
alpha <- 1
initNCluster <- 5
niterations <- 5000
burnin <- 1000

r0 <- 1

samples <- pp_Gibbs(X_grid = simudata$X_grid, n = n, index_grid = index_grid, 
                    a = a, b = b, alpha = alpha, A = A, r = r0,
                    initNCluster = initNCluster, niterations = niterations)
ri <- RI(z_true, samples$Iterates)
plot(ri, type = "l")
samples_post <- samples$Iterates[(burnin+1):niterations]
Dahl_index <- dahl(samples_post)
## fitted number of clusters
K <- length(samples_post[[Dahl_index$cLS]]$lambda0)

bic <- BIC(samples_post, Dahl_index$cLS, simudata$X_grid, index_grid)

MSE <- MSE_count(simudata, samples_post, Dahl_index$cLS)

## fitted beta
beta_sample <- matrix(0, nrow = niterations, ncol = ncol(simudata$X_grid))
gam_sample <- matrix(0, nrow = niterations, ncol = ncol(simudata$X_grid))
ar_beta <- matrix(0, nrow = niterations, ncol = ncol(simudata$X_grid))
for(i in 1:niterations) {
  beta_sample[i, ] <- samples$Iterates[[i]]$beta
  ar_beta[i, ] <- samples$Iterates[[i]]$ar_beta
  gam_sample[i, ] <- samples$Iterates[[i]]$gam
}

colMeans(ar_beta)

beta_sample <- beta_sample[(burnin+1):niterations, ]
gam_sample <- gam_sample[(burnin+1):niterations, ]
plot(beta_sample[, 1], type = "l") ## trace plot

## significant variables
select_index <- c(1:ncol(simudata$X_grid))[colSums(gam_sample) > colSums(1-gam_sample)]
beta_fitted <- colMeans(beta_sample)[colSums(gam_sample) > colSums(1-gam_sample)]
beta_sd <- colSds(beta_sample)[colSums(gam_sample) > colSums(1-gam_sample)]
HPDinterval(as.mcmc(beta_sample))[colSums(gam_sample) > colSums(1-gam_sample), ]


## fitted intensity heat plot
lambda0 <- samples_post[[Dahl_index$cLS]]$lambda0[samples_post[[Dahl_index$cLS]]$index_cluster]
lambda0_fit <- matrix(lambda0, nrow = grid[1], byrow = TRUE)
colnames(lambda0_fit) <- as.character(seq(1, by = 1, length = 20))
rownames(lambda0_fit) <- as.character(seq(1, by = 1, length = 20))
melted_lambda0fit <- melt(lambda0_fit)
p1 <- ggplot(data = melted_lambda0fit, aes(y=X1, x=X2, fill=value, alpha=0)) + 
  coord_fixed(ratio = 1, ylim = c(0, 20), xlim = c(0, 20))+
  scale_fill_distiller(palette = "Spectral",limits=range(lambda0_fit))+
  geom_tile() + xlab("x") + ylab("y") + scale_alpha(guide = 'none') + 
  labs(title = "Baseline Intensity") +
  theme(legend.position = "bottom", legend.title = element_text(size=10), 
        legend.text = element_text(size = 5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


lambda <- lambda0 * drop(exp(simudata$X_grid %*% colMeans(beta_sample)))
lambda_fit <- matrix(lambda, nrow = grid[1], byrow = TRUE)
colnames(lambda_fit) <- as.character(seq(1, by = 1, length = 20))
rownames(lambda_fit) <- as.character(seq(1, by = 1, length = 20))
melted_lambdafit <- melt(lambda_fit)
p2 <- ggplot(data = melted_lambdafit, aes(y=X1, x=X2, fill=value, alpha=0)) + 
  coord_fixed(ratio = 1, ylim = c(0, 20), xlim = c(0, 20))+
  scale_fill_distiller(palette = "Spectral",limits=range(lambda_fit)) +
  geom_tile() + xlab("x") + ylab("y") + scale_alpha(guide = 'none') + 
  labs(title = "Intensity") +
  theme(legend.position = "bottom", 
        legend.title = element_text(size=10), legend.text = element_text(size = 5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


