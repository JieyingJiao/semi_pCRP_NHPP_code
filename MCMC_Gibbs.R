# library(MHadaptive)
##' model: lambda(s) = lambda0(s)exp(X * beta)
##' data: (x, y), length N 
##' grid: n

## variable selection, spike-slab for beta

##' function for Collapsed sampler;
##' input: X_grid (covariates value for every gird, n by p matrix);
##' input: n, number of grid;
##' input: index_grid, index showing which grid the data belongs to, N by 1 vector,
##'        integer valued from 1 to n;
##' input: hyperparameters, a, b, alpha;
##' input: A, area of every grid (equal area now);
##' inpur: r, power of the pCRP;
##' input: initNCluster, initial number of cluster
##' input: niterations, posterior sample size;
##' output: sample, matrix, with column: lambda0 (nCluster), index_cluster (n)
pp_Gibbs <- function(X_grid, n, index_grid, a, b, alpha, A, r, 
                     initNCluster, niterations) {
  N <- length(index_grid)
  count_grid0 <- as.matrix(table(index_grid)) # may have missing grid
  count_grid <- rep(0, n)
  count_grid[as.numeric(rownames(count_grid0))] <- count_grid0
  
  X <- as.matrix(X_grid[index_grid, ])
  ## initial posterior sample
  ## index_cluster: cluster index for each grid, corresponds to "z" in draft
  ## from 1:initNCluster, each cluster has elements (grids), no wasting number for cluster labeling.
  index_cluster <- c(sample(1:initNCluster, size = initNCluster, replace = FALSE),
                     sample(1:initNCluster, size = n-initNCluster, replace = TRUE))
  lambda0 <- rgamma(initNCluster, shape = a, rate = b)
  nCluster <- initNCluster
  beta <- rep(0, ncol(X))
  gam <- rep(0, ncol(X))
  
  History <- vector("list", niterations)
  
  ## start Gibb's sampling
  for (iter in 1:niterations) {
    Lambda <- A*exp(X_grid %*% beta)
    
    ## permutate n grids
    np <- sample(1:n, size = n)
    ## update index_cluster
    for (i in np) {
      count_cluster <- table(index_cluster)
      ## determine if i-th grid is a singleton
      if (count_cluster[index_cluster[i]] > 1) {
        ## if not a singleton, then have nCluster + 1 choice
        count_cluster[index_cluster[i]] <- count_cluster[index_cluster[i]] - 1
        ## goes to an existing cluster: z_i = 1:nCluster
        logclusterProb <- sapply(1:nCluster, function(x) {
          r * log(count_cluster[x]) + count_grid[i] * log(lambda0[x]) - 
            Lambda[i]*lambda0[x]
        })
        ## goes to a new cluster: z_i = nCluster+1
        logclusterProb[nCluster+1] <- log(alpha) + a * log(b) + 
          log(gamma(min(count_grid[i] + a, 150))) - 
          (count_grid[i] + a) * log(b+Lambda[i]) - 
          log(gamma(a))
        ## get the posterior sample for Z_i
        clusterProb <- exp(logclusterProb - max(logclusterProb))
        
        index_i <- sample(1:(nCluster+1), size = 1,
                          prob = clusterProb/sum(clusterProb))
        ## if the i-th grid really goes to a new cluster
        if (index_i > nCluster) {
          lambda0_new <- rep(0, nCluster + 1)
          lambda0_new[1:nCluster] <- lambda0
          lambda0_new[nCluster+1] <- rgamma(1, shape = a, rate = b)
          lambda0 <- lambda0_new
          
          index_cluster[i] <- index_i
          nCluster <- nCluster + 1
        } else { ## if i-th grid goes to an existing cluster
          index_cluster[i] <- index_i
        }
      } else { ## if grid is a singleton, then has nCluster choices
        ## move all the cluster index that are greater than index_cluster[i] 
        ## foward by one to fill in the blank, also change the count_cluster 
        ## and lambda0 coresspondingly
        index_cluster[index_cluster > index_cluster[i]] <- 
          index_cluster[index_cluster > index_cluster[i]] - 1
        lambda0[index_cluster[i] : (nCluster-1)] <- lambda0[(index_cluster[i]+1):nCluster]
        count_cluster <- table(index_cluster)
        count_cluster[index_cluster[i]] <- count_cluster[index_cluster[i]] - 1
        ## only nCluster-1 clusters remained after removing i-th grid
        lambda0 <- lambda0[1:(nCluster-1)]
        ## goes to existing cluster : 1:(nCluster-1)
        logclusterProb <- sapply(1:(nCluster-1), function(x) {
          r * log(count_cluster[x]) + count_grid[i] * log(lambda0[x]) - 
            Lambda[i]*lambda0[x]
        })
        ## goes to new cluster: nCluster
        logclusterProb[nCluster] <- log(alpha) + a * log(b) + 
          log(gamma(min(count_grid[i] + a, 150))) - 
          (count_grid[i] + a) * log(b+Lambda[i]) - 
          log(gamma(a))
        ## get the posterior sample for Z_i
        clusterProb <- exp(logclusterProb - max(logclusterProb))
        index_i <- sample(1:nCluster, size = 1, 
                          prob = clusterProb/sum(clusterProb))
        index_cluster[i] <- index_i
        if (index_i < nCluster) {
          nCluster <- nCluster-1
        } else {
          lambda0_new <- rep(0, nCluster)
          lambda0_new[1:(nCluster-1)] <- lambda0
          lambda0_new[nCluster] <- rgamma(1, shape = a, rate = b)
          lambda0 <- lambda0_new
        }
      }
    }
    
    ## update lambda
    for (i in 1:nCluster) {
      lambda0[i] <- rgamma(1, shape = a+sum(count_grid[index_cluster == i]), 
                           rate = b + sum(Lambda[index_cluster == i]))
    }
    
    
    ## update gam
    gam <- rbinom(length(gam), size = 1, prob = (1+dnorm(beta, sd = 0.1)/dnorm(beta, sd = 10))^(-1))
    
    ## update beta, using MH
    loglikelihood <- function(beta) {
      Lambda <- A*exp(X_grid %*% beta)
      l_grid <- rep(0, n)
      for (i in 1:n) {
        l_grid[i] <- count_grid[i] * log(lambda0[index_cluster[i]]) +  
          sum(as.matrix(X[index_grid == i, ]) %*% beta) - 
          Lambda[i]*lambda0[index_cluster[i]]
      }
      return(sum(l_grid))
    }
    
    ## gam = 0 means not significant
    logprior <- function(beta, gam) {
      return(sum(dnorm(beta, mean = 0, sd = sqrt(0.01*(1-gam) + 100*gam), log = T)))
    }
    
    logposterior <- function(beta, gam) {
      return(loglikelihood(beta) + logprior(beta, gam))
    }
    
    proposalfunction <- function(beta, sd) {
      return(rnorm(length(beta), mean = beta, sd = sd))
    }
    
    step_beta <- rep(0.05, length(beta))
    ar_beta <- rep(0, length(beta))
    for (i in 1:length(beta)) {
      ## adaptive step
      beta_prop <- proposalfunction(beta[i], sd = step_beta[i])
      beta_proposal <- beta
      beta_proposal[i] <- beta_prop
      probab <- exp(logposterior(beta_proposal, gam) - logposterior(beta, gam))
      if (runif(1) < probab) {
        beta <- beta_proposal
        ar_beta[i] <- 1
      }
    }
    
    
    History[[iter]] <- list(gam = gam, beta = beta, lambda0 = lambda0, 
                            index_cluster = index_cluster, 
                            ar_beta = ar_beta)
    cat(" iteration:", iter,"\n")
  }
  
  list(Iterates = History)
}