## Functions to sample the second order tilted knockoffs

# Estimate the covariance matrix for the tilted distribution
# given the population distribution of X and the probability of being selected 
# INTPUTS
# y - (numeric) the response
# sample_x --  function to sample X (takes the no. of samples and outputs sampled obs.)
# pselect -- function to compute the selection probability (takes a matrix of covariates X and a response of Y and returns a matrix of selection probability)
# B - (numeric) number of samples (* 5000) used to estimate the mean/covariance of the tilted distribution (if B = 10, then use 50k samples) (default = 50)

estimate_mucov_tilted <- function(y, sample_x_fun, pselect_fun, B = 50){
  
  for(j in 1:B){
    x_new <- sample_x_fun(5000)
    probs <- pselect_fun(x_new, rep(y, 5000)) # this computes p1, ..., pd
    
    n_cat <- ncol(probs) # number of categories
    
    # initialize the arrays to store estimates
    if(j == 1){
      p <- ncol(x_new)
      mu <- array(NA, dim = c(n_cat, B, p))
      mu2 <- array(NA, dim = c(n_cat, B, p, p))
      probs_sums <- matrix(NA, n_cat, B) 
    }
    
    probs_sums[ ,j] <- colSums(probs) # sum of pselect
    mu[, j, ] <- t(apply(probs, 2, function(t) colSums(x_new * t) / sum(t)))
    
    for(k in 1:n_cat){
      P <- diag(probs[,k])
      mu2[k, j,,] <- as.matrix( t(x_new) %*% (P %*% x_new) / sum(probs[,k]))
    }
  }
  
  mu_hat <- matrix(NA, nrow = n_cat, ncol = p)
  mu2_hat <- array(NA, dim = c(n_cat, p, p))
  sigma_hat <- array(NA, dim = c(n_cat, p, p))
  for(k in 1:n_cat){
    # aggregate estimates
    mu_hat[k, ] <- t(probs_sums[k, ]) %*% mu[k, ,] / sum(probs_sums[k, ])
    mu2_hat[k, ,] <- apply(mu2[k, , ,] * probs_sums[k, ], c(2, 3), sum) / sum(probs_sums[k, ])
    sigma_hat[k, ,] <- mu2_hat[k, ,] - mu_hat[k, ] %*% t(mu_hat[k, ]) 
  }
  
  
  return(list(mu_hat = mu_hat, sigma_hat = sigma_hat))
}


