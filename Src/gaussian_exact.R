## Function to sample exact knockoffs for a mixture of Gaussian
## See section 3.2 and section 4 (see formula in the Appendix)
# INPUTS
# X -- n * p covariate matrix
# Y -- vector of length n, responses
# n -- a list with the following elements: 
#      npos -- # D = 1
#      nneg -- # D = 0
#      ncase -- # cases in the sample 
#      ncontrol -- # controls in the sample, ncase + ncontrol = nrow(X)
# params -- a list with the following elements:
#           Sigma -- p * p matrix, covariance of X 
#           mu -- vector of length p, mean vector of X 
#           gamma, gamma_y, gamma_0 -- P(D = 1 | X, Y) = exp(-v^2/2), v = X gamma + gamma_y Y + gamma_0
exact_knockoff <- function(X, Y, n, params){
  Sigma2 <- params$Sigma
  Theta2 <- solve(Sigma2)
  mu2 <- params$mu 
  
  Theta1 <- Theta2 + params$gamma %*% t(params$gamma)
  Sigma1 <- solve(Theta1)
  # each *row* corresponds to the mean of each obs in the first mixture
  mu1 <- - params$gamma_y * (Y %*% t(Sigma1 %*% params$gamma))  
  
  # Compute the diagonal terms (to avoid repeatedly solving for this term for each obs.)
  s1 <- create.solve_sdp(Sigma1)
  s2 <- create.solve_sdp(Sigma2)
  
  ntot <- nrow(X)
  p <- ncol(X)
  # w1[i] is the relative weight of the first mixture of the i-th obs
  logw1 <- log(n$ncase / n$npos - n$ncontrol / n$nneg) + (1/2 * c(t(params$gamma) %*% (Sigma1 %*% params$gamma) - 1) * params$gamma_y^2 * Y^2) 
  w2 <- n$ncontrol / n$nneg
  
  Xtilde <- matrix(NA, ntot, p)
  for(i in 1:ntot){
    # Compute the probability to be in each of the mixtures
    logd2 <-  - t(X[i,] - mu2) %*%  Theta2 %*% (X[i,] - mu2) / 2 
    logd1 <- - t(X[i,] - mu1[i, ]) %*%  Theta1 %*% (X[i,] - mu1[i, ]) / 2 
    q1 <- exp(logw1[i]  + logd1) / (exp(logw1[i]+ logd1) + w2 * exp(logd2)) 
    
    a <- rbinom(1, 1, q1) # if a = 1 then it's in the *first* mixture
    if(a == 0){
      Xtilde[i, ] <- create.gaussian(X[i,,drop = F], mu2, Sigma2,  diag_s = s2)
    }else{
      Xtilde[i, ] <- create.gaussian(X[i,,drop = F], mu1[i, ], Sigma1,  diag_s = s1)
    }
    
    if(i %% 50 == 0) cat(i, ",")
  }
  
  return(Xtilde)
}
