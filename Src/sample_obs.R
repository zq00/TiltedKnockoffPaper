# functions to sample observations 

## Sample covariates X from a multivariate Gaussian 
## the mean is 0 and users inpu the cholesky decomposition of the covariance
# INPUTS
# n -- number of obs.
# mu -- vector of the mean (default = zero vector)
# R -- choleskey decomposition of the covariance matrix
sample_gaussian <- function(n, mu = 0, R){
  p <- nrow(R)
  
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% R
  
  t(t(X) + mu)
}

## Sample covariates X from a  Markov chain
# the covariates are centered to have mean 0
# INPUTS
# n -- number of obs.
# pInit -- vector of the initial probability of each state
# Q -- an array of transition matrix (see the documentation of sampleDMC function from SNPknock package)
# OUTPUTS
# a matrix of size n * p 
sample_mc <- function(n, pInit, Q){ # sample X from a multivariate Gaussian
   as.matrix(sampleDMC(pInit, Q, n)) 
}

## Sample response Y from a linear model
## Y = X beta + beta0 + eps
## the error esp ~ N(0, 1)
# INPUTS
# X -- n * p covariate matrix
# beta -- coef for X
# beta0 -- intercept 
sample_y_linear <- function(X, beta, beta0){
  N <- nrow(X)
  Y <- X %*% beta + beta0 + rnorm(N, 0, 1)
  
  return(Y)
}

## Sample response Y from a logistic model 
## P(Y = 1) = 1-P(Y = -1) = 1 / (1 + exp(X beta + beta0))
# INPUTS
# X -- n * p centered covariate matrix
# beta -- coef for X
# beta0 -- intercept 
# Outputs
# A vector of +/- 1
sample_y_logistic <- function(X, beta, beta0){ 
  nx <- nrow(X);
  v <- X %*% beta + beta0
  nu <- 1 / (1 + exp(-v))
  
  y0 <- rbinom(nx, 1, prob = nu)
  2 * y0 - 1
}

## selection probability
## This returns P(S = 1 | X, Y) and P(S = 0 | X, Y)
# If more than one category, then returns P(D = d | X, Y) (this assumes that D includes only case categories, and add a last column for the probability of begin the control
# INPUTS
# X -- n * p covariate matrix
# Y -- vector of length n, responses 
# params -- a list with the following elements
#   model -- gaussian (probability is exp(-v^2 / 2)) or logistic (probability is 1/(1+exp(-v)))
#            v = gamma_0 + X gamma + Y gamma_y
#   gamma_0, gamma, gamma_y -- see above. they can be a vector or a matrix. 
pselect <- function(X, Y, params){
  
  v <- Y %*% t(params$gamma_y) + X %*% t(params$gamma_x)
  v <- t(t(v) + gamma_0)
  
  if(params$model == "logistic"){
    odds_ratio <- exp(v)
    
    # return the probability of being selected in each of the categories (include an additional column for being in the control)
    odds_ratio_addcontrol <- cbind(odds_ratio, 1) 
    probs <- t(apply(odds_ratio_addcontrol, 1, function(t) t / sum(t)))
    colnames(probs)[ncol(probs)] <- "Control"
  }else if(params$model == "gaussian"){
    probs <- exp(-1/2*v^2) 
    probs <- cbind(probs, 1-probs)
    colnames(probs) <- c("Case", "Control")
  }

  probs 
}


## generate a case-control sample 
# INPUTS 
# X -- n * p covariate matrix
# Y -- vector of length n, responses
# ncase -- number of cases
# ncontrol -- number of controls 
# probs - probability of being a case 
case_control_sample <- function(X, Y, ncase, ncontrol, N, probs){ 
  N <- nrow(X)
  Z <- rbinom(N, 1, probs)
  
  ind <- c(sample(which(Z == 1), ncase, replace = F),
           sample(which(Z == 0), ncontrol, replace = F))
    
  list(X = as.matrix(X[ind, ]), 
       Y = Y[ind], 
       probs = probs[ind], 
       D = Z[ind], 
       npos = sum(Z)
       )
}














