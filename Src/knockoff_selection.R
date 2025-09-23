# Functions related to knockoff selections 

## Compute knockoff selections (for single knockoffs)
## return the selected variables at a vector of target FDR levels
## return number of false selections if the true model coefficients are provided  
# INPUTS
# X -- n * p covariate matrix 
# XTilde -- n * p knockoff matrix 
# Y -- vector of length n, responses
# stat_fun -- a function to compute the knockoff scores (examples are functions starting with "stat" in the knockoff package)
# q -- a vector of target FDR levels 
# offset -- either 0 or 1 (see "knockoff.filter" in the knockoff package) (default = 1, the knockoff+ procedure)
# beta -- (optional) vector of length p, true model coefficients
knockoff_selection <- function(X, XTilde, Y, stat_fun, q, offset = 1, beta = NULL){
  W <- stat_fun(X, XTilde, Y) # knockoff scores
  
  selected_var <- list()  # a list containing indices of the selected variables at each FDR level specified in q
  result <- matrix(NA, nrow = length(q), ncol = 2) # columns contain the FDR and power
  
  for(i in 1:length(q)){
    thres <-  knockoff.threshold(W, fdr=q[i], offset=offset)
    selected <-  which(W >= thres)
    selected_var[[i]] <- selected 
    
    if(!is.null(beta)){ # if true coef is provided
      nfd <- sum(beta[selected] == 0)
      result[i, ] <- c(length(selected), nfd)
    }
  }
  
  return(
    list(selected_var = selected_var,
         result = result)
  )
}