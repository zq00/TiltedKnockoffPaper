# The metro sampling algorithm
# Input: 
# x --  vector of length p. Observed variables.
# xprop -- vector of length p. Proposals.
# y -- numeric. Observed outcome.
# probsX -- vector of length p. q(x) probability under the proposal distribution
# probsXProp -- vector of length p. q(xprop) probability under the proposal distribution for the proposals.
# neighbors -- list of length p. Neighborhood of Xj.
# Output:
# xtilde -- vector of length p. Knockoffs. 
metro <- function(x, xprop,  probsX, probsXProp, y, d){
  p <- length(x)
  acc <- numeric(p) # vector recording whether each proposal is accepted
  
  for(j in ordered_nodes$order){
    cat(j, ",")
    paccept <- Compute_AcceptProb(x, xprop, probsX, probsXProp, affected_vars, j, index = rep(0, p) ,acc, y, d)
    acc[j] <- rbinom(1, 1, as.numeric(paccept))
  }
  
  xtilde <- x
  xtilde[which(acc == 1)] <- xprop[which(acc == 1)]
  return(xtilde)
}

# Compute the acceptance probability of an individual proposal
# Input
# x --  vector of length p. Observed variables.
# xprop -- vector of length p. Proposals.
# y -- numeric. Observed outcome.
# probsX -- vector of length p. q(x) probability under the proposal distribution
# probsXProp -- vector of length p. q(xprop) probability under the proposal distribution for the proposals.
# affected variables - list of length p. Which variable's acceptance probability is affected by the j-th proposal? 
# j - Numeric. should we accept the j-th proposal?
# index -- vector of length p. if 1 then that value is the proposal otherwise it is xj
# y -- Numeric (1 or -1). Observed outcome.
# d -- Numeric (0 or 1). Case-control status
# acc -- Binary vector of length p. Is the j-th proposal accepted? 
Compute_AcceptProb <- function(x, xprop,probsX, probsXProp, affected_vars, j, index, acc, y, d){
  
  if(x[j] == xprop[j]) return(0) # do not accept if the proposal and the observations are equal 
  
  nbdj <- computed_probs[[j]] 
  # check if the acceptance probability has been computed 
  if(is.null(nbdj$nbd)){
    # have we computed an acceptance probability? 
    if(!is.na(nbdj$acc_prob)){
      return(nbdj$acc_prob)
    }
  }else{
    ind <- which(nbdj$values[,1] == paste0(index[nbdj$nbd], collapse = ""))
    if(!is.na(nbdj$values[ind, 2])){
      return(nbdj$values[ind, 2])
    }
  }
  
  index0 <- index1 <- index
  index0[j] <- 0 # the denominator
  index1[j] <- 1 # the numerator
  if(j != ordered_nodes$order[1]){
    ii <- which(ordered_nodes$order == j)
    index1[ ordered_nodes$order[1:(ii-1)]] <- 0; 
    index0[ ordered_nodes$order[1:(ii-1)]] <- 0;
  }
  
  p <- length(x)
  
  x1 <- x0 <- x
  x1[which(index1 == 1)] <- xprop[which(index1 == 1)]
  if(sum(index0 == 1) > 0){
    x0[which(index0 == 1)] <- xprop[which(index0 == 1) ]
  }

  LogProb1 <- ComputeLogDensity(x1, y, d) + log(probsX[j])
  LogProb0 <- ComputeLogDensity(x0, y, d) + log(probsXProp[j])
  
  if(length(affected_vars[[j]]) > 0){
    for(k in affected_vars[[j]]){
      pnew1 <- as.numeric(Compute_AcceptProb( x, xprop, probsX, probsXProp, affected_vars, k, index1, acc, y, d))
      pnew0 <- as.numeric(Compute_AcceptProb( x, xprop, probsX, probsXProp, affected_vars, k, index0, acc,y, d))
      
      LogProb1 <- LogProb1 + ifelse(acc[k] == 1, log(pnew1 + 1e-10), log(1 + 1e-10 -pnew1))
      LogProb0 <- LogProb0 + ifelse(acc[k] == 1, log(pnew0 + 1e-10), log(1 + 1e-10 -pnew0))
    }
  }
  
  paccept <- min(1, exp(LogProb1 - LogProb0))
  
  if(is.null(nbdj$nbd)){
    computed_probs[[j]]$acc_prob <<- paccept
  }else{
    ind <- which(nbdj$values[,1] == paste0(index[nbdj$nbd], collapse = ""))
    
    computed_probs[[j]]$values[ind, 2] <<- paccept
  }
  
  return(paccept)
}
