# because the case is very rare, I generate case-control samples in multiple batches 
# sample ncases and ncontrols out of a population of size N 
# the parameters for the selection probability is gamma 
sample_obs_s <- function(ncase, ncontrol, sample_x_fun, sample_y_fun, selection_fun, N){
  X0 <- sample_x_fun(N) 
  Y0 <- sample_y_fun(X0)
  
  probs <- selection_fun(X0, Y0)[,1]
  D <- which(rbinom(N, 1, probs) == 1)
  
  cases <- sample(D, ncase, replace = F)
  controls <- sample((1:N)[-D], ncontrol, replace = F)
  
  Xcase <- X0[cases, ]
  Xcontrol <- X0[controls, ]
  Ycase <- Y0[cases]
  Ycontrol <- Y0[controls]
  
  list(X = rbind(Xcase, Xcontrol), 
       Y = c(Ycase, Ycontrol),
       D = rep(c(1, 0), time = c(ncase, ncontrol)))
}

sample_obs <- function(ncase, ncontrol, sample_x_fun, sample_y_fun, selection_fun){
  n <- ncase + ncontrol
  X <- matrix(NA, n, p)
  Y <- numeric(n)
  D <- numeric(n)
  times <- max(ncase, ncontrol) %/% 50 # each time sample ~50 cases and controls
  
  for(b in 1:times){
    cat(b, ",")
    newObs <- sample_obs_s(ncase/times, ncontrol/times, sample_x_fun, sample_y_fun, selection_fun, 20000)
    X[((b-1) * n / times + 1) : (b * n / times), ] <- newObs$X
    Y[((b-1) * n / times + 1) : (b * n / times)] <- newObs$Y
    D[((b-1) * n / times + 1) : (b * n / times)] <- newObs$D
  }
  
  list(X = X, 
       Y = Y,
       D = D)
}
