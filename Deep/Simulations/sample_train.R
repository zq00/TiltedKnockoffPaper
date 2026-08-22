# sample training data

# Inputs:
# size of the training data 
# Y values for constructing the training set 
# Output:
# A training set of a given size

sample_train <- function(n_train, Y, method, is_case, X = NULL, empirical = FALSE, px = NULL){
  
  # estimate distribution for Y
  density_est <- density(Y)
  # estimated cdf
  t <- seq(min(Y)-0.1, max(Y)+0.1, by = 0.1)
  cdf_fun <- CDF(density_est)
  cdf_est <- cdf_fun(t)
  
  # sample y values for the training set 
  y_train <- sapply(runif(n_train), function(s) {
    x <- which(cdf_est >= s)
    if(length(x) == 0){
      max(Y) + 0.1
    }else{
      t[min(x)]
    }
  })
  
  # Sample X 
  x_train <- matrix(NA, nrow = n_train, ncol = p)
  bad_ind <- NULL # indices where we cannot sample x 
  n0 <- 2000 # sample n0 x at each step 
  i <- 1
  
  while(i <= n_train){
    if(i %% 100 == 0) cat(i / 100, ",")
    if(covariate == "gaussian" & empirical == FALSE){
      x0 <- sample_x(n0, R = R)
    }else if(empirical == TRUE){
      n <- nrow(X)
      ind_new <- sample(1:n, n0, prob = px, replace = T)
      x0 <- X[ind_new, ]
    }
    
    # compute selection probability
    if(method == "logistic"){
      v <- cbind(x0[,selected_x], y_train[i]) %*% coef(fit)[-1] + coef(fit)[1] - adj
      probs <- logistic(v) 
    }else if(method == "xgboost"){
      probs <-pselect(x0, rep(y_train[i], n0), 
                      params = list(model = "xgboost"), fit = fit_xgboost,
                      pi0 = pi0)[,1]
    }
    
    if(is_case == FALSE) probs <- 1 - probs
    
    U <- runif(n0)
    ind <- mapply(function(t1, t2) t1 > t2, probs, U)
    if(sum(ind) == 0){
      k <- k + 1
      if(k == 10){
        bad_ind <- c(bad_ind, i)
        i <- i+1
        cat("# bad obs = ", length(bad_ind), "\n")
      }
    }else{
      x_train[i, ] <- x0[sample(which(ind == T), 1), ]
      i <- i + 1
      k <- 0
    }
  }
  
  if(length(bad_ind) > 0){
    list(x_train = x_train[-bad_ind, ],
         y_train = y_train[-bad_ind])
  }else{
    list(x_train = x_train,
         y_train = y_train)
  }
   
}
  
  