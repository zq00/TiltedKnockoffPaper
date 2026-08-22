# Simulation in Supplement section 2 and paper section 3.2
rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
ncase <- ncontrol <- as.integer(args[1])
batch_num <- as.integer(args[2]) 

cat("=== n = ", ncase, "===\n")
cat("=== Batch ", batch_num, " ===\n")


# 1 -- Setup
# directories
out_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/correct/"
coef_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/correct/coef/"

# load functions
source_code_dir <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/src/"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# sample coefficients
p <- 200
beta <- numeric(p)
gamma_x <- numeric(p)

n_nonnull_gamma <- 40
n_nonnull_beta <- 40

nonnull_loc_beta <- sample(1:p, n_nonnull_beta, replace = F)
nonnull_loc_gamma <- sample(1:p, n_nonnull_gamma, replace = F)

beta[nonnull_loc_beta] <- rnorm(n_nonnull_beta, 0, 0.4)
gamma_x[nonnull_loc_gamma] <- rnorm(n_nonnull_gamma, 0, 0.4)

beta0 <- 0

gamma_y <- 2
gamma_0 <- -6

# number of cases and controls
n_vals <- c(200, 400, 800, 1600)

# MC
K=3;  # Number of possible states for each variable
# Marginal distribution for the first variable
pInit = rep(1/K,K)
# Create p-1 transition matrices
Qsmall <- matrix(c(0.5, 0.2, 0.3, 0.3, 0.5, 0.2, 0.2, 0.3, 0.5), nrow = K, ncol = K)
Q = array(NA,c(p-1, K, K))
for(j in 1:(p-1)) {
  Q[j,,] = Qsmall
}
mu <- sum(0:(K-1) * pInit)

# Knockoff selections
stat_fun <- get("stat.lasso_lambdasmax_bin")
q <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
offset <- 1 # offset = 1 means the knockoff+ procedure

# second order knockoff master function
SecondOrderKO <- function(X, Y, D, sample_x_fun, selection_fun, B){
  d <- c(1, 0)
  XTilde <- matrix(NA, nrow(X), ncol(X))
  for(y in c(-1, 1)){
    est <- estimate_mucov_tilted(y, sample_x_fun, pselect_fun = selection_fun, B = 50) # estimate the mean and the covariance matrix of the tilted distribution 
    for(i in 1:2){
      index <- which(Y == y & D == d[i])
      
      if(sum(index) == 0) next; 
      
      # cat(y, d[i], min(eigen(est$sigma_hat[i,,])$values), "\n")
      # check if the covariances are symmetric
      Sigma <- (est$sigma_hat[i,,] + t(est$sigma_hat[i,,])) / 2
      XTilde[index, ] <- create.gaussian(X[index,,drop = F], mu = est$mu_hat[i,], Sigma = Sigma, method="sdp")
    }
  }
  
  XTilde
}


sample_x_fun <- function(n){
  sample_mc(n, pInit, Q)
}

sample_y_fun <- function(X){
  XCentered <- X - mu
  sample_y_logistic(XCentered, beta, beta0)
}

selection_fun <- function(X, Y){
  XCentered <- X - mu
  pselect(XCentered, Y, params = list(
    model = "logistic",
    gamma_0 = gamma_0,
    gamma_x = t(gamma_x), 
    gamma_y = gamma_y
  ))
}

B <- 25 # number of repeatitions

# store results
# standard KO
nselect_ko <- matrix(0, nrow = B, ncol = length(q))
nfd_ko <- matrix(0, nrow = B, ncol = length(q))
# case/control knockoff
nselect_ko_case_control <- matrix(0, nrow = B, ncol = length(q))
nfd_ko_case_control <- matrix(0, nrow = B, ncol = length(q))

# SO tilted 
nselect_so <- matrix(0, nrow = B, ncol = length(q))
nfd_so <- matrix(0, nrow = B, ncol = length(q))
# SO tilted with estimated 
nselect_est <- matrix(0, nrow = B, ncol = length(q))
nfd_est <- matrix(0, nrow = B, ncol = length(q))
# ko + ipw 
nselect_ko_ipw <- matrix(0, nrow = B, ncol = length(q))
nfd_ko_ipw <- matrix(0, nrow = B, ncol = length(q))
# SO tilted with missing covariates 
nselect_missing <- matrix(0, nrow = B, ncol = length(q))
nfd_missing <- matrix(0, nrow = B, ncol = length(q))

# SO tilted with incorrect model
nselect_probit <- matrix(0, nrow = B, ncol = length(q))
nfd_probit <- matrix(0, nrow = B, ncol = length(q))

# SO tilted with xgboost model
nselect_xgboost <- matrix(0, nrow = B, ncol = length(q))
nfd_xgboost <- matrix(0, nrow = B, ncol = length(q))

for(b in 1:B){
  cat(b, ":\n")
  # ------
  # 1. sample obs
  # ------
  obs <- sample_obs(ncase, ncontrol, sample_x_fun, 
                    sample_y_fun = sample_y_fun, 
                    selection_fun = selection_fun)
  
  # ----
  # 2. Knockoff with no adjustment
  # ----
  XtildeKO <- knockoffDMC(obs$X, pInit, Q)
  resultKO <- knockoff_selection(obs$X, XtildeKO, obs$Y, stat_fun, q, offset = offset, beta = beta)
  
  # print result
  cat("Knockoff without adjustment: fdp (q = 0.3) = ", 
      round(resultKO$result[6,2] / max(resultKO$result[6,1], 1), 2), 
      ", power = ", round((resultKO$result[6,1] - resultKO$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  
  nselect_ko[b, ] <- resultKO$result[ ,1]
  nfd_ko[b, ] <- resultKO$result[ ,2]
  
  # ----
  # 3. Second order knockoff with known selection probability
  # ----
  XtildeSO <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun, B = 50)
  resultSO <- knockoff_selection(obs$X, XtildeSO, obs$Y, stat_fun, q, offset = offset, beta = beta)
  
  # print result
  cat("SO Knockoff (known model): fdp (q = 0.3) = ", 
      round(resultSO$result[6,2] / max(resultSO$result[6,1], 1), 2), 
      ", power = ", round((resultSO$result[6,1] - resultSO$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  
  nselect_so[b, ] <- resultSO$result[ ,1]
  nfd_so[b, ] <- resultSO$result[ ,2]
  
  # ----
  # 4. Case/control knockoff 
  # ----
  
  X_case <- obs$X[which(obs$D == 1), ]
  X_control <- obs$X[which(obs$D == 0), ]
  Y_case <- obs$Y[which(obs$D == 1)]
  Y_control <- obs$Y[which(obs$D == 0)]
  
  Xtilde_case <- create.second_order(X = X_case, method = "sdp")
  Xtilde_control <- create.second_order(X =X_control, method = "sdp")
  
  result_ko_case_control <- knockoff_selection(rbind(X_case, X_control), 
                                               rbind(Xtilde_case, Xtilde_control), 
                                               c(Y_case, Y_control), 
                                               stat_fun, q, offset, beta = beta)
  
  nselect_ko_case_control[b, ] <- result_ko_case_control$result[ ,1]
  nfd_ko_case_control[b, ] <- result_ko_case_control$result[ ,2]
  
  cat("Case/control Knockoff: fdp (q = 0.3) = ", 
      round(result_ko_case_control$result[6,2] / max(result_ko_case_control$result[6,1], 1), 2), 
      ", power = ", round((result_ko_case_control$result[6,1] - result_ko_case_control$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  
  
  # ----
  # 5. Second order knockoff with estimated selection probability
  # ----
  
  # estimate prevalence 
  D <- numeric(100*1000)
  for(k in 1:100){
    X <- sample_x_fun(1000)
    Y <- sample_y_fun(X)
    probs <- selection_fun(X, Y)
    
    D[((k-1) * 1000 + 1) : (k * 1000)] <- rbinom(1000, 1, probs)
    
    # cat(k, ",")
  }
  pi0 <- mean(D == 1)
  
  # adjustment for intercept: for the estimated intercept, we *subtract* adj to get back to the intercept in the population
  adj <- log(ncase * (1-pi0) / ncontrol / pi0)
  
  XCentered <- obs$X - mu
  cvfit <- cv.glmnet(cbind(XCentered, obs$Y), obs$D, family = "binomial", type.measure = "class")
  coef_lasso   <- coef(cvfit, s = "lambda.min")
  selected_idx <- which(coef_lasso[-c(1, p+2)] != 0)  
  
  fit <- glm(obs$D ~ cbind(XCentered[,selected_idx], obs$Y), family = "binomial")
  
  gamma_0_lambdamin <- coef(fit)[1] - adj # adjusted coefficient
  gamma_x_lambdamin <- numeric(p)
  gamma_x_lambdamin[selected_idx] <- as.vector(coef(fit)[2:(length(selected_idx) + 1)])
  gamma_y_lambdamin <- tail(coef(fit), 1)   
  
  selection_fun_lambdamin <- function(X, Y){
    XCentered <- X - mu
    pselect(XCentered, Y, params = list(
      model = "logistic",
      gamma_0 = gamma_0_lambdamin, 
      gamma_x = t(gamma_x_lambdamin), 
      gamma_y = gamma_y_lambdamin
    ))
  }
  
  XtildeLambdaMin <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_lambdamin, B = 50)
  resultLambdaMin <- knockoff_selection(obs$X, XtildeLambdaMin, obs$Y, stat_fun, q, offset = offset, beta= beta)
  
  # print result
  cat("SO Knockoff (estimated model): fdp (q = 0.3) = ", 
      round(resultLambdaMin$result[6,2] / max(resultLambdaMin$result[6,1], 1), 2), 
      ", power = ", round((resultLambdaMin$result[6,1] - resultLambdaMin$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  # store result
  nselect_est[b, ] <- resultLambdaMin$result[ ,1]
  nfd_est[b, ] <- resultLambdaMin$result[ ,2]
  
  # ------
  # 6. Standard knockoffs + IPW (with estimated selection probability)
  # ------
  
  est_selection_probs <- selection_fun_lambdamin(obs$X, obs$Y)
  w_ipw <- ifelse(obs$D == 1, pi0 / est_selection_probs, (1-pi0) / (1-est_selection_probs))
  
  # knockoff selection
  result_ko_ipw <- knockoff_selection(obs$X, XtildeKO, obs$Y, stat_fun, q, offset, beta = beta, weights = w_ipw)
  
  # print result
  cat("Knockoff + IPW (estimated model): fdp (q = 0.3) = ", 
      round(result_ko_ipw$result[6,2] / max(result_ko_ipw$result[6,1], 1), 2), 
      ", power = ", round((result_ko_ipw$result[6,1] - result_ko_ipw$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  # store result
  nselect_ko_ipw[b, ] <- result_ko_ipw$result[ ,1]
  nfd_ko_ipw[b, ] <- result_ko_ipw$result[ ,2]
  
  # ----
  # 7. Correct logistic model but with missing covariates
  # ----
  
  cvfit <- cv.glmnet(cbind(XCentered[,1:(p/2)], obs$Y), obs$D, family = "binomial", type.measure = "class")
  
  coef_missing   <- coef(cvfit, s = "lambda.min")
  selected_idx_missing <- which(coef_missing[-c(1, p/2+2)] != 0)  
  
  fit_missing <- glm(obs$D ~ cbind(XCentered[,selected_idx_missing], obs$Y), family = "binomial")
  
  gamma_0_missing <- coef(fit_missing)[1] - adj # adjusted coefficient
  gamma_x_missing <- numeric(p)
  gamma_x_missing[selected_idx_missing] <- as.vector(coef(fit_missing)[2:(length(selected_idx_missing) + 1)])
  gamma_y_missing <- tail(coef(fit_missing), 1)   
  
  selection_fun_missing <- function(X, Y){
    XCentered <- X - mu
    pselect(XCentered, Y, params = list(
      model = "logistic",
      gamma_0 = gamma_0_missing, 
      gamma_x = t(gamma_x_missing), 
      gamma_y = gamma_y_missing
    ))
  }
  
  XtildeMissing <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_missing, B = 50)
  resultMissing <- knockoff_selection(obs$X, XtildeMissing, obs$Y, stat_fun, q, offset = offset, beta= beta)
  
  # print result
  cat("SO Knockoff (missing covariates): fdp (q = 0.3) = ", 
      round(resultMissing$result[6,2] / max(resultMissing$result[6,1], 1), 2), 
      ", power = ", round((resultMissing$result[6,1] - resultMissing$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  # store result
  nselect_missing[b, ] <- resultMissing$result[ ,1]
  nfd_missing[b, ] <- resultMissing$result[ ,2]
  
  # ----
  # 8. Incorrect model (Probit)
  # ----
  # weighted GLM to estimate probit regression
  w <- ifelse(obs$D == 1, pi0, 1-pi0)
  
  cvfit_probit <- cv.glmnet(cbind(XCentered, obs$Y), obs$D, 
                            weights = w, 
                            family =  binomial(link = "probit"))
  
  coef_probit   <- coef(cvfit_probit, s = "lambda.min")
  selected_idx_probit <- which(coef_probit[-c(1, p+2)] != 0)  
  
  fit_probit <- glm(obs$D ~ cbind(XCentered[,selected_idx_probit], obs$Y),
                    family =  binomial(link = "probit"), 
                    weights = w)
  
  gamma_0_probit <- coef(fit_probit)[1] 
  gamma_x_probit <- numeric(p)
  gamma_x_probit[selected_idx_probit] <- as.vector(coef(fit_probit)[2:(length(selected_idx_probit) + 1)])
  gamma_y_probit <- tail(coef(fit_probit), 1)   
  
  selection_fun_lambdamin_probit <- function(X, Y){
    XCentered <- X - mu
    pselect(XCentered, Y, params = list(
      model = "probit",
      gamma_0 = gamma_0_probit,
      gamma_x = t(gamma_x_probit), 
      gamma_y = gamma_y_probit
    ))
  }
  
  XtildeLambdaMin_probit <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_lambdamin_probit, B = 50)
  resultLambdaMin_probit <- knockoff_selection(obs$X, XtildeLambdaMin_probit, obs$Y, stat_fun, q, offset = offset, beta= beta)
  
  # print result
  cat("Incorrect model (probit): fdp (q = 0.3) = ", 
      round(resultLambdaMin_probit$result[6,2] / max(resultLambdaMin_probit$result[6,1], 1), 2), 
      ", power = ", round((resultLambdaMin_probit$result[6,1] - resultLambdaMin_probit$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  # store result
  nselect_probit[b, ] <- resultLambdaMin_probit$result[ ,1]
  nfd_probit[b, ] <- resultLambdaMin_probit$result[ ,2]
  
  
  # ----
  # 9. Flexible model: XGBoost 
  # ----
  
  # Fit an XGBoost model and adjust for the intercept by transforming the estimated probability
  dat <- cbind(XCentered, obs$Y)
  dimnames(dat)[[2]] <- c(paste0("X_", 1:p), "Y")
  
  ind_test <- sample(1:nrow(dat), ncase * 0.4, replace = F)
  ind_train <- (1:nrow(dat))[-ind_test]
  
  dtrain <- xgb.DMatrix(data = dat[ind_train, ], label = obs$D[ind_train])
  dtest <- xgb.DMatrix(data = dat[ind_test, ], label = obs$D[ind_test])
  
  params <- list(
    max_depth = 7,
    eta = 0.1,
    subsample = 0.5,
    objective = "binary:logistic"  
  )
  
  fit_xgboost <- xgb.train(
    data = dtrain,
    nrounds = 400,
    early_stopping_rounds = 10,
    evals = list(train = dtrain, eval = dtest),  
    params = params,
    verbose = 0
  )
  
  fitted_probs_xgboost <- predict(fit_xgboost, newdata = dat)
  # adjust the fitted probs
  fitted_probs_xgboost_adj <- 1 / (1 + (1 - pi0) / pi0 * (1-fitted_probs_xgboost) / fitted_probs_xgboost)
  
  selection_fun_xgboost <- function(X, Y){
    XCentered <- X - mu
    pselect(XCentered, Y, params = list(model = "xgboost"),
            fit = fit_xgboost,
            pi0 = pi0
    )
  }
  
  est_selection_probs_xgboost <- selection_fun_xgboost(obs$X, obs$Y)
  # true_selection_probs <- selection_fun(obs$X, obs$Y)
  # plot(true_selection_probs[,1], fitted_probs_xgboost_adj)
  
  XtildeXGBoost <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_xgboost, B = 50)
  resultXGBoost <- knockoff_selection(obs$X, XtildeXGBoost, obs$Y, stat_fun, q, offset = offset, beta= beta)
  
  # print result
  cat("XGBoost: fdp (q = 0.3) = ", 
      round(resultXGBoost$result[6,2] / max(resultXGBoost$result[6,1], 1), 2), 
      ", power = ", round((resultXGBoost$result[6,1] - resultXGBoost$result[6,2]) / sum(beta!=0), 2) ,
      "\n")
  # store result
  nselect_xgboost[b, ] <- resultXGBoost$result[ ,1]
  nfd_xgboost[b, ] <- resultXGBoost$result[ ,2]
}

# store result
save.image(file = paste0(out_dir, "n_", ncase, "_batch_", batch_num, "_V2.RData"))




