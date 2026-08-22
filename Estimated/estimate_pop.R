# Simulation in supplement 3.4
# Approximating the population distribution by multivariate Gaussian
rm(list = ls())

library(rlang, lib.loc = "/home/qianzhao_umass_edu/R/x86_64-pc-linux-gnu-library/4.5/")

# input index  
args <- commandArgs(trailingOnly = T)
ind <- as.numeric(args[1])

cat(ind, "\n")

# 1 -- Setup
# directories
out_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/estimated/"
coef_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/correct/coef/"

# load functions
source_code_dir <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/src/" 
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# model coef
beta <- scan(paste0(coef_dir, "beta.txt"))
beta0 <- 0

gamma_x <-  scan(paste0(coef_dir, "gamma.txt"))
gamma_x <- gamma_x 
gamma_y <- 2
gamma_0 <- -6

p <- length(beta)

# number of cases and controls
ncase <- 2000
ncontrol <- 2000

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
      
      cat(y, d[i], min(eigen(est$sigma_hat[i,,])$values), "\n")
      # check if the covariances are symmetric
      Sigma <- (est$sigma_hat[i,,] + t(est$sigma_hat[i,,])) / 2
      XTilde[index, ] <- create.gaussian(X[index,,drop = F], mu = est$mu_hat[i,], Sigma = Sigma, method="sdp")
    }
  }
  
  XTilde
}


# ----
# 1. Sample obs
# ----
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

obs <- sample_obs(ncase, ncontrol, sample_x_fun, 
                  sample_y_fun = sample_y_fun, 
                  selection_fun = selection_fun
)

write.table(obs$X, paste0(out_dir, "data/X_", ind, ".txt"))
write.table(obs$Y, paste0(out_dir, "data/Y_", ind,".txt"))
write.table(obs$D, paste0(out_dir, "data/D_", ind, ".txt"))

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

# store result
fdpKO <- data.frame(
  nselect = resultKO$result[,1],
  nfd = resultKO$result[,2],
  q = q,
  method = "ko"
)

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

# store result
fdpSO <- data.frame(
  nselect = resultSO$result[,1],
  nfd = resultSO$result[,2],
  q = q,
  method = "so"
)

# ----
# 4. second order knockoff: 
# - Estimated selection probability 
# - Estimating the distribution of X with Gaussian
# ----

# estimate prevalence 
D <- numeric(100*1000)
for(b in 1:100){
  X <- sample_x_fun(1000)
  Y <- sample_y_fun(X)
  probs <- selection_fun(X, Y)
  
  D[((b-1) * 1000 + 1) : (b * 1000)] <- rbinom(1000, 1, probs)
  
  cat(b, ",")
}
pi0 <- mean(D == 1)
# adjustment for intercept: for the estimated intercept, we *subtract* adj to get back to the intercept in the population
adj <- log(ncase * (1-pi0) / ncontrol / pi0)

fit_case <- CVglasso(obs$X[which(obs$D == 1), ], trace = 'none')
Sigma_hat_case <- fit_case$Sigma

fit_control <- CVglasso(obs$X[which(obs$D == 0), ], trace = 'none')
Sigma_hat_control <- fit_control$Sigma

Sigma_hat <- Sigma_hat_case * (1 - pi0) + Sigma_hat_control * pi0
R_hat <- chol(Sigma_hat)

mu_hat_case <- colMeans(obs$X[which(obs$D == 1), ])
mu_hat_control <- colMeans(obs$X[which(obs$D == 0), ])
mu_hat <- mu_hat_case * (1 - pi0) + mu_hat_control * pi0


sample_x_fun_gaussian <- function(n){
  # input the gaussian covariance matrix 
  sample_gaussian(n, mu = mu_hat, R = R_hat)
}


# Estimate selection probability
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

# Use the double second order approximation to generate knockoffs 
XtildeLambdaMin <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun_gaussian, selection_fun_lambdamin, B = 50)
resultLambdaMin <- knockoff_selection(obs$X, XtildeLambdaMin, obs$Y, stat_fun, q, offset = offset, beta= beta)

# print result
cat("SO Knockoff (est + gaussian approx of Px): fdp (q = 0.3) = ", 
    round(resultLambdaMin$result[6,2] / max(resultLambdaMin$result[6,1], 1), 2), 
    ", power = ", round((resultLambdaMin$result[6,1] - resultLambdaMin$result[6,2]) / sum(beta!=0), 2) ,
    "\n")

# store result
fdpLambdaMin <- data.frame(
  nselect = resultLambdaMin$result[,1],
  nfd = resultLambdaMin$result[,2],
  q = q,
  method = "gaussian"
)


# ----
# 4. Store results in a data frame 
# ----

fdp <- rbind(fdpKO, fdpSO, fdpLambdaMin)

write.table(fdp, paste0(out_dir, "FDP_", ind, ".txt"), row.names = F, col.names = F)









