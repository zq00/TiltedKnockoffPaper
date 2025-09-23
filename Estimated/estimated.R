# Section 3.2 estimated selection probability

# input index  
args <- commandArgs(trailingOnly = T)
ind <- as.numeric(args[1])

# 1 -- Setup
# directories
out_dir <- "/result/estimated/"
coef_dir <- "/param/estimated/"

# load functions
source_code_dir <- "/Src/"  #The directory where all source code files are saved.
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

# 2 -- Functions
SecondOrderKO <- function(X, Y, D, sample_x_fun, selection_fun, B){
  d <- c(1, 0)
  XTilde <- matrix(NA, nrow(X), ncol(X))
  for(y in c(-1, 1)){
    est <- estimate_mucov_tilted(y, sample_x_fun, pselect_fun = selection_fun, B = 50) # estimate the mean and the covariance matrix of the tilted distribution 
    for(i in 1:2){
      index <- which(Y == y & D == d[i])
      
      if(sum(index) == 0) next; 
      
      cat(y, d[i], min(eigen(est$sigma_hat[i,,])$values), "\n")
      XTilde[index, ] <- create.gaussian(X[index,,drop = F], mu = est$mu_hat[i,], Sigma = est$sigma_hat[i,,], method="sdp")
    }
  }
  
  XTilde
}

## Extra: no confounding
X <- sample_mc(ncase + ncontrol, pInit, Q)
X_centered <- X - mu
Y <- sample_y_logistic(X_centered, beta, beta0)
XtildeKO <- knockoffDMC(X, pInit, Q)
resultKONoConfounding <- knockoff_selection(X, XtildeKO, Y, stat_fun, q, offset = offset, beta = beta)

# print result
cat("Knockoff without confounding: fdp (q = 0.3) = ", 
    round(resultKONoConfounding$result[6,2] / max(resultKONoConfounding$result[6,1], 1), 2), 
    ", power = ", round((resultKONoConfounding$result[6,1] - resultKONoConfounding$result[6,2]) / sum(beta!=0), 2) ,
    "\n")

# store result
fdpKONoConfounding <- data.frame(
  nselect = resultKONoConfounding$result[,1],
  nfd = resultKONoConfounding$result[,2],
  q = q,
  method = "ko_NoConfounding"
)


## 1 -- Sample obs
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

## 2 -- Knockoff with no adjustment
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

## 3 -- Second order knockoff with known gamma
XtildeSO <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun, B = 50)
resultSO <- knockoff_selection(obs$X, XtildeSO, obs$Y, stat_fun, q, offset = offset, beta = beta)

# print result
cat("SO Knockoff (known gamma): fdp (q = 0.3) = ", 
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

## 4 -- Second order knockoff with estimated gamma 
# 4.0 compute the prevalence 
D <- numeric(100*10000)
for(b in 1:100){
  X <- sample_x_fun(10000)
  Y <- sample_y_fun(X)
  probs <- selection_fun(X, Y)
  
  D[((b-1) * 10000 + 1) : (b * 10000)] <- rbinom(10000, 1, probs)
  
  cat(b, ",")
}
pi0 <- mean(D == 1)
# adjustment: for the estimated intercept, we *subtract* adj to get back to the intercept in the population
adj <- log(ncase * (1-pi0) / ncontrol / pi0)

# 4.1 lasso-lambda.min
XCentered <- obs$X - mu
cvfit <- cv.glmnet(cbind(XCentered, obs$Y), obs$D, family = "binomial", type.measure = "class")

# 4.1.1 lasso-lambda.min
gammaHat_lambdamin <- as.vector(coef(cvfit, s =  "lambda.min"))
gammaHat_lambdamin[1] <- gammaHat_lambdamin[1] - adj

selection_fun_lambdamin <- function(X, Y){
  XCentered <- X - mu
  pselect(XCentered, Y, params = list(
    model = "logistic",
    gamma_0 = gammaHat_lambdamin[1],
    gamma_x = t(gammaHat_lambdamin[2:(p+1)]), 
    gamma_y = gammaHat_lambdamin[p+2]
  ))
}
XtildeLambdaMin <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_lambdamin, B = 50)
resultLambdaMin <- knockoff_selection(obs$X, XtildeLambdaMin, obs$Y, stat_fun, q, offset = offset, beta= beta)

# print result
cat("SO Knockoff (lasso-lambda.min): fdp (q = 0.3) = ", 
    round(resultLambdaMin$result[6,2] / max(resultLambdaMin$result[6,1], 1), 2), "\n")
# store result
fdpLambdaMin <- data.frame(
  nselect = resultLambdaMin$result[,1],
  nfd = resultLambdaMin$result[,2],
  q = q,
  method = "lambda-min"
)

# 4.1.2 lasso-lambda.1se 
gammaHat_lambda1se <- as.vector(coef(cvfit, s =  "lambda.1se"))
gammaHat_lambda1se[1] <- gammaHat_lambda1se[1] - adj

selection_fun_lambdamin <- function(X, Y){
  XCentered <- X - mu
  pselect(XCentered, Y, params = list(
    model = "logistic",
    gamma_0 = gammaHat_lambda1se[1],
    gamma_x = t(gammaHat_lambda1se[2:(p+1)]), 
    gamma_y = gammaHat_lambda1se[p+2]
  ))
}
XtildeLambda1se <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_lambdamin, B = 50)
resultLambda1se <- knockoff_selection(obs$X, XtildeLambda1se, obs$Y, stat_fun, q,offset = offset, beta= beta)

# print result
cat("SO Knockoff (lasso-lambda.1se): fdp (q = 0.3) = ", 
    round(resultLambda1se$result[6,2] / max(resultLambda1se$result[6,1], 1), 2), "\n")
# store result
fdpLambda1se <- data.frame(
  nselect = resultLambda1se$result[,1],
  nfd = resultLambda1se$result[,2],
  q = q,
  method = "lambda-1se"
)

# 4.3 knockoff selection
# construct knockoff variables for X
XtildeKO <- knockoffDMC(obs$X, pInit, Q)

# insert the original Y for Y 
result_ko_gamma <- knockoff_selection(X = cbind(obs$X, obs$Y), 
                                      XTilde = cbind(XtildeKO, obs$Y), 
                                      Y = obs$D, stat_fun, q, offset = offset, beta = c(gamma_x, gamma_y))

# fit a logistic regression for each set of selected variable 
# I am making this too complicated, I only need to consider the case of (1) no X is selected (2) some X are selected
for(l in 1:6){
  gamma_name <- paste0("gammaHat_KO_", q[l])
  fdp_name <- paste0("fdpKO_", q[l])
  # load the result
  selected <- result_ko_gamma$selected_var[[l]]
  
  gamma_hat <- numeric(p + 2)
  
  # remove p+1 from selection
  if(sum(selected == (p+1)) == 1) selected <- selected[-which(selected == (p+1))]
  if(length(selected) == 0){
    fit <- glm(obs$D ~ obs$Y, family = binomial)
    gamma_hat[1] <- fit$coef[1]- adj
    gamma_hat[p + 2] <- fit$coef[2]
  }else{
    fit <- glm(obs$D ~  XCentered[, selected] + obs$Y, family = binomial)
    gamma_hat[1] <- fit$coef[1]- adj
    gamma_hat[selected + 1] <- fit$coef[2:(length(selected) + 1)]
    gamma_hat[p+2] <- tail(fit$coef, 1)
  }
  
  assign(gamma_name, gamma_hat)
  
  selection_fun_ko <- function(X, Y){
    XCentered <- X - mu
    pselect(XCentered, Y, params = list(
      model = "logistic",
      gamma_0 = gamma_hat[1],
      gamma_x = t(gamma_hat[2:(p+1)]), 
      gamma_y = gamma_hat[p+2]
    ))
  }
  Xtilde <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun, selection_fun_ko, B = 50)
  result <- knockoff_selection(obs$X, Xtilde, obs$Y, stat_fun, q,offset = offset, beta)
  
  # print result
  cat("SO Knockoff (KO), q = ",q[l] ,": fdp (q = 0.3) = ", 
      round(result$result[6,2] / max(result$result[6,1], 1), 2), "\n")
  # store result
  assign(fdp_name, data.frame(
    nselect = result$result[,1],
    nfd = result$result[,2],
    q = q,
    method = paste0("lambda-ko-", q[l])
  ))
}

# 4.4 logistic regression 
fit <- glm(obs$D ~ XCentered + obs$Y, family = binomial, x = T, y = T)

gammaHat_logistic <- fit$coef
gammaHat_logistic[1] <- gammaHat_logistic[1] - adj

selection_fun_logistic <- function(X, Y){
  XCentered <- X - mu
  pselect(XCentered, Y, params = list(
    model = "logistic",
    gamma_0 = gammaHat_logistic[1],
    gamma_x = t(gammaHat_logistic[2:(p+1)]), 
    gamma_y = gammaHat_logistic[p+2]
  ))
}
XtildeLogistic <- SecondOrderKO(obs$X, obs$Y, obs$D, sample_x_fun,selection_fun_logistic , B = 50)
resultLogistic <- knockoff_selection(obs$X, XtildeLogistic, obs$Y, stat_fun, q, offset = offset, beta)
# print result
cat("SO Knockoff (logistic): fdp (q = 0.3) = ", 
    round(resultLogistic$result[6,2] / max(resultLogistic$result[6,1], 1), 2), "\n")
# store result
fdpLogistic <- data.frame(
  nselect = resultLogistic$result[,1],
  nfd = resultLogistic$result[,2],
  q = q,
  method = "logistic"
)

## 6 -- Combining all the results and store in a data frame 
GammaHats <- cbind(gammaHat_lambdamin, 
                   gammaHat_lambda1se,
                   gammaHat_KO_0.05,
                   gammaHat_KO_0.1,
                   gammaHat_KO_0.15,
                   gammaHat_KO_0.2,
                   gammaHat_KO_0.25,
                   gammaHat_KO_0.3,
                   gammaHat_logistic)
colnames(GammaHats) <- c("lambda.min", "lambda.1se", "ko-0.05", "ko-0.1",
                         "ko-0.15", "ko-0.2", "ko-0.25", "ko-0.3", 
                         "logistic")

fdp <- rbind(fdpKONoConfounding, fdpKO, fdpSO, 
             fdpLambdaMin, fdpLambda1se,
             fdpKO_0.05, fdpKO_0.1, fdpKO_0.15, fdpKO_0.2, 
             fdpKO_0.25, fdpKO_0.3, fdpLogistic
)

write.table(GammaHats, paste0(out_dir, "result0514/GammaHat_", ind, ".txt"))
write.table(fdp, paste0(out_dir, "result0514/FDP_", ind, ".txt"), row.names = F, col.names = F)














