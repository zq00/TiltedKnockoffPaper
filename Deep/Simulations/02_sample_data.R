## 02_sample_data.R
## For simulation replicate b:
##   - generate population (X, Y) from the model
##   - apply case-control selection
##   - estimate selection model
##   - generate training data for the deep knockoff machine
##
## Usage: Rscript 02_sample_data.R <b>
##   where b = 1, 2, ..., n_sim

rm(list = ls())
source("/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/src_V2/00_config.R")

# ============================================================
# 0. Parse command-line argument
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 02_sample_data.R <b>")
b <- as.integer(args[1])
cat("=== Simulation replicate:", b, "===\n")

# ============================================================
# 1. Load shared coefficients
# ============================================================
load(coef_file)

# ============================================================
# 2. Load required packages and source files
# ============================================================
source(paste0(path_src, "load_files.R"))
source(paste0(path_src, "sample_train.R"))

# ============================================================
# 3. Set up block covariance structure 
# ============================================================
if (covariate == "gaussian") {
  rho   <- 0.6
  bsize <- 20
  Sigma_small <- toeplitz(rho^(0:(bsize - 1)))
  R_small     <- chol(Sigma_small)
  
  dd_sigma <- list()
  dd       <- list()
  for (i in 1:(p %/% bsize)) {
    dd[[i]]       <- R_small
    dd_sigma[[i]] <- Sigma_small
  }
  if (p %% bsize != 0) {
    dd[[length(dd) + 1]]       <- R_small[1:(p %% bsize), 1:(p %% bsize)]
    dd_sigma[[length(dd_sigma) + 1]] <- Sigma_small[1:(p %% bsize), 1:(p %% bsize)]
  }
  R     <- as.matrix(bdiag(dd))
  Sigma <- as.matrix(bdiag(dd_sigma))
  
  sample_x <- sample_gaussian
}

# ============================================================
# 4. Sample population data (X0, Y0)
# ============================================================
set.seed(236 * b)   # different seed per replicate

X0 <- sample_x(N, R = R)

if (response == "mean") {
  v <- X0[, non_null_mean] %*% beta_mean
  Y0 <- v + rnorm(N, 0, 1)
} else if (response == "interaction") {
  v <- X0[, non_null_mean] %*% beta_mean +
    (X0[, non_null_interaction[, 1]] * X0[, non_null_interaction[, 2]]) %*% beta_interaction
  Y0 <- v + rnorm(N, 0, 1)
}

# print SNR
cat("SNR = ", sd(v))

# ============================================================
# 5. Compute selection probabilities and sample observed data
# ============================================================
gamma_y <- 1
gamma_0 <- -6

if (selection == "interaction") {
  mean_effect <- gamma_0 + Y0 * gamma_y + X0[, non_null_select_mean] %*% gamma_mean
  interaction_effect <-
    colSums(t(X0[, non_null_select_interaction_y] * as.numeric(Y0)) * gamma_interaction_y) +
    colSums(t(X0[, non_null_select_interaction_x[, 1]] * X0[, non_null_select_interaction_x[, 2]]) * gamma_interaction_x)
  v <- as.numeric(mean_effect) + interaction_effect
} else if (selection == "mean") {
  mean_effect <- gamma_0 + Y0 * gamma_y + X0[, non_null_select_mean] %*% gamma_mean
  v <- as.numeric(mean_effect)
}

probs   <- logistic(v)
is_case <- rbinom(N, 1, probs)
cat("Total cases in population:", sum(is_case), "\n")

ind_case    <- sample(which(is_case == 1), n_case, replace = FALSE)
ind_control <- sample(which(is_case == 0), n_control, replace = FALSE)
ind         <- c(ind_case, ind_control)
status      <- rep(c("case", "control"), times = c(n_case, n_control))

X_obs <- X0[ind, ]
Y_obs <- Y0[ind]
d     <- ifelse(status == "case", 1, 0)

# ============================================================
# 6. Estimate selection model (logistic regression)
# ============================================================
cvfit <- cv.glmnet(cbind(X_obs, Y_obs), y = status,
                   family = "binomial", type.measure = "deviance")

selected_x <- setdiff(which(coef(cvfit, s = cvfit$lambda.min)[-1] != 0), p + 1)

cat("# selected =", length(selected_x), "\n",
    "FPP =", sum(selected_x %in% non_null_gamma) / max(length(selected_x), 1), "\n",
    "Power =", sum(selected_x %in% non_null_gamma) / length(non_null_gamma), "\n")

fit <- glm(d ~ cbind(X_obs[, selected_x], Y_obs), family = "binomial")

# update the intercept 
pi0 <- mean(is_case)
adj <- log(n_case * (1-pi0) / n_control / pi0)

# ============================================================
# 7. Generate training data for the deep knockoff machine
# ============================================================
cat("Generating training data (case) logistic...\n")
training_data_logistic_case <- sample_train(n_train, Y = Y_obs[which(d == 1)],
                                            method = "logistic", is_case = TRUE,
                                            empirical = FALSE)

cat("Generating training data (control) logistic...\n")
training_data_logistic_control <- sample_train(n_train, Y = Y_obs[which(d == 0)],
                                               method = "logistic", is_case = FALSE,
                                               empirical = FALSE)

# ============================================================
# 8. Estimate selection model (XGBoost)
# ============================================================

dat <- cbind(X_obs, Y_obs)
dimnames(dat)[[2]] <- c(paste0("X_", 1:p), "Y")

ind_test <- sample(1:nrow(dat), 2000, replace = F)
ind_train <- (1:nrow(dat))[-ind_test]

dtrain <- xgb.DMatrix(data = dat[ind_train, ], label = d[ind_train])
dtest <- xgb.DMatrix(data = dat[ind_test, ], label = d[ind_test])
watchlist <- list(train = dtrain, eval = dtest)

fit_xgboost <- xgb.train(
  data = dtrain,
  max.depth = 8,
  eta = 0.1,
  nrounds = 400,
  watchlist = watchlist,
  early_stopping_rounds = 10,
  objective = "binary:logistic", 
  # eval_metric = "error",   
  subsample = 0.5
)

# ============================================================
# 9. Generate training data for the deep knockoff machine (not run)
# ============================================================
# cat("Generating training data (case) XGBoost...\n")
# training_data_xgboost_case <- sample_train(n_train, Y = Y_obs[which(d == 1)], 
#                                             method = "xgboost", is_case = TRUE)
# 
# cat("Generating training data (control) XGBoost ...\n")
# training_data_xgboost_control <- sample_train(n_train, Y = Y_obs[which(d == 0)], 
#                                                method = "xgboost", is_case = FALSE)

# ============================================================
# 8. Save everything needed by the training and result scripts
# ============================================================
out_dir <- sim_dir(b)

save(
  # observed data
  X_obs, Y_obs, d, status, Sigma,
  # selection model objects
  fit, selected_x, cvfit, adj, 
  # training data (logistic)
  training_data_logistic_case, training_data_logistic_control,
  # training data (xgboost)
  # training_data_xgboost_case, training_data_xgboost_control,
  # coefficients (for convenience in result script)
  non_null_loc, non_null_gamma, R,
  file = file.path(out_dir, "sample_data.RData")
)

cat("Sample data saved to:", file.path(out_dir, "sample_data.RData"), "\n")

