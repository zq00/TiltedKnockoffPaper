# 3/1/2023
# Getting training data for the deep knockoff machine 
# Each time I generate a small sample 
args <- commandArgs(trailingOnly = T)

# directories
coef_dir <- "/param/deep/"
out_dir <- "/result/deep/training_data/"

# load functions
source_code_dir <- "/Src/"  
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# model coef
beta <- scan(paste0(coef_dir, "beta.txt"))
beta0 <- 0

gamma_x <-  scan(paste0(coef_dir, "gamma.txt"))
gamma_y <- 2
gamma_0 <- -6

p <- length(beta)

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

selection_fun <- function(X, Y){
  pselect(X, Y, params = list(
    model = "logistic",
    gamma_0 = gamma_0,
    gamma_x = t(gamma_x), 
    gamma_y = gamma_y
  ))
}

# Functions
# rejection sampling 
# generate N samples from the tilted distribution where Y = y and D = d
# INPUTS
# N - (numeric) total number of obs. to sample from
# Y - (numeric) response (either 1 or -1)
# D - (numeric) case/control status (either 1 or 0)
rejection_sample <- function(N, y, d){
  X0 <- sample_mc(N, pInit, Q)
  X0centered <- X0 - mu
  
  # compute the selection probability (here Y is binary automatically)
  probs <- selection_fun(X0centered, rep(y, N))
  if(d == 1){
    Z <- rbinom(N, 1, probs[,1]) 
  }else{
    Z <- rbinom(N, 1, probs[,2]) 
  }
  ind <- which(Z == 1)
  
  # return the samples
  return(X0[ind, , drop = F])
}


# Sample a number of X and Y (in total generate 10,000 samples)
ntot <- 10000
X <- list()
for(i in 1:4) X[[i]] <- matrix(NA, nrow = ntot, ncol = p)

names <- c("control_neg", "case_neg", "control_pos", "case_pos")
i <- 0
for(y in c(-1, 1)){
  for(d in c(0, 1)){
    i <- i + 1
    
    flag <- FALSE
    xnew <- NULL
    while(!flag){
      new_sample <- rejection_sample(N = 100000, y = y, d = d)
      
      xnew <- rbind(xnew, new_sample)
      cat(y, d, nrow(xnew), "\n")
      if(nrow(xnew) >= ntot) flag <- T
    }
    X[[i]] <- xnew[1:ntot, ]
    write.table(X[[i]], paste0(out_dir, names[i], ".txt"),
                row.names = F, col.names = F)
  }
}

# Concatenate data correponding to each category (y, d combination) into a single data frame 
# B <- 105
# control_neg <- matrix(NA, nrow = ntot * B, ncol = p)
# case_neg <-  matrix(NA, nrow = ntot * B, ncol = p)
# control_pos <-  matrix(NA, nrow = ntot * B, ncol = p)
# case_pos <-  matrix(NA, nrow = ntot * B, ncol = p)
# 
# for(i in 1:B){
#   cat(i, ",")
#   new_control_neg <- read.table(paste0(out_dir, "control_neg_", i, ".txt"))
#   new_case_neg <- read.table(paste0(out_dir, "case_neg_", i, ".txt"))
#   new_control_pos <- read.table(paste0(out_dir, "control_pos_", i, ".txt"))
#   new_case_pos <- read.table(paste0(out_dir, "case_pos_", i, ".txt"))
#   
#   control_neg[((i-1)*ntot + 1): (i * ntot), ] <- as.matrix(new_control_neg)
#   case_neg[((i-1)*ntot + 1): (i * ntot), ] <- as.matrix(new_case_neg)
#   control_pos[((i-1)*ntot + 1): (i * ntot), ] <- as.matrix(new_control_pos)
#   case_pos[((i-1)*ntot + 1): (i * ntot), ] <- as.matrix(new_case_pos)
# }
# 
# # store data 
# write.table(control_neg, paste0(out_dir, "control_neg.txt"), col.names = F, row.names = F)
# write.table(case_neg, paste0(out_dir, "case_neg.txt"), col.names = F, row.names = F)
# write.table(control_pos, paste0(out_dir, "control_pos.txt"), col.names = F, row.names = F)
# write.table(case_pos, paste0(out_dir, "case_pos.txt"), col.names = F, row.names = F)



# sample the coefficients 
# p <- 50
# beta <- numeric(p)
# gamma <- numeric(p)
# nonnull_beta <- sample(1:p, 20, replace = F)
# nonnull_gamma <- sample(1:p, 20, replace = F)
# 
# beta[nonnull_beta] <- rnorm(20, 0, 0.2)
# gamma[nonnull_gamma] <- rnorm(20, 0, 0.2)
# write.table(beta, "/param/deep/beta.txt", row.names = F, col.names = F)
# write.table(gamma, "/param/deep/gamma.txt", row.names = F, col.names = F)











