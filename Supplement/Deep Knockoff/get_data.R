# 3/1/2024
# Sample observations for the deep knockoff 
# 1 -- Setup

# directories
out_dir <- "/result/metro/data/"
coef_dir <- "/param/metro/"

# load functions
source_code_dir <- "/Src/"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# sample model coefficients
p <- 50
# beta <- numeric(p)
# gamma <- numeric(p)
# non_null_beta <- sample(1:p, 15)
# non_null_gamma <- sample(1:p[-non_null_beta], 15) # remove the overlapping ones
# 
# beta[non_null_beta] <- rnorm(15, 0, 0.25)
# gamma[non_null_gamma] <- rnorm(15, 0, 0.25)
# write.table(beta, paste0(coef_dir, "beta.txt"), row.names = F, col.names = F)
# write.table(gamma, paste0(coef_dir, "gamma.txt"), row.names = F, col.names = F)


# model coef
beta <- scan(paste0(coef_dir, "beta.txt"))
beta0 <- 0

gamma_x <-  scan(paste0(coef_dir, "gamma.txt"))
gamma_y <- 2
gamma_0 <- -4

p <- length(beta)

# number of cases and controls
ncase <- 1000
ncontrol <- 1000

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


# Sample data
sample_x_fun <- function(n){
  sample_mc(n, pInit, Q)
}

sample_y_fun <- function(X){
  XCentered <- X - mu
  sample_y_logistic(XCentered, beta, beta0)
}

selection_fun <- function(X, Y){
  XCentered <- X - mu
  pselect(selection_fun, Y, params = list(
    model = "logistic",
    gamma_0 = gamma_0,
    gamma_x = t(gamma_x), 
    gamma_y = gamma_y
  ))
}


B <- 250 # 500 repetitions 
# each time get a sample, and store the scenarios for each combination of y and d separately 
# print out whether each of the instances exist 

names <- c("control_neg", "case_neg", "control_pos", "case_pos")

dat <- NULL 
for(b in c(401:450)){
  cat("--", b,":\n")
  
  obs <- sample_obs(ncase, ncontrol, sample_x_fun, sample_y_fun, selection_fun)
  
  i <- 0
  for(y in c(-1, 1)){
    for(d in c(0, 1)){
      i <- i + 1
      ind <- which(obs$Y == y & obs$D == d)
      
      cat(y, d, length(ind), "\n")
      if(length(ind) != 0) { # store the data table 
        write.table(obs$X[ind, , drop = F], file = paste0(out_dir, names[i], "_", b , ".txt"), row.names = F, col.names = F)
      }else{# store the situation where there's no instance 
        dat <- rbind(dat, c(b, y, d))
      }
    }
  }
}


# combine data together
names <- c("control_neg", "case_neg", "control_pos", "case_pos")
for(b in 251:400){
  cat(b, ",")
  n <- numeric(4)
  X <- NULL
  for(i in 1:4){
    xnew <- as.matrix(read.table(paste0(out_dir ,names[i], "_", b, ".txt")))
    X <- rbind(X, xnew)
    n[i] <- nrow(xnew)
  }
  Y <- rep(c( -1, 1), time = c(n[1] + n[2], n[3] + n[4]))
  D <- rep(c(0, 1), time = c(n[1] + n[3], n[2] + n[4]))
  
  write.table(X, paste0(out_dir, "X_", b, ".txt"), row.names = F, col.names = F)
  write.table(Y, paste0(out_dir, "Y_", b, ".txt"), row.names = F, col.names = F)
  write.table(D, paste0(out_dir, "D_", b, ".txt"), row.names = F, col.names = F)
}

