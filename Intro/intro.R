# Simulation in Section 2.1
rm(list = ls())
# load source functions
source_code_dir <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/src/" # source code directory
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# function: BH
bh_thresh <- function(pval, q){
p <- length(pval)
sorted_pval <- sort(pval)
i0 <- max(sum(sorted_pval < q * (1:p) / p), 1)
tau <- q * i0 / p # BH-q threshold
return(tau)
}

# set parameters
N <- 40000 # population size
p <- 400 # no. variables
nnonnull_beta <- p * 0.1
nnonnull_gamma <- p * 0.2
ncase <- 2000
ncontrol <- 2000

# set the covariance matrix
rho <- 0.5 # each block is a toeplitz matrix Sigma_(i,j) = rho^|i - j|
bsize <- 10 # block size
bdiag <- toeplitz(rho^(0:(bsize - 1)))
dd <- list()
for(i in 1:(p %/% bsize)){
dd[[i]] <- bdiag
}
if(p%%bsize != 0){
dd[[i + 1]] <- bdiag[1:(p%%bsize), 1:(p%%bsize)]
}
Sigma <- as.matrix(bdiag(dd)) # covariance matrix
R <- chol(Sigma)

# sample coefficients
beta <- numeric(p)
nonnull_beta <- sample(1:p, nnonnull_beta, replace= F)
beta[nonnull_beta] <- runif(nnonnull_beta, 0.05, 0.15)
beta[nonnull_beta] <- beta[nonnull_beta] * sample(c(1, -1), nnonnull_beta, replace = T)

# sample coefficients
beta <- numeric(p)
nonnull_beta <- sample(1:p, nnonnull_beta, replace= F)
beta[nonnull_beta] <- runif(nnonnull_beta, 0.05, 0.15)
beta[nonnull_beta] <- beta[nonnull_beta] * sample(c(1, -1), nnonnull_beta, replace = T)
beta_0 <- 0
gamma <- numeric(p)
nonnull_gamma <- sample(1:p, nnonnull_gamma, replace= F)
gamma[nonnull_gamma] <- runif(nnonnull_gamma, 0.1, 0.3)
gamma[nonnull_gamma] <- gamma[nonnull_gamma] * sample(c(1, -1), nnonnull_gamma, replace = T)
gamma_y <- 2
gamma_0 <- 0

# parameters for knockoff selections
q <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
offset <- 1
stat_fun <- stat.glmnet_lambdasmax

# Store the FDR and power
B <- 200 # number of repetitions
nselect_ko <- matrix(0, nrow = B, ncol = length(q))
nfd_ko <- matrix(0, nrow = B, ncol = length(q))
nselect_tilted <- matrix(0, nrow = B, ncol = length(q))
nfd_tilted <- matrix(0, nrow = B, ncol = length(q))
# ipw
nselect_ipw <- matrix(0, nrow = B, ncol = length(q))
nfd_ipw <- matrix(0, nrow = B, ncol = length(q))
# ko + ipw
nselect_ko_ipw <- matrix(0, nrow = B, ncol = length(q))
nfd_ko_ipw <- matrix(0, nrow = B, ncol = length(q))
# case/control ko
nselect_ko_case_control <- matrix(0, nrow = B, ncol = length(q))
nfd_ko_case_control <- matrix(0, nrow = B, ncol = length(q))

for(b in 1:B){
  # sample covariates
  X0 <- sample_gaussian(N, R = R)
  # sample responses 
  Y0 <- sample_y_linear(X0, beta, beta_0)
  # compute case-control probability
  pcase <- pselect(X0, Y0, 
                   params  = list(model = "gaussian", gamma_0 = gamma_0, gamma_y = gamma_y, 
                                  gamma_x = matrix(gamma, nrow = 1)))[,1]
  
  # sample cases and controls
  obs <- case_control_sample(X0, Y0, ncase, ncontrol, N, pcase)
  X <- obs$X; Y <- obs$Y
  
  # ------
  # 1. Inverse probability weighting (ipw)
  # ------
  pi0 <- obs$npos / N
  
  w0 <- pselect(X, Y, 
                params  = list(model = "gaussian", gamma_0 = gamma_0, gamma_y = gamma_y, 
                               gamma_x = matrix(gamma, nrow = 1)))[,1]
  w <- c(pi0 / w0[1:ncase], (1- pi0) / (1-w0[(ncase + 1):(ncase + ncontrol)]))
  
  fit <- lm(Y ~ X, weights = w)
  
  p_val <- summary(fit)$coef[-1,4]
  
  # use bh to select variables at levels of q 
  for(k in 1:length(q)){
    thresh <- bh_thresh(p_val, q[k])
    selected <- which(p_val <= thresh)
    nselect_ipw[b, k] <- length(selected)
    nfd_ipw[b, k] <- sum(!selected %in% nonnull_beta)
  }
  
  # ------
  # 2. Standard knockoffs
  # ------
  
  # create standard knockoffs 
  Xtilde <- create.gaussian(X, mu = rep(0, p), Sigma = Sigma, method = "sdp")
  # knockoff selection
  selection_ko <- knockoff_selection(X, Xtilde, Y, stat_fun, q, offset, beta = beta)
  
  nselect_ko[b, ] <- selection_ko$result[ ,1]
  nfd_ko[b, ] <- selection_ko$result[ ,2]
  
  # ------
  # 3. Tilted knockoffs
  # ------
  
  # create exact tilted knockoffs (time the run)
  start_time <- Sys.time()
  XtildeTilted <- exact_knockoff(X, Y, n = list(npos = obs$npos, nneg = N - obs$npos, ncase = ncase, ncontrol = ncontrol),
                                 params = list(Sigma = Sigma, mu = rep(0, p),
                                               gamma = gamma, gamma_0 = gamma_0, gamma_y = gamma_y))
  end_time <- Sys.time()
  cat(end_time - start_time, "\n")
  selection_tilted_ko <- knockoff_selection(X, XtildeTilted, Y, stat_fun, q, offset, beta = beta)
  
  # store results
  nselect_tilted[b, ] <- selection_tilted_ko$result[ ,1]
  nfd_tilted[b, ] <- selection_tilted_ko$result[ ,2]
  
  # print fdp and power in each run 
  cat(b, ":\n ", "# ko true discovery = ", selection_ko$result[4, 1] - selection_ko$result[4, 2], 
      ", FDP = ",  selection_ko$result[4, 2] / max(selection_ko$result[4, 1], 1), "\n")
  cat("# Tilted knockoff true discovery = ", selection_tilted_ko$result[4, 1] - selection_tilted_ko$result[4, 2], 
      ", FDP = ",  selection_tilted_ko$result[4, 2] / max(selection_tilted_ko$result[4, 1], 1), "\n")  
  
  # ------
  # 4. Standard knockoffs with IPW
  # ------
  
  # knockoff selection
  selection_ko_ipw <- knockoff_selection(X, Xtilde, Y, stat_fun, q, offset, beta = beta, weights = w)
  
  nselect_ko_ipw[b, ] <- selection_ko_ipw$result[ ,1]
  nfd_ko_ipw[b, ] <- selection_ko_ipw$result[ ,2]
  
  # print fdp and power in each run 
  cat("# ko + IPW true discovery  =", selection_ko_ipw$result[4, 1] - selection_ko_ipw$result[4, 2], 
      ", FDP = ",  selection_ko_ipw$result[4, 2] / max(selection_ko_ipw$result[4, 1], 1), "\n")  
  
  
  # ------
  # 5. Case-control specific knockoffs
  # ------

  # create standard knockoffs 
  X_case <- X[which(obs$D == 1), ]
  X_control <- X[which(obs$D == 0), ]
  Y_case <- Y[which(obs$D == 1)]
  Y_control <- Y[which(obs$D == 0)]
  
  Xtilde_case <- create.second_order(X = X_case, method = "sdp")
  Xtilde_control <- create.second_order(X =X_control, method = "sdp")
  # knockoff selection
  selection_ko_case_control <- knockoff_selection(rbind(X_case, X_control), 
                                     rbind(Xtilde_case, Xtilde_control), 
                                     c(Y_case, Y_control), 
                                     stat_fun, q, offset, beta = beta)
  
  nselect_ko_case_control[b, ] <- selection_ko_case_control$result[ ,1]
  nfd_ko_case_control[b, ] <- selection_ko_case_control$result[ ,2]
  
  # print fdp and power in each run 
  cat("# case-control knockoff true discovery  =", selection_ko_case_control$result[4, 1] - selection_ko_case_control$result[4, 2], 
      ", FDP = ",  selection_ko_case_control$result[4, 2] / max(selection_ko_case_control$result[4, 1], 1), "\n")  
  
  
  }

# ------
# 5. Store results
# ------

outdir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/intro/"

save.image(file = paste0(outdir, "intro_V2.RData"))




