# Figure 2 
# Simulation study of a secondary phenotype using a case-control sample
# The tilted distribution is a mixture of two Gaussians

# load source functions 
source_code_dir <- "/Src/"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

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
beta[nonnull_beta] <- rnorm(nnonnull_beta, 0, 0.5)
beta_0 <- 0

gamma <- numeric(p)
nonnull_gamma <- sample(1:p, nnonnull_gamma, replace= F)
gamma[nonnull_gamma] <- rnorm(nnonnull_gamma, 0, 0.5)
gamma_y <- 2
gamma_0 <- 0

# parameters for knockoff selections
q <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
offset <- 1
stat_fun <- stat.glmnet_lambdasmax

# Store the FDR and power
B <- 500 # number of repetitions
nselect_ko <- matrix(0, nrow = B, ncol = length(q))
nfd_ko <- matrix(0, nrow = B, ncol = length(q))
nselect_tilted <- matrix(0, nrow = B, ncol = length(q))
nfd_tilted <- matrix(0, nrow = B, ncol = length(q))

for(b in 1:B){
  # sample covariates
  X0 <- sample_gaussian(N, R)
  # sample responses 
  Y0 <- sample_y_linear(X0, beta, beta_0)
  # compute case-control probability
  pcase <- pselect(X0, Y0, params  = list(model = "gaussian", gamma_0 = gamma_0, gamma_y = gamma_y, gamma = gamma))[,1]
  # sample cases and controls
  obs <- case_control_sample(X0, Y0, ncase, ncontrol, N, pcase)
  X <- obs$X; Y <- obs$Y
  
  # create standard knockoffs 
  Xtilde <- create.gaussian(X, mu = rep(0, p), Sigma = Sigma, method = "sdp")
  # knockoff selection
  selection_ko <- knockoff_selection(X, Xtilde, Y, stat_fun, q, offset, beta = beta)
  
  # create exact tilted knockoffs (time the run)
  start_time <- Sys.time()
  XtildeTilted <- exact_knockoff(X, Y, n = list(npos = obs$npos, nneg = N - obs$npos, ncase = ncase, ncontrol = ncontrol),
                                 params = list(Sigma = Sigma, mu = rep(0, p),
                                               gamma = gamma, gamma_0 = gamma_0, gamma_y = gamma_y))
  end_time <- Sys.time()
  cat(end_time - start_time, "\n")
  selection_tilted_ko <- knockoff_selection(X, XtildeTilted, Y, stat_fun, q, offset, beta = beta)
  
  
  # print fdp and power in each run 
  cat(b, ":\n ", "# ko selection = ", selection_ko$result[4, 1], ", # false discovery = ",  selection_ko$result[4, 2], "\n")
  cat("# tilted ko selection = ", selection_tilted_ko$result[4, 1], ", # false discovery = ",  selection_tilted_ko$result[4, 2], "\n")
  
  # store results
  nselect_ko[b, ] <- selection_ko$result[ ,1]
  nfd_ko[b, ] <- selection_ko$result[ ,2]
  nselect_tilted[b, ] <- selection_tilted_ko$result[ ,1]
  nfd_tilted[b, ] <- selection_tilted_ko$result[ ,2]
}

outdir <- "/scratch/users/qzhao1/paisa/simV2/result/gaussian/"

write.table(nselect_ko, file = paste0(outdir, "nselect_ko.txt"), row.names = F, col.names = F)
write.table(nfd_ko, file = paste0(outdir, "nfd_ko.txt"), row.names = F, col.names = F)
write.table(nselect_tilted, file = paste0(outdir, "nselect_tilted.txt"), row.names = F, col.names = F)
write.table(nfd_tilted, file = paste0(outdir, "nfd_tilted.txt"), row.names = F, col.names = F)






