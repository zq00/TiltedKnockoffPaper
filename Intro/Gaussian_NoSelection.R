# Supplement Fig 3
# Benchmark: knockoff FDR when there's no selection 

# load source functions 
source_code_dir <- "/home/users/qzhao1/paisa/collider/SimV2/Src/"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# set parameters
N <- 4000 # population size
p <- 400 # no. variables
nnonnull_beta <- p * 0.1
nnonnull_gamma <- p * 0.2

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

# parameters for knockoff selections
q <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
offset <- 1
stat_fun <- stat.glmnet_lambdasmax

# Store the FDR and power
B <- 500 # number of repetitions
nselect_ko <- matrix(0, nrow = B, ncol = length(q))
nfd_ko <- matrix(0, nrow = B, ncol = length(q))

for(b in 1:B){
  # sample covariates
  X <- sample_gaussian(N, R)
  # sample responses 
  Y <- sample_y_linear(X, beta, beta_0)
  
  # create standard knockoffs 
  Xtilde <- create.gaussian(X, mu = rep(0, p), Sigma = Sigma, method = "sdp")
  # knockoff selection
  selection_ko <- knockoff_selection(X, Xtilde, Y, stat_fun, q, offset, beta = beta)
  
  # print fdp and power in each run 
  cat(b, ":\n ", "# ko selection = ", selection_ko$result[4, 1], ", # false discovery = ",  selection_ko$result[4, 2], "\n")
  
  # store results
  nselect_ko[b, ] <- selection_ko$result[ ,1]
  nfd_ko[b, ] <- selection_ko$result[ ,2]
}

outdir <- "/scratch/users/qzhao1/paisa/simV2/result/gaussian/"

# store FDR and power 
write.table(nselect_ko, file = paste0(outdir, "nselect_ko_NoSelection.txt"), row.names = F, col.names = F)
write.table(nfd_ko, file = paste0(outdir, "nfd_ko_NoSelection.txt"), row.names = F, col.names = F)

# store coefficients
write.table(beta, file = paste0(outdir, "beta_NoSelection.txt"), row.names = F, col.names = F)




