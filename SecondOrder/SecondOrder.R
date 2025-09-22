# Section 3.1 Second order tilted knockoff
# Simulation with a selected sample (selection with a logistic model based on both X and Y)

# input index  
args <- commandArgs(trailingOnly = T)
ind <- as.numeric(args[1])

# directories
source_code_dir <- "/Src/"  #The directory where all source code files are saved.
params_dir <- "/param/SecondOrder/"
out_dir <- paste0("/result/second_order/")

# load source functions 
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# set parameters
N <- 5000 # population size
p <- 200 # no. variables

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

# read coefficients
beta <- scan(paste0(params_dir, "beta.txt"))
beta_0 <- 0

gamma_x <- scan(paste0(params_dir, "gamma.txt"))
gamma_y <- 2
gamma_0 <- -4

# parameters for knockoff selections
q <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
offset <- 1
stat_fun <- stat.glmnet_lambdasmax

# store relevant quantities
result <- matrix(NA, ncol = 4, nrow = length(q))
rownames(result) = q
colnames(result) = c("nselect_ko", "ndf_ko", "nselect_so", "nfd_so")

# sample X and Y from the population
X0 <- sample_gaussian(N, R = R)
Y0 <- sample_y_linear(X0, beta, beta_0)

# compute selection probability
probs <- pselect(X0, Y0, params  = list(model = "logistic", gamma_0 = gamma_0, gamma_y = gamma_y, gamma_x = matrix(gamma_x,nrow = 1)))[,1]

# selected sample
S <- rbinom(N, 1, probs)

X <- X0[which(S == 1), ]
Y <- Y0[which(S == 1)]
ntot <- nrow(X) # no obs.
cat("tot # obs = ", ntot, "\n") 

# create standard knockoffs 
Xtilde <- create.gaussian(X, mu = rep(0, p), Sigma = Sigma, method = "sdp")
# knockoff selection
selection_ko <- knockoff_selection(X, Xtilde, Y, stat_fun, q, offset, beta = beta)

# create tilted knockoffs
XtildeTilted <- matrix(NA, nrow = ntot, ncol = p)
for(i in 1:ntot){
  if(i == 1) start_time <- Sys.time()
  # estimate the mean and covariance of the tilted distribution
  cov_hat <- estimate_mucov_tilted(Y[i], sample_x_fun = function(n) sample_gaussian(n = n, R = R), 
                                   pselect_fun = function(x, y) pselect(x, y, params  = list(model = "logistic", gamma_0 = gamma_0, gamma_y = gamma_y, 
                                                                                          gamma_x = matrix(gamma_x,nrow = 1))),
                                   B = 50)
  # compute the minimum eigenvalue
  cat(min(eigen(cov_hat$sigma_hat[1, ,])$value), "\n")
  # sample knockoffs
  XtildeTilted[i, ] <- create.gaussian(X[i, ,drop = F], mu = cov_hat$mu_hat[1, ], Sigma = cov_hat$sigma_hat[1, ,], method = "sdp")
  # report time 
  if(i == 10){
    end_time <- Sys.time()
    cat("Time to generate 10 tilted knockoffs = ", end_time - start_time, "\n")
  }
  if(i %% 50 == 0) cat(i, ",")
}

# tilted knockoff selection
selection_tilted_ko <- knockoff_selection(X, XtildeTilted, Y, stat_fun, q, offset, beta = beta)


# store results
result[,1:2] <- selection_ko$result
result[,3:4] <- selection_tilted_ko$result

write.table(result, file = paste0(out_dir, "result", ind, ".txt"), row.names = T, col.names = T)

# print fdp and power in each run 
cat("# ko selection = ", selection_ko$result[4, 1], ", # false discovery = ",  selection_ko$result[4, 2], "\n")
cat("# tilted ko selection = ", selection_tilted_ko$result[4, 1], ", # false discovery = ",  selection_tilted_ko$result[4, 2], "\n")





