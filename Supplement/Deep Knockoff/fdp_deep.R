# 3/1/2024
# Compute the FDP using the deep + metro knockoff 

# directories
coef_dir <- "/param/deep/"
data_dir <- "/result/deep/data/"
metro_ko_dir <- "/result/metro/ko/"
out_dir <- "/result/deep/result/"

# load functions
source_code_dir <- "/Src/"  #The directory where all source code files are saved.
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

# compute parameters for the tilted knockoff 
sample_x_fun <- function(n){
  sample_mc(n, pInit, Q)
}

selection_fun <- function(X, Y){
  pselect(X, Y, params = list(
    model = "logistic",
    gamma_0 = gamma_0,
    gamma_x = t(gamma_x), 
    gamma_y = gamma_y
  ))
}

vcov0 <- list()

vcov0[[1]] <- estimate_mucov_tilted(y = -1, sample_x_fun = sample_x_fun, pselect_fun = selection_fun, B = 50)
vcov0[[2]] <- estimate_mucov_tilted(y = 1, sample_x_fun = sample_x_fun, pselect_fun = selection_fun, B = 50)

vcov <- list()
vcov[[1]] <- list(mu = vcov0[[1]]$mu_hat[2, ], sigma = vcov0[[1]]$sigma_hat[2, ,])
vcov[[2]] <- list(mu = vcov0[[1]]$mu_hat[1, ], sigma = vcov0[[1]]$sigma_hat[1, ,])
vcov[[3]] <- list(mu = vcov0[[2]]$mu_hat[2, ], sigma = vcov0[[2]]$sigma_hat[2, ,])
vcov[[4]] <- list(mu = vcov0[[2]]$mu_hat[1, ], sigma = vcov0[[2]]$sigma_hat[1, ,])

# Store relevant quantities
B <- 500

nselect_ko <- matrix(0, nrow = B, ncol = length(q))
nfd_ko <- matrix(0, nrow = B, ncol = length(q))
nselect_so <- matrix(0, nrow = B, ncol = length(q))
nfd_so <- matrix(0, nrow = B, ncol = length(q))
nselect_deep <- matrix(0, nrow = B, ncol = length(q))
nfd_deep <- matrix(0, nrow = B, ncol = length(q))

names <- c("control_neg", "case_neg", "control_pos", "case_pos")
for(b in 1:B){
  cat(b, ":")
  # Load X, Y and XTilde
  X <- NULL
  XTildeDeep <- NULL
  XTildeKo <- NULL
  XTildeSo <- NULL
  XTildeMetro <- NULL
  
  n <- numeric(4)
  for(i in 1:4){
    xnew <- as.matrix(read.table(paste0(data_dir,names[i], "_", b, ".txt")))
    # deep knockoff
    xkdeepnew <- as.matrix(read.table(paste0(data_dir, "xk_", names[i], "_", b, ".txt")))
    # standard knockoff
    xknew <- knockoffDMC(xnew, pInit, Q)
    # SO knockoff
    xksonew <- create.gaussian(xnew, mu = vcov[[i]]$mu, Sigma = vcov[[i]]$sigma, method="sdp")
      
    # generate the SO knockoff 
    X <- rbind(X, xnew)
    XTildeKo <- rbind(XTildeKo, xknew)
    XTildeDeep <- rbind(XTildeDeep, xkdeepnew)
    XTildeSo <- rbind(XTildeSo, xksonew)
    
    
    n[i] <- nrow(xnew)
  }
  
  # generate the standard knockoff 
  Y <- rep(c( -1, 1), time = c(n[1] + n[2], n[3] + n[4]))
  # load the metro knockoff
  for(i in 1:10){
    xkmetronew <- as.matrix(read.table(paste0(metro_ko_dir, "xtildeMetro_", b, "_", i, ".txt")))
    XTildeMetro <- rbind(XTildeMetro, xkmetronew)
  }
  
  
  # compute result for each type of knockoffs
  resultKO <- knockoff_selection(X, XTildeKo, Y, stat_fun, q, offset = offset, beta = beta)
  resultSO <- knockoff_selection(X, XTildeSo, Y, stat_fun, q, offset = offset, beta = beta)
  resultDeep <- knockoff_selection(X, XTildeDeep, Y, stat_fun, q, offset = offset, beta = beta)
  resultMetro <- knockoff_selection(X, XTildeMetro, Y, stat_fun, q, offset = offset, beta = beta)
  
  # store the number of selections and false discoveries 
  nselect_ko[b, ] <- resultKO$result[,1]
  nfd_ko[b, ] <- resultKO$result[,2]
  
  nselect_so[b, ] <- resultSO$result[,1]
  nfd_so[b, ] <- resultSO$result[,2]
  
  nselect_deep[b, ] <- resultDeep$result[,1]
  nfd_deep[b, ] <- resultDeep$result[,2]
  
  nselect_metro[b, ] <- resultMetro$result[,1]
  nfd_metro[b, ] <- resultMetro$result[,2]
  
  cat("ko:", resultKO$result[5,2] / resultKO$result[5,1], "\n")
  cat("so:", resultSO$result[5,2] / resultSO$result[5,1], "\n")
  cat("deep:", resultDeep$result[5,2] / resultDeep$result[5,1], "\n")
  cat("metro:", resultMetro$result[5,2] / resultMetro$result[5,1], "\n")
}


# store results 
write.table(nselect_ko, file = paste0(out_dir, "nselect_ko.txt"), row.names = F, col.names = F)
write.table(nfd_ko, file = paste0(out_dir, "nfd_ko.txt"), row.names = F, col.names = F)
write.table(nselect_so, file = paste0(out_dir, "nselect_so.txt"), row.names = F, col.names = F)
write.table(nfd_so, file = paste0(out_dir, "nfd_so.txt"), row.names = F, col.names = F)
write.table(nselect_deep, file = paste0(out_dir, "nselect_deep.txt"), row.names = F, col.names = F)
write.table(nfd_deep, file = paste0(out_dir, "nfd_deep.txt"), row.names = F, col.names = F)
write.table(nselect_metro, file = paste0(out_dir, "nselect_metro.txt"), row.names = F, col.names = F)
write.table(nfd_metro, file = paste0(out_dir, "nfd_metro.txt"), row.names = F, col.names = F)





