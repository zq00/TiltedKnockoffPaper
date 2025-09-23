# Supplement 3.2 Metro knockoff
library(knockoff)
library(SNPknock)
library(graph)
library(gRbase)
library(igraph) 

args <- commandArgs(trailingOnly = T)
aa <- as.numeric(args[1])

# directories
coef_dir <- "/param/metro/"
data_dir <- "/result/metro/data/"
out_dir <- "/result/metro/result/"

# 0 -- Source functions 
source_code_dir <- "/Src/"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# 0.1 model coefficients
beta <- scan(paste0(coef_dir, "beta.txt"))
beta0 <- 0

gamma_x <-  scan(paste0(coef_dir, "gamma.txt"))
gamma_y <- 2
gamma_0 <- -4
nonnull_gamma <- which(gamma_x != 0)

p <- length(beta)

# 0.2 coefficients of the MC
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


# 0.3 -- selection probability 
sample_x_fun <- function(n){
  sample_mc(n, pInit, Q)
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

# 0.4 -- Compute the log-likelihood 
# Input:
# x - Vector of length p. The observed covariates 
# y - Numeric. Outcome
# d - case control status
ComputeLogDensity <- function(x, y, d){
  logProbsX <- log(pInit[x[1] + 1]) + sum(log(Qsmall[cbind(x[1:(p-1)]+1, x[2:p]+1)]))
  
  pselect <- selection_fun(x, y)
  logProbS <- ifelse(d == 1, log(pselect[1]), log(pselect[2]))
  
  return(logProbsX + logProbS)
}

# 0.5 knockoff selections
stat_fun <- get("stat.lasso_lambdasmax_bin")
q <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
offset <- 1 # offset = 1 means the knockoff+ procedure


# 1 -- preparation
# 1.1 -- compute a node ordering 
V <- 1:p
edL <- vector("list", length=p)
names(edL) <- as.character(V)
for(i in 1:p){
  if(i == 1){
    edges <- c(2)
  }else if(i == p){
    edges <- c(p-1)
  }else{
    edges <- c(i-1, i+1)
  }
  if(i %in% nonnull_gamma){
    edges <- c(edges, nonnull_gamma[which(nonnull_gamma!=i)])
  }
  edges <- as.character(edges)
  edL[[as.character(i)]] <- list(edges = edges)
}

gR <- graphNEL(nodes=as.character(V), edgeL=edL)
junction_tree <- get_junction_tree(gR)
ordered_nodes <- get_node_order(junction_tree$junction_tree, junction_tree$cliques, p)

# 2 -- load samples
X <- as.matrix(read.table(paste0(data_dir, "X_", aa, ".txt")))
Y <- scan(paste0(data_dir, "Y_", aa, ".txt"))
D <- scan(paste0(data_dir, "D_", aa, ".txt"))

n <- nrow(X)

# 3 -- generate knockoffs
# 3.1 -- sample proposals
Xprop <- matrix(sample(c(0:(K-1)), n * p,  replace = T), n, p)
probX <- matrix(1, n, p) / K # equal probability to be in each state
probXProp <- matrix(1,n , p ) / K

# 3.2 -- sample knockoffs
BB <- 1000
XTildeMetro <- matrix(NA, BB, p)
myindices <- 1:BB

for(index in 1:BB){
  i <- myindices[index]
  cat("\n", i, ":\n")
  
  # If a variable is the same as its proposal, then it does not affect anyone
  ind <- which(X[i, ] - Xprop[i, ] == 0)
  affected_vars <- list() # which variable's acceptance probability is affected by Xj?
  for(j in 1:p){
    affected_vars[[j]] <- numeric(length = 0L)
  }
  for(j in 1:p){
    if(j %in% ind) next; # if x[i] and its proposals are equal, then it does not affect any variable, move on;
    new_nbd <- setdiff(ordered_nodes$nbd[[j]]$nbd, ind) # if a variable is equal to its proposal, then it cannot be affected either
    if(length(new_nbd) > 0){
      for(k in 1:length(new_nbd)){
        affected_vars[[new_nbd[k]]] <- c(affected_vars[[new_nbd[k]]], j)
      }
    }
  }
  
  #Store the computed acceptance probabilities
  computed_probs <- list()
  for(j in 1:p){
    computed_probs[[j]] <- list(nbd = NULL, values = NULL, acc_prob = NULL)
  }
  for(j in 1:p){
    if(j %in% ind) next; # the acceptance probability doesn't depend on anyone 
    
    new_vertex <- ordered_nodes$nbd[[j]]$vertex
    nnbd <- length(ordered_nodes$nbd[[j]]$nbd) 
    
    if(nnbd > 0){
      new_nbd <- setdiff(sort(ordered_nodes$nbd[[j]]$nbd), ind)
      nnbd <- length(new_nbd)
    }
    if(nnbd == 0){
      computed_probs[[new_vertex]]$nbd <- NULL
      computed_probs[[new_vertex]]$acc_prob <- NA
    }else{
      new_value <- matrix(NA, nrow = 2^nnbd, ncol = 2)
      
      if(nnbd == 1){
        newcomb <- c("0", "1")
      }else{
        newcomb <- apply(sapply(0:(2^nnbd - 1), function(x) as.integer(intToBits(x)) )[1:nnbd,], 2, function(t) paste0(t, collapse = ""))}
      new_value[,1] <- newcomb
      
      computed_probs[[new_vertex]]$nbd <- new_nbd
      computed_probs[[new_vertex]]$values <- new_value
    }
  }
  XTildeMetro[index, ] <- metro(x = X[i, ], xprop = Xprop[i, ],  probsX = probX[i, ], probsXProp = probXProp[i, ], y = Y[i], d = D[i])
}

# 4 -- store the knockoffs
write.table(XTildeMetro, paste0(data_dir, "xtildeMetro_", aa, ".txt"), row.names = F, col.names = F)

# compute the FDR 
resultMetro <- knockoff_selection(X, XTildeMetro, Y, stat_fun, q, offset = offset, beta = beta)
cat("\n FDP = ", resultMetro$result[5,2] / resultMetro$result[5,1], "\n")
print(resultMetro$result) 