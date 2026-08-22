## 03_train.R
## For simulation replicate b, train the deep tilted knockoff machine.
## Trains BOTH the case and control machines.
##
## Usage: Rscript 03_train.R <b>
##   where b = 1, 2, ..., n_sim
##
## Alternatively, to train only case or control (e.g., on separate GPUs):
##   Rscript 03_train.R <b> case
##   Rscript 03_train.R <b> control

rm(list = ls())
source("/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/src_V2/00_config.R")

# ============================================================
# 0. Parse command-line arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 03_train.R <b> [case|control|both]")
b    <- as.integer(args[1])
mode <- if (length(args) >= 2) args[2] else "both"   # "case", "control", or "both"
method <- args[3]
cat("=== Training replicate:", b, "| mode:", mode, "===\n")

# ============================================================
# 1. Load functions and training data from Step 02
# ============================================================
out_dir <- sim_dir(b)
load(file.path(out_dir, "sample_data.RData"))

source("~/Research/TiltedKnockoff/Deep/src_V2/load_files.R")
# ============================================================
# 3. Set machine parameters
# ============================================================
my_pars <- list(
  epochs        = 50L,
  epoch_length  = 25L,
  batch_size    = 4096L,
  lr            = 0.01,
  lr_milestones = c(30,40),
  lr_decay      = 0.1,
  hidden_dim    = 10L * p,
  alphas        = c(1, 2, 4, 8, 16, 32, 64, 128),
  family        = "continuous",
  test_size     = 0.0,
  y_dim         = 1L
)



# ============================================================
# 4. Train CASE machine
# ============================================================
if (mode %in% c("both", "case")) {
  cat("--- Training CASE machine ---\n")
  
  Y_train_case <- get(paste0("training_data_", method, "_case"))$y_train
  X_train_case <- get(paste0("training_data_", method, "_case"))$x_train
  
  X_case_mean <- colMeans(X_train_case)
  X_case_sd   <- apply(X_train_case, 2, sd)
  Y_case_mean <- mean(Y_train_case)
  Y_case_sd   <- sd(Y_train_case)
  
  X_case_s <- scale(X_train_case, center = X_case_mean, scale = X_case_sd)
  Y_case_s <- (Y_train_case - Y_case_mean) / Y_case_sd
  
  
  # reload functions (ensures clean state)
  source(paste0(path_deep_tilted, "filmnet.R"))
  source(paste0(path_deep_tilted, "machine_V5.R"))
  source(paste0(path_deep_tilted, "stat_V3.R"))
  source(paste0(path_deep_tilted, "utils_V3.R"))
  source(paste0(path_deep_tilted, "mmd_eval_V2.R"))
  
  machine <- KnockoffMachine(p = p, pars = my_pars)
  machine_case <- train(machine, X_case_s, Y_case_s, verbose = TRUE)
  
  # save
  torch::torch_save(machine_case$net, 
                    file.path(out_dir, paste0(method, "_machine_case_net.pt")))
  save(machine_case, X_case_mean, X_case_sd, Y_case_mean, Y_case_sd,
       file = file.path(out_dir, paste0(method, "_machine_case.RData")))
  cat("Case machine saved.\n")
}

# ============================================================
# 5. Train CONTROL machine
# ============================================================
if (mode %in% c("both", "control")) {
  cat("--- Training CONTROL machine ---\n")
  
  Y_train_control <- training_data_logistic_control$y_train
  X_train_control <- training_data_logistic_control$x_train
  
  X_control_mean <- colMeans(X_train_control)
  X_control_sd   <- apply(X_train_control, 2, sd)
  Y_control_mean <- mean(Y_train_control)
  Y_control_sd   <- sd(Y_train_control)
  
  X_control_s <- scale(X_train_control, center = X_control_mean, scale = X_control_sd)
  Y_control_s <- (Y_train_control - Y_control_mean) / Y_control_sd
  
  
  # reload functions
  source(paste0(path_deep_tilted, "filmnet.R"))
  source(paste0(path_deep_tilted, "machine_V5.R"))
  source(paste0(path_deep_tilted, "stat_V3.R"))
  source(paste0(path_deep_tilted, "utils_V3.R"))
  source(paste0(path_deep_tilted, "mmd_eval_V2.R"))
  
  machine <- KnockoffMachine(p = p, pars = my_pars)
  machine_control <- train(machine, X_control_s, Y_control_s, verbose = TRUE)
  
  # save
  torch::torch_save(machine_control$net, 
                    file.path(out_dir, paste0(method, "_machine_control_net.pt")))
  save(machine_control, X_control_mean, X_control_sd, Y_control_mean, Y_control_sd,
       file = file.path(out_dir,  paste0(method, "_machine_control.RData")))
  cat("Control machine saved.\n")
}

cat("=== Training complete for replicate", b, ", method ",method ,"===\n")
