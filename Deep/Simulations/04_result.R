## 04_result.R
## For simulation replicate b:
##   - load trained knockoff machines
##   - generate knockoff copies of observed data
##   - compute FDR for the deep tilted knockoff method
##
## Usage: Rscript 04_result.R <b>
##   where b = 1, 2, ..., n_sim

rm(list = ls())
source("/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/src_V2/00_config.R")

# ============================================================
# 0. Parse command-line argument
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 04_result.R <b>")
b <- as.integer(args[1])
method <- args[2]
cat("=== Training replicate:", b, "| method:", method, "===\n")

cat("=== Results for replicate:", b, "===\n")

# ============================================================
# 1. Load packages and functions
# ============================================================
source("~/Research/TiltedKnockoff/Deep/src_V2/load_files.R")

# ============================================================
# 2. Load data and trained machines
# ============================================================
out_dir <- sim_dir(b)

load(file.path(out_dir, "sample_data.RData"))
load(file.path(out_dir, paste0(method, "_machine_case.RData")))
load(file.path(out_dir, paste0(method, "_machine_control.RData")))

# Rebuild machine objects and load network weights
machine_case <- KnockoffMachine(p = machine_case$p, pars = machine_case$pars)
machine_case$net <- torch::torch_load(file.path(out_dir, paste0(method, "_machine_case_net.pt")))

machine_control <- KnockoffMachine(p = machine_control$p, pars = machine_control$pars)
machine_control$net <- torch::torch_load(file.path(out_dir, paste0(method, "_machine_control_net.pt")))

# ============================================================
# 3. Generate deep tilted knockoffs
# ============================================================

# --- Case ---
id_case     <- which(status == "case")
Y_obs_case  <- Y_obs[id_case]
X_obs_case  <- X_obs[id_case, ]

X_case_mean <- colMeans(X_obs_case)
X_case_sd <- apply(X_obs_case, 2, sd)
Y_case_mean <- mean(Y_obs_case)
Y_case_sd <- sd(Y_obs_case)

X_obs_case_s <- scale(X_obs_case, center = X_case_mean, scale = X_case_sd)
Y_obs_case_s <- (Y_obs_case - Y_case_mean) / Y_case_sd

Xk_case_s <- generate(machine_case, X_obs_case_s, Y_obs_case_s)
Xk_case   <- t(t(Xk_case_s) * X_case_sd + X_case_mean)

# --- Control ---
id_control     <- which(status == "control")
Y_obs_control  <- Y_obs[id_control]
X_obs_control  <- X_obs[id_control, ]

X_control_mean <- colMeans(X_obs_control)
X_control_sd <- apply(X_obs_control, 2, sd)
Y_control_mean <- mean(Y_obs_control)
Y_control_sd <- sd(Y_obs_control)

X_obs_control_s <- scale(X_obs_control, center = X_control_mean, scale = X_control_sd)
Y_obs_control_s <- (Y_obs_control - Y_control_mean) / Y_control_sd

Xk_control_s <- generate(machine_control, X_obs_control_s, Y_obs_control_s)
Xk_control   <- t(t(Xk_control_s) * X_control_sd + X_control_mean)

# ============================================================
# 4. Knockoff selection settings
# ============================================================

stat_fun <- stat.lasso_lambdasmax

# true beta indicator
beta <- numeric(p)
beta[non_null_loc] <- 1

# target FDR levels
q <- seq(0.05, 0.5, by = 0.05)

# ============================================================
# 5. Deep tilted knockoff selections
# ============================================================

# combine case + control
X_all         <- rbind(X_obs_control, X_obs_case)
Xk_tilted_all <- rbind(Xk_control, Xk_case)
Y_all         <- c(Y_obs_control, Y_obs_case)

selection_tilted_ko <- knockoff_selection(X_all, Xk_tilted_all, Y_all,
                                          stat_fun, q, offset = 1, beta = beta)

result_tilted <- tibble(
  sim     = b,
  q       = q,
  nselect = selection_tilted_ko$result[, 1],
  nfd     = selection_tilted_ko$result[, 2]
) %>%
  mutate(fdp = round(nfd / pmax(nselect, 1), 4))

cat("Deep Tilted Knockoff FDP:\n")
print(result_tilted)

# ============================================================
# 6. Standard knockoff 
# ============================================================
Xk_standard <- create.gaussian(X_all, mu = rep(0, p), Sigma = Sigma, method = "asdp")
selection_ko <- knockoff_selection(X_all, Xk_standard, Y_all,
                                   stat_fun, q, offset = 1, beta = beta)

result_standard <- tibble(
  sim     = b,
  q       = q,
  nselect = selection_ko$result[, 1],
  nfd     = selection_ko$result[, 2]
) %>%
  mutate(fdp = round(nfd / pmax(nselect, 1), 4))

cat("\nStandard Knockoff FDP:\n")
print(result_standard)

# ============================================================
# 7. Second order knockoff (binned)
# ============================================================

n_bins <- 25

# estimate mean and covariances within each bin for case/control separately
mean_list_case <- list()
cov_list_case <- list()

breaks_case <- quantile(Y_obs_case, probs = seq(0, 1, length.out = n_bins + 1))

case_bin_ids <- cut(get(paste0("training_data_", method, "_case"))$y_train,
                    breaks         = breaks_case,
                    include.lowest = TRUE,
                    labels         = FALSE)

for(i in 1:(length(breaks_case) - 1)){
  id <- which(case_bin_ids == i)
  cat(i, ":", length(id), "\n")
  x_small <- get(paste0("training_data_", method, "_case"))$x_train[id, ]
  mean_list_case[[i]] <- colMeans(x_small)
  cov_list_case[[i]] <- cov(x_small)
}

mean_list_control <- list()
cov_list_control <- list()

breaks_control <- quantile(Y_obs_control, probs = seq(0, 1, length.out = n_bins + 1))
control_bin_ids <- cut(get(paste0("training_data_", method, "_control"))$y_train,
                       breaks         = breaks_control,
                       include.lowest = TRUE,
                       labels         = FALSE)

for(i in 1:(length(breaks_control) - 1)){
  id <- which(control_bin_ids == i)
  cat(i, ":", length(id), "\n")
  x_small <- get(paste0("training_data_", method, "_control"))$x_train[id, ]
  mean_list_control[[i]] <- colMeans(x_small)
  cov_list_control[[i]] <- cov(x_small)
}

# sample binned second order tilted knockoffs
Xk_so_case <- matrix(NA, nrow = nrow(X_obs_case), ncol = p)
Xk_so_control <- matrix(NA, nrow = nrow(X_obs_control), ncol = p)

case_bin_ids_obs <- cut(Y_obs_case,
                        breaks         = breaks_case,
                        include.lowest = TRUE,
                        labels         = FALSE)

for(i in 1:(length(breaks_case) - 1)){
  
  id <- which(case_bin_ids_obs == i)
  
  x_small <- X_obs_case[id, ]
  
  Xk_so_case[id,] <- create.gaussian(x_small, mu = mean_list_case[[i]], Sigma = cov_list_case[[i]], method = "asdp")
}


control_bin_ids_obs <- cut(Y_obs_control,
                           breaks         = breaks_control,
                           include.lowest = TRUE,
                           labels         = FALSE)

for(i in 1:(length(breaks_control) - 1)){
  id <- which(control_bin_ids_obs == i)
  x_small <- X_obs_control[id, ]
  
  Xk_so_control[id,] <- create.gaussian(x_small, mu = mean_list_control[[i]], Sigma = cov_list_control[[i]], method = "asdp")
}

# compute selections
X <- rbind(X_obs_control,X_obs_case)
Xk_tilted_so <- rbind(Xk_so_control, Xk_so_case)
Y <- c(Y_obs_control, Y_obs_case)

selection_tilted_ko_so <- knockoff_selection(X, Xk_tilted_so, Y, stat_fun, q, offset = 1, beta = beta)

result_so <- tibble(
  sim     = b,
  q       = q,
  nselect = selection_tilted_ko_so$result[, 1],
  nfd     = selection_tilted_ko_so$result[, 2]
) %>%
  mutate(fdp = round(nfd / pmax(nselect, 1), 4))

cat("\n Second order Knockoff FDP:\n")
print(result_so)


# ============================================================
# 8. Save results
# ============================================================
save(result_tilted, result_standard, result_so, 
     file = file.path(out_dir, "results.RData"))

cat("\nResults saved to:", file.path(out_dir, "results.RData"), "\n")
