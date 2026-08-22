# Shared configuration for the deep tilted knockoff simulation.
# Source this file at the top of every other script.

# ============================================================
# Paths 
# ============================================================
path_src   <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/src_V2/" # source code for simulation
path_data  <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/Deep/data" # dir to store data

# ============================================================
# Simulation parameters
# ============================================================
n_sim      <- 50L        # total number of simulation replicates

p          <- 100L       # number of variables
N          <- 150000L     # population size

n_case     <- 2000L      # case-control sample sizes
n_control  <- 2000L

n_train    <- 100000L    # training data size for deep knockoff

# ============================================================
# Model settings
# ============================================================

setting    <- "case_control"
covariate  <- "gaussian"
response   <- "mean"
selection  <- "mean"

# ============================================================
# Helper functions
# ============================================================

logit     <- function(t) log(t / (1 - t))
logistic  <- function(t) 1 / (1 + exp(-t))

# ============================================================
# Build directory for b-th simulation replicate
# ============================================================

sim_dir <- function(b) {
  d <- file.path(path_data, sprintf("sim_%03d", b))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

# Path to the shared coefficient file (same for all sims)
coef_file <- file.path(path_data, "coef.RData")
