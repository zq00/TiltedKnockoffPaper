## 01_sample_coef.R
## Run this script ONCE before all simulations.
## It samples the coefficients and saves them to a shared file.

rm(list = ls())
source("/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/src_V2/00_config.R")

set.seed(0)   # fix seed so coefficients are reproducible

# ============================================================
# 1. Response model coefficients
# ============================================================
if (response == "mean") { # only has mean effect 
  p_mean <- p * 0.2
  
  non_null_mean <- sample(1:p, p_mean, replace = FALSE)
  beta_mean     <- rnorm(p_mean, 0, 0.4)
  
  non_null_loc  <- non_null_mean
  non_null_interaction <- NULL
  beta_interaction     <- NULL
  
} else if (response == "interaction") { # mean and interaction effects 
  p_interaction <- p * 0.05
  p_mean        <- p * 0.1
  
  non_null_mean <- sample(1:p, p_mean, replace = FALSE)
  non_null_interaction <- matrix(NA, nrow = p_interaction * 2, ncol = 2)
  non_null_interaction[, 1] <- sample(1:p, p_interaction * 2, replace = FALSE)
  non_null_interaction[, 2] <- sample(1:p, p_interaction * 2, replace = FALSE)
  non_null_interaction <- unique(non_null_interaction)
  non_null_interaction <- non_null_interaction[1:p_interaction, ]
  
  beta_mean        <- rnorm(p_mean, 0, 0.5)
  beta_interaction <- rnorm(p_interaction, 0, 0.25) # slightly weaker interaction effects
  
  non_null_loc <- unique(c(non_null_mean, non_null_interaction))
}

# ============================================================
# 2. Selection model coefficients
# ============================================================
if (selection == "interaction") {
  p_select_interaction_y <- p * 0.04
  p_select_interaction_x <- p * 0.08
  p_select_mean          <- p * 0.2
  
  non_null_select_mean <- sample(1:p, p_select_mean, replace = FALSE)
  non_null_select_interaction_x <- matrix(NA, nrow = p_select_interaction_x * 2, ncol = 2)
  non_null_select_interaction_x[, 1] <- sample(1:p, p_select_interaction_x * 2, replace = FALSE)
  non_null_select_interaction_x[, 2] <- sample(1:p, p_select_interaction_x * 2, replace = FALSE)
  non_null_select_interaction_x <- unique(non_null_select_interaction_x)
  non_null_select_interaction_x <- non_null_select_interaction_x[1:p_select_interaction_x, ]
  
  non_null_select_interaction_y <- sample(1:p, p_select_interaction_y, replace = FALSE)
  
  gamma_mean          <- rnorm(p_select_mean, 0, 0.5)
  gamma_interaction_x <- rnorm(nrow(non_null_select_interaction_x), 0, 0.25)
  gamma_interaction_y <- rnorm(length(non_null_select_interaction_y), 0, 0.25)
  
  non_null_gamma <- unique(c(non_null_select_mean,
                             non_null_select_interaction_x,
                             non_null_select_interaction_y))
  
} else if (selection == "mean") {
  p_select_mean          <- p * 0.4
  non_null_select_mean   <- sample(1:p, p_select_mean, replace = FALSE)
  gamma_mean             <- rnorm(p_select_mean, 0, 0.4)
  non_null_gamma         <- non_null_select_mean
  
  # placeholders so save does not fail
  non_null_select_interaction_x <- NULL
  non_null_select_interaction_y <- NULL
  gamma_interaction_x <- NULL
  gamma_interaction_y <- NULL
}

# ============================================================
# 3. Save
# ============================================================
if (!dir.exists(path_data)) dir.create(path_data, recursive = TRUE)

save(
  # response model
  non_null_mean, non_null_interaction, beta_mean, beta_interaction, non_null_loc,
  # selection model
  non_null_select_mean, non_null_select_interaction_x, non_null_select_interaction_y,
  gamma_mean, gamma_interaction_x, gamma_interaction_y, non_null_gamma,
  file = coef_file
)

cat("Coefficients saved to:", coef_file, "\n")
