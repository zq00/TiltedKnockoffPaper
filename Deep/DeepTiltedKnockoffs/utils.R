#' Default Hyperparameters for KnockoffMachine
#' @param p Integer. Number of features.
#' @return A named list of hyperparameters.
#' @export
default_params <- function(p) {
  p <- as.integer(p)
  list(
    
    # training parameter
    epochs        = 100L,
    epoch_length  = 100L,
    batch_size    = 256L,
    lr            = 0.01,
    lr_milestones = NULL,
    lr_decay      = 0.1,
    max_norm      = 1.0,
    
    # data type
    family        = "continuous",
    test_size     = 0.0,
    y_dim         = 1L,
    
    # network architecture 
    num_layers    = 6L,
    hidden_dim    = 10L * p,
    noise_dim     = p,
    noise_std     = 1.0,
    
    # weights for losses 
    lambda_swap        = 2,  # mmd swap loss = alpha 
    lambda_swap_b      = 0,  # second moment swap loss: not used 
    lambda_cov         = 1,  # covariance loss
    lambda_cov_diag    = 2.0,  # weight on correlation between Xj and its knockoffs: lambda_cov * lambda_cov_diag = gamma 
    lambda_cov_offdiag = 4,  # weight on cross-covariance: lambda_cov * lambda_cov_offdiag = beta * lambda_3
    lambda_marginal    = 4,  # covariance loss = lambda_2 * beta
    lambda_mean        = 1,  # mean loss = lambda_1 * beta 
    target_corr        = NULL,
    
    # MMD kernel bandwidth
    alphas        = c(1, 2, 4, 8, 16, 32, 64, 128),
    
    # evaluate loss in stratum 
    n_bins        = 40L,
    bin_type      = "quantile",
    min_bin_size  = 10L,
  )
}
