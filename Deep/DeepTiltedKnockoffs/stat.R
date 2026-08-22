
# K(xi, yj) = exp(-gamma * ||xi - yj||^2)
gaussian_kernel <- function(x, y, gamma) {
  xx <- torch::torch_sum(x * x, dim = 2L, keepdim = TRUE)
  yy <- torch::torch_sum(y * y, dim = 2L, keepdim = TRUE)
  xy <- torch::torch_mm(x, y$t())
  d2 <- xx + yy$t() - 2.0 * xy
  torch::torch_exp(-gamma * d2)
}


# Squared MMD Loss with mixture of Gaussian kernels
mmd_loss <- function(x, y, alphas, W = NULL) {
  n_x <- x$shape[1]
  n_y <- y$shape[1]
  dev <- x$device
  
  # initialize 
  Kxx <- torch::torch_zeros(n_x, n_x, device = dev)
  Kxy <- torch::torch_zeros(n_x, n_y, device = dev) 
  Kyy <- torch::torch_zeros(n_y, n_y, device = dev)
  
  # add up kernels across all bandwidth parameters
  for (alpha in alphas) {
    g   <- 1.0 / (2.0 * alpha^2)
    Kxx <- Kxx + gaussian_kernel(x, x, g)
    Kxy <- Kxy + gaussian_kernel(x, y, g)
    Kyy <- Kyy + gaussian_kernel(y, y, g)
  }
  
  # apply optional weight matrix
  if (!is.null(W)) {
    Kxx <- Kxx * W
    Kxy <- Kxy * W
    Kyy <- Kyy * W
  }
  
  # Create masks to exclude diagonal elements (only for within-sample kernels)
  diag_mask_x <- 1.0 - torch::torch_eye(n_x, device = dev)
  diag_mask_y <- 1.0 - torch::torch_eye(n_y, device = dev)
  
  # MMD estimator 
  Kxx_off_diag <- torch::torch_sum(Kxx * diag_mask_x) / (n_x * (n_x - 1))
  Kxy_mean     <- torch::torch_sum(Kxy) / (n_x * n_y)
  Kyy_off_diag <- torch::torch_sum(Kyy * diag_mask_y) / (n_y * (n_y - 1))
  
  # Return mean(K_xx off-diag) - 2*mean(K_xy) + mean(K_yy off-diag)
  ll <- Kxx_off_diag - 2.0 * Kxy_mean + Kyy_off_diag
  
  mmd_final <- torch::torch_clamp(ll, min = 0.0)^2
  
  mmd_final
}


# Second Moment Loss (not used)
second_moment_loss <- function(x, y) {
  n  <- x$shape[1]
  Mx <- torch::torch_mm(x$t(), x) / n
  My <- torch::torch_mm(y$t(), y) / n
  torch::torch_mean((Mx - My)^2)
}
