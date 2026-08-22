# Train a tilted knockoff machine 

train    <- function(object, ...) UseMethod("train")

generate <- function(object, ...) UseMethod("generate")

# Create a new tilted knockoff machine 
KnockoffMachine <- function(p, pars = list()) {
  p      <- as.integer(p)
  params <- utils::modifyList(default_params(p), pars)
  
  params$y_dim        <- as.integer(params$y_dim)
  params$n_bins       <- as.integer(params$n_bins)
  params$min_bin_size <- as.integer(params$min_bin_size)
  params$epoch_length <- as.integer(params$epoch_length)
  params$batch_size   <- as.integer(params$batch_size)
  params$num_layers   <- as.integer(params$num_layers)
  params$hidden_dim   <- as.integer(params$hidden_dim)
  params$noise_dim    <- as.integer(params$noise_dim)
  
  params$lr              <- as.numeric(params$lr)
  params$lr_decay        <- as.numeric(params$lr_decay)
  params$noise_std       <- as.numeric(params$noise_std)
  params$lambda_swap     <- as.numeric(params$lambda_swap)
  params$lambda_swap_b   <- as.numeric(params$lambda_swap_b)
  params$lambda_cov      <- as.numeric(params$lambda_cov)
  params$lambda_marginal <- as.numeric(params$lambda_marginal)
  params$lambda_mean     <- as.numeric(params$lambda_mean)
  params$lambda_cov_diag    <- as.numeric(
    if (is.null(params$lambda_cov_diag))    1.0
    else params$lambda_cov_diag
  )
  params$lambda_cov_offdiag <- as.numeric(
    if (is.null(params$lambda_cov_offdiag)) 1.0
    else params$lambda_cov_offdiag
  )
  if (!is.null(params$max_norm)) {
    params$max_norm <- as.numeric(params$max_norm)
  }
 
  # target_corr_by_bin is target correlation in each bin
  # (1) A matrix of shape (n_bins, p) - rows are bins, columns are features
  # (2) A list of length n_bins, where each element is a numeric vector of length p
  if (!is.null(params$target_corr_by_bin)) {
    
    if (is.matrix(params$target_corr_by_bin)) {
      n_bins_spec <- nrow(params$target_corr_by_bin)
      p_spec      <- ncol(params$target_corr_by_bin)
      if (p_spec != p) {
        stop(sprintf(
          "target_corr_by_bin matrix must have %d columns (features). Got %d.",
          p, p_spec
        ))
      }
      params$target_corr_by_bin <- as.matrix(params$target_corr_by_bin)
      
    } else if (is.list(params$target_corr_by_bin)) {
      n_bins_spec <- length(params$target_corr_by_bin)
      for (i in seq_len(n_bins_spec)) {
        if (!is.numeric(params$target_corr_by_bin[[i]]) || 
            length(params$target_corr_by_bin[[i]]) != p) {
          stop(sprintf(
            "Each element of target_corr_by_bin must be a numeric vector of length %d. Element %d has length %d.",
            p, i, length(params$target_corr_by_bin[[i]])
          ))
        }
      }
      params$target_corr_by_bin <- lapply(params$target_corr_by_bin, as.numeric)
      
    } else {
      stop("target_corr_by_bin must be either a matrix (n_bins x p) or a list of length n_bins.")
    }
  }
  
  if (!is.null(params$lr_milestones)) {
    params$lr_milestones <- as.integer(params$lr_milestones)
  }
  
  if (any(params$alphas <= 0)) {
    stop("All values in alphas must be positive.")
  }
  params$alphas <- as.numeric(params$alphas)
  
  if (!params$family %in% c("continuous", "binary")) {
    stop(sprintf(
      "family must be 'continuous' or 'binary'. Got '%s'.",
      params$family
    ))
  }
  
  if (params$test_size < 0 || params$test_size >= 1) {
    stop(sprintf(
      "test_size must be in [0, 1). Got %.3f.",
      params$test_size
    ))
  }
  
  lambda_names <- c("lambda_swap", "lambda_swap_b", "lambda_cov",
                    "lambda_cov_diag", "lambda_cov_offdiag",
                    "lambda_marginal", "lambda_mean")
  for (nm in lambda_names) {
    if (params[[nm]] < 0) {
      stop(sprintf("'%s' must be non-negative. Got %.3f.", nm, params[[nm]]))
    }
  }
  
  # Build network
  net <- knockoff_net(
    p          = p,
    noise_dim  = params$noise_dim,
    hidden_dim = params$hidden_dim,
    num_layers = params$num_layers,
    y_dim      = params$y_dim
  )
  
  structure(
    list(p = p, pars = params, net = net),
    class = "KnockoffMachine"
  )
}


print.KnockoffMachine <- function(x, ...) {
  cat("── KnockoffMachine (S3) ──────────────────────\n")
  
  # ── Data ───────────────────────────────────────────────────────────
  cat("  [Data]\n")
  cat(sprintf("    Features (p)       : %d\n",   x$p))
  cat(sprintf("    Family             : %s\n",   x$pars$family))
  cat(sprintf("    Y dimension        : %d\n",   x$pars$y_dim))
  
  # ── Training ───────────────────────────────────────────────────────
  cat("  [Training]\n")
  cat(sprintf("    Epochs             : %d\n",   x$pars$epochs))
  cat(sprintf("    Epoch length       : %d\n",   x$pars$epoch_length))
  cat(sprintf("    Batch size         : %d\n",   x$pars$batch_size))
  cat(sprintf("    Learning rate      : %.5f\n", x$pars$lr))
  cat(sprintf("    LR milestones      : %s\n",
              ifelse(is.null(x$pars$lr_milestones), "none",
                     paste(x$pars$lr_milestones, collapse = ", "))))
  cat(sprintf("    LR decay           : %.3f\n", x$pars$lr_decay))
  cat(sprintf("    Max norm           : %s\n",
              ifelse(is.null(x$pars$max_norm), "none (no clipping)",
                     as.character(x$pars$max_norm))))
  cat(sprintf("    Test size          : %.2f\n", x$pars$test_size))
  
  # ── Network architecture ───────────────────────────────────────────
  cat("  [Network]\n")
  cat(sprintf("    Hidden dim         : %d\n",   x$pars$hidden_dim))
  cat(sprintf("    Num layers         : %d\n",   x$pars$num_layers))
  cat(sprintf("    Noise dim          : %d\n",   x$pars$noise_dim))
  cat(sprintf("    Noise std          : %.3f\n", x$pars$noise_std))
  
  # ── Loss weights ───────────────────────────────────────────────────
  cat("  [Loss Weights]\n")
  cat(sprintf("    lambda_swap        : %.3f  (MMD swap loss)\n",
              x$pars$lambda_swap))
  cat(sprintf("    lambda_swap_b      : %.3f  (second moment swap loss)\n",
              x$pars$lambda_swap_b))
  cat(sprintf("    lambda_cov         : %.3f  (covariance loss — outer weight)\n",
              x$pars$lambda_cov))
  cat(sprintf("    lambda_cov_diag    : %.3f  (cov loss — diagonal weight)\n",
              ifelse(is.null(x$pars$lambda_cov_diag),    1.0,
                     x$pars$lambda_cov_diag)))
  cat(sprintf("    lambda_cov_offdiag : %.3f  (cov loss — off-diagonal weight)\n",
              ifelse(is.null(x$pars$lambda_cov_offdiag), 1.0,
                     x$pars$lambda_cov_offdiag)))
  cat(sprintf("    lambda_marginal    : %.3f  (marginal second moment loss)\n",
              ifelse(is.null(x$pars$lambda_marginal),    1.0,
                     x$pars$lambda_marginal)))
  cat(sprintf("    lambda_mean        : %.3f  (mean matching loss)\n",
              ifelse(is.null(x$pars$lambda_mean),        1.0,
                     x$pars$lambda_mean)))
  
  cat("  [Covariance Loss Details]\n")
  
  # bin-specific target correlations 
  if (!is.null(x$pars$target_corr_by_bin)) {
    if (is.matrix(x$pars$target_corr_by_bin)) {
      n_bins_display <- nrow(x$pars$target_corr_by_bin)
    } else {
      n_bins_display <- length(x$pars$target_corr_by_bin)
    }
    cat(sprintf("    Bin-specific targets: YES (%d bins x %d features)\n", 
                n_bins_display, x$p))
  } else {
    cat(sprintf("    Bin-specific targets: NO\n"))
  }
  
  # ── MMD kernel ────────────────────────────────────────────────────
  cat("  [MMD Kernel]\n")
  cat(sprintf("    Alphas             : %s\n",
              paste(x$pars$alphas, collapse = ", ")))
  
  # ── Evaluation ────────────────────────────────────────────────────
  cat("  [Evaluation]\n")
  cat(sprintf("    N bins             : %d\n",   x$pars$n_bins))
  cat(sprintf("    Bin type           : %s\n",   x$pars$bin_type))
  cat(sprintf("    Min bin size       : %d\n",   x$pars$min_bin_size))
  
  cat("──────────────────────────────────────────────\n")
  invisible(x)
}

summary.KnockoffMachine <- function(object, ...) print.KnockoffMachine(object)

.km_loss_stratified_all <- function(object, X, Xk, Y) {
  
  n_bins       <- object$pars$n_bins
  bin_type     <- object$pars$bin_type
  min_bin_size <- object$pars$min_bin_size
  
  Y_vec <- .to_y_vec(Y)
  
  if (bin_type == "quantile") {
    breaks <- quantile(Y_vec, probs = seq(0, 1, length.out = n_bins + 1))
  } else if (bin_type == "equal_width") {
    breaks <- seq(min(Y_vec), max(Y_vec), length.out = n_bins + 1)
  } else {
    stop(sprintf("Unknown bin_type '%s'.", bin_type))
  }
  
  breaks        <- unique(breaks)
  n_bins_actual <- length(breaks) - 1
  
  bin_ids <- cut(Y_vec,
                 breaks         = breaks,
                 include.lowest = TRUE,
                 labels         = FALSE)
  
  total_loss_swap       <- torch::torch_zeros(1, device = X$device)
  total_loss_swap_b     <- torch::torch_zeros(1, device = X$device)
  total_loss_cov_diag   <- torch::torch_zeros(1, device = X$device)
  total_loss_cov_offdiag <- torch::torch_zeros(1, device = X$device)
  total_loss_marginal   <- torch::torch_zeros(1, device = X$device)
  total_loss_mean       <- torch::torch_zeros(1, device = X$device)
  total_loss            <- torch::torch_zeros(1, device = X$device)
  
  n_valid       <- 0L
  
  for (k in seq_len(n_bins_actual)) {
    
    idx <- which(bin_ids == k)
    if (length(idx) < min_bin_size) next
    
    # Subset to bin k
    idx_t <- torch::torch_tensor(idx,
                                 dtype  = torch::torch_long(),
                                 device = X$device)
    X_k  <- X$index_select(1L,  idx_t)   
    Xk_k <- Xk$index_select(1L, idx_t)  
    
    # Extract bin-specific target correlation for bin k
    target_corr_k <- NULL
    if (!is.null(object$pars$target_corr_by_bin)) {
      if (is.matrix(object$pars$target_corr_by_bin)) {
        target_corr_k <- object$pars$target_corr_by_bin[k, ]
      } else {
        target_corr_k <- object$pars$target_corr_by_bin[[k]]
      }
    }
    
    # ── Loss 1: MMD swap loss (within bin k) ────────────────────────
    L_swap <- .km_loss_mmd_single(object, X_k, Xk_k)
    
    # ── Loss 2: Second moment swap loss (within bin k) ──────────────
    L_swap_b <- .km_loss_swap_b_single(object, X_k, Xk_k)
    
    # ── Loss 3: Cross covariance and correlation betwen Xj and XTilde-j (within bin k) ─────────
    L_cov_result <- .km_loss_cov_single(object, X_k, Xk_k, target_corr = target_corr_k)
    L_cov_diag <- L_cov_result$diag
    L_cov_offdiag <- L_cov_result$offdiag
    
    # ── Loss 4: Marginal second moment loss (within bin k) ───────────
    L_marginal <- .km_loss_marginal_sm_single(object, X_k, Xk_k)
    
    # ── Loss 5: Mean matching loss (within bin k) ────────────────────
    L_mean <- .km_loss_mean_single(object, X_k, Xk_k)
    
    # add all loss components (for monitoring progress)
    total_loss_swap       <- total_loss_swap       + L_swap
    total_loss_swap_b     <- total_loss_swap_b     + L_swap_b
    total_loss_cov_diag   <- total_loss_cov_diag   + L_cov_diag
    total_loss_cov_offdiag <- total_loss_cov_offdiag + L_cov_offdiag
    total_loss_marginal   <- total_loss_marginal   + L_marginal
    total_loss_mean       <- total_loss_mean       + L_mean
    
    # Weighted sum for bin k
    bin_loss <- object$pars$lambda_swap     * L_swap     +
      object$pars$lambda_swap_b   * L_swap_b   +
      object$pars$lambda_cov      * (L_cov_diag + L_cov_offdiag) +
      object$pars$lambda_marginal * L_marginal +
      object$pars$lambda_mean     * L_mean
    
    total_loss <- total_loss + bin_loss 
    n_valid    <- n_valid + 1L
    
  }
  
  # Average across bins
  if (n_valid > 0L) {
    total_loss            <- total_loss            / n_valid
    total_loss_swap       <- total_loss_swap       / n_valid
    total_loss_swap_b     <- total_loss_swap_b     / n_valid
    total_loss_cov_diag   <- total_loss_cov_diag   / n_valid
    total_loss_cov_offdiag <- total_loss_cov_offdiag / n_valid
    total_loss_marginal   <- total_loss_marginal   / n_valid
    total_loss_mean       <- total_loss_mean       / n_valid
  }
  
  # Return list with total loss + individual components (cov split)
  list(
    total = total_loss,
    swap = total_loss_swap,
    swap_b = total_loss_swap_b,
    cov_diag = total_loss_cov_diag,
    cov_offdiag = total_loss_cov_offdiag,
    marginal = total_loss_marginal,
    mean = total_loss_mean
  )
}

train.KnockoffMachine <- function(object, X, Y, verbose = TRUE, ...) {
  
  # ── Input validation ────────────────────────────────────────────────
  if (!is.matrix(X)) {
    stop(sprintf(
      "X must be a matrix. Got: '%s'. Try X <- as.matrix(X).",
      paste(class(X), collapse = ", ")
    ))
  }
  if (is.null(Y)) {
    stop("Y is required. Pass your response variable as Y = your_vector.")
  }
  if (ncol(X) != object$p) {
    stop(sprintf(
      "ncol(X) = %d does not match machine$p = %d.",
      ncol(X), object$p
    ))
  }
  if (nrow(X) != length(Y)) {
    stop(sprintf(
      "nrow(X) = %d does not match length(Y) = %d.",
      nrow(X), length(Y)
    ))
  }
  if (!is.numeric(Y)) {
    stop(sprintf(
      "Y must be numeric. Got: '%s'. Try Y <- as.numeric(Y).",
      paste(class(Y), collapse = ", ")
    ))
  }
  
  # ── Convert to tensors ──────────────────────────────────────────────
  X_t <- torch::torch_tensor(X,            dtype = torch::torch_float())
  Y_t <- torch::torch_tensor(as.matrix(Y), dtype = torch::torch_float())
  
  # ── Validation set ────────────────────────────────────────────────
  n_total <- X_t$shape[1]
  has_val <- object$pars$test_size > 0
  
  if (has_val) {
    n_val       <- floor(object$pars$test_size * n_total)
    n_train     <- n_total - n_val
    all_idx     <- seq_len(n_total)
    val_idx_r   <- sample(all_idx, n_val)
    train_idx_r <- setdiff(all_idx, val_idx_r)
    
    val_idx_t   <- torch::torch_tensor(val_idx_r,   dtype = torch::torch_long())
    train_idx_t <- torch::torch_tensor(train_idx_r, dtype = torch::torch_long())
    
    X_val <- X_t$index_select(1L, val_idx_t)
    Y_val <- Y_t$index_select(1L, val_idx_t)
    X_t   <- X_t$index_select(1L, train_idx_t)
    Y_t   <- Y_t$index_select(1L, train_idx_t)
    
    if (verbose) {
      cat(sprintf("  Train: %d obs  |  Val: %d obs\n", n_train, n_val))
    }
  } else {
    n_train <- n_total
  }
  
  # ── Optimizer ───────────────────────────────────────────────────────
  optimizer <- torch::optim_adam(object$net$parameters, lr = object$pars$lr)
  
  object$net$train()
  
  # ── Epoch loop ──────────────────────────────────────────────────────
  for (epoch in seq_len(object$pars$epochs)) {
    
    running_loss_swap       <- 0.0
    running_loss_swap_b     <- 0.0
    running_loss_cov_diag   <- 0.0
    running_loss_cov_offdiag <- 0.0
    running_loss_marginal   <- 0.0
    running_loss_mean       <- 0.0
    running_loss_total      <- 0.0
    
    for (iter in seq_len(object$pars$epoch_length)) {
      
      # Random mini-batch with replacement
      idx_r   <- sample(n_train, object$pars$batch_size, replace = TRUE)
      idx_t   <- torch::torch_tensor(idx_r, dtype = torch::torch_long())
      X_batch <- X_t$index_select(1L, idx_t)
      Y_batch <- Y_t$index_select(1L, idx_t)
      
      optimizer$zero_grad()
      
      Xk_batch <- .km_generate(object, X_batch, Y_batch)
      loss_result <- .km_compute_loss(object, X_batch, Xk_batch, Y_batch)
      loss <- loss_result$total
      
      # NaN guard
      if (is.nan(loss$item())) {
        warning(sprintf(
          paste0("NaN loss at epoch %d, iter %d.\n",
                 "  Try: (1) lower lr, (2) scale your data, ",
                 "(3) check for NaN/Inf in X or Y."),
          epoch, iter
        ))
        return(invisible(object))
      }
      
      loss$backward()
      
      if (!is.null(object$pars$max_norm)) {
        torch::nn_utils_clip_grad_norm_(
          object$net$parameters,
          max_norm = object$pars$max_norm
        )
      }
      
      optimizer$step()
      
      # Accumulate individual loss components
      running_loss_total        <- running_loss_total        + loss_result$total$item()
      running_loss_swap         <- running_loss_swap         + loss_result$swap$item()
      running_loss_swap_b       <- running_loss_swap_b       + loss_result$swap_b$item()
      running_loss_cov_diag     <- running_loss_cov_diag     + loss_result$cov_diag$item()
      running_loss_cov_offdiag  <- running_loss_cov_offdiag  + loss_result$cov_offdiag$item()
      running_loss_marginal     <- running_loss_marginal     + loss_result$marginal$item()
      running_loss_mean         <- running_loss_mean         + loss_result$mean$item()
    }
    
    # ── LR milestone decay ──────────────────────────────────────────
    if (!is.null(object$pars$lr_milestones)) {
      if (epoch %in% object$pars$lr_milestones) {
        current_lr <- optimizer$param_groups[[1]]$lr
        new_lr     <- current_lr * object$pars$lr_decay
        optimizer$param_groups[[1]]$lr <- new_lr
        if (verbose) {
          cat(sprintf(
            "  [Epoch %d] LR decayed: %.6f -> %.6f\n",
            epoch, current_lr, new_lr
          ))
        }
      }
    }
    
    # Print progress
    if (verbose) {
      # Calculate averages for this epoch
      avg_loss_total      <- running_loss_total      / object$pars$epoch_length
      avg_loss_swap       <- running_loss_swap       / object$pars$epoch_length
      avg_loss_swap_b     <- running_loss_swap_b     / object$pars$epoch_length
      avg_loss_cov_diag   <- running_loss_cov_diag   / object$pars$epoch_length
      avg_loss_cov_offdiag <- running_loss_cov_offdiag / object$pars$epoch_length
      avg_loss_marginal   <- running_loss_marginal   / object$pars$epoch_length
      avg_loss_mean       <- running_loss_mean       / object$pars$epoch_length
      
      if (has_val) {
        val_loss_result <- .km_val_loss(object, X_val, Y_val)
        val_loss_total      <- val_loss_result$total$item()
        val_loss_swap       <- val_loss_result$swap$item()
        val_loss_swap_b     <- val_loss_result$swap_b$item()
        val_loss_cov_diag   <- val_loss_result$cov_diag$item()
        val_loss_cov_offdiag <- val_loss_result$cov_offdiag$item()
        val_loss_marginal   <- val_loss_result$marginal$item()
        val_loss_mean       <- val_loss_result$mean$item()
        
        # Print with individual components (cov split into diag and offdiag)
        cat(sprintf(
          "Epoch [%3d / %d]  TOTAL: train=%.4f  val=%.4f\n",
          epoch, object$pars$epochs, avg_loss_total, val_loss_total
        ))
        cat(sprintf(
          "  ├─ swap       : train=%.4f  val=%.4f\n",
          avg_loss_swap, val_loss_swap
        ))
        cat(sprintf(
          "  ├─ swap_b     : train=%.4f  val=%.4f\n",
          avg_loss_swap_b, val_loss_swap_b
        ))
        cat(sprintf(
          "  ├─ cov_diag   : train=%.4f  val=%.4f\n",
          avg_loss_cov_diag, val_loss_cov_diag
        ))
        cat(sprintf(
          "  ├─ cov_offdiag: train=%.4f  val=%.4f\n",
          avg_loss_cov_offdiag, val_loss_cov_offdiag
        ))
        cat(sprintf(
          "  ├─ marginal   : train=%.4f  val=%.4f\n",
          avg_loss_marginal, val_loss_marginal
        ))
        cat(sprintf(
          "  └─ mean       : train=%.4f  val=%.4f\n",
          avg_loss_mean, val_loss_mean
        ))
      } else {
        # Print without validation
        cat(sprintf(
          "Epoch [%3d / %d]  TOTAL: %.4f\n",
          epoch, object$pars$epochs, avg_loss_total
        ))
        cat(sprintf(
          "  ├─ swap       : %.4f\n",
          avg_loss_swap
        ))
        cat(sprintf(
          "  ├─ swap_b     : %.4f\n",
          avg_loss_swap_b
        ))
        cat(sprintf(
          "  ├─ cov_diag   : %.4f\n",
          avg_loss_cov_diag
        ))
        cat(sprintf(
          "  ├─ cov_offdiag: %.4f\n",
          avg_loss_cov_offdiag
        ))
        cat(sprintf(
          "  ├─ marginal   : %.4f\n",
          avg_loss_marginal
        ))
        cat(sprintf(
          "  └─ mean       : %.4f\n",
          avg_loss_mean
        ))
      }
      
      flush.console()
    }
  }
  
  invisible(object)
}

generate.KnockoffMachine <- function(object, X, Y, ...) {
  
  if (!is.matrix(X)) stop("X must be a matrix.")
  if (is.null(Y))     stop("Y is required.")
  if (ncol(X) != object$p)
    stop(sprintf("ncol(X) = %d != machine$p = %d", ncol(X), object$p))
  if (nrow(X) != length(Y))
    stop(sprintf("nrow(X) = %d != length(Y) = %d", nrow(X), length(Y)))
  
  object$net$eval()
  
  X_t  <- torch::torch_tensor(X,            dtype = torch::torch_float())
  Y_t  <- torch::torch_tensor(as.matrix(Y), dtype = torch::torch_float())
  
  Xk_t <- torch::with_no_grad(.km_generate(object, X_t, Y_t))
  
  as.matrix(Xk_t$cpu())
}

.km_generate <- function(object, X, Y) {
  n     <- X$shape[1]
  noise <- torch::torch_randn(
    n, object$pars$noise_dim, device = X$device
  ) * object$pars$noise_std
  object$net(X, noise, Y)
}

.km_compute_loss <- function(object, X, Xk, Y) {
  .km_loss_stratified_all(object, X, Xk, Y)
}



# ── Loss 1: MMD swap loss 
# This function divides observations into two disjoint subsets and computes:
# (1) MMD between Z1=(X1,Xk1) and Z2=(Xk2,X2)  [full swap]
# (2) MMD between Z1=(X1,Xk1) and Z3=(X2,Xk2)_swap  [partial random swap]
.km_loss_mmd_single <- function(object, X, Xk) {
  p_feat <- X$shape[2]
  n_total <- X$shape[1]
  
  # Divide observations into two disjoint batches
  n_half <- as.integer(n_total / 2)
  
  # Subset indices for each half 
  idx1 <- seq_len(n_half)
  idx2 <- (n_half + 1L):(n_total)
  
  idx1_tensor <- torch::torch_tensor(idx1, 
                                     dtype = torch::torch_long(),
                                     device = X$device)
  idx2_tensor <- torch::torch_tensor(idx2,
                                     dtype = torch::torch_long(),
                                     device = X$device)
  
  # Extract subsets (dim=1L means rows)
  X1   <- X$index_select(1L, idx1_tensor)
  Xk1  <- Xk$index_select(1L, idx1_tensor)
  X2   <- X$index_select(1L, idx2_tensor)
  Xk2  <- Xk$index_select(1L, idx2_tensor)
  
  # Create joint variables (concatenate along columns, dim=2L)
  # Z1 = [X1 | Xk1]  (concatenate column-wise)
  Z1 <- torch::torch_cat(list(X1, Xk1), dim = 2L)
  
  # Z2 = [Xk2 | X2]  (full swap)
  Z2 <- torch::torch_cat(list(Xk2, X2), dim = 2L)
  
  # Z3 = [X2 | Xk2] with random column swap
  Z3_base <- torch::torch_cat(list(X2, Xk2), dim = 2L)
  
  # Apply random swap: for each feature dimension, with prob 0.5, swap columns
  swap_inds <- which(runif(p_feat) < 0.5)
  
  if (length(swap_inds) > 0) {
    # Create a copy of Z3_base to modify
    Z3 <- Z3_base$clone()
    
    # Swap selected columns: move left half to right and vice versa
    for (idx in swap_inds) {
      col_left  <- idx            # Column in first p features (X2)
      col_right <- idx + p_feat   # Column in second p features (Xk2)
      
      # Swap columns
      temp <- Z3[, col_left]$clone()
      Z3[, col_left] <- Z3[, col_right]
      Z3[, col_right] <- temp
    }
  } else {
    Z3 <- Z3_base
  }
  
  # Compute MMD losses between:
  # (1) Z1 versus Z2 (discrepancy between (X,Xk) and (Xk,X))
  mmd_full <- mmd_loss(Z1, Z2, object$pars$alphas)
  
  # (2) Z1 versus Z3_swap (discrepancy between (X,Xk) and (X,Xk)_swap)
  mmd_swap <- mmd_loss(Z1, Z3, object$pars$alphas)
  
  # Return sum of both MMD losses
  mmd_full + mmd_swap
}


# ── Loss 2: Second moment swap loss (not used) ──────────────────────────────────────────
.km_loss_swap_b_single <- function(object, X, Xk) {
  p_feat <- X$shape[2]
  T      <- torch::torch_cat(list(X, Xk), dim = 2L)
  
  # Random partial swap
  col_order_random <- .km_swap_cols(p_feat)
  T_swap_random    <- T[, col_order_random]
  loss_random      <- second_moment_loss(T, T_swap_random)
  
  # Full swap
  col_order_full <- c((p_feat + 1L):(2L * p_feat), 1L:p_feat)
  T_swap_full    <- T[, col_order_full]
  loss_full      <- second_moment_loss(T, T_swap_full)
  
  loss_random + loss_full
}



# ── Loss 3: Cross covariance and de-correlation ─────────────────────────────────────
.km_loss_cov_single <- function(object, X, Xk, target_corr = NULL) {
  n <- X$shape[1]
  p <- X$shape[2]
  
  # Centre X and Xk
  X_c  <- X  - torch::torch_mean(X,  dim = 1L, keepdim = TRUE)
  Xk_c <- Xk - torch::torch_mean(Xk, dim = 1L, keepdim = TRUE)
  
  # Covariance matrices
  Sigma_X    <- torch::torch_mm(X_c$t(),  X_c)  / (n - 1L)
  Sigma_Xk   <- torch::torch_mm(Xk_c$t(), Xk_c) / (n - 1L)
  Sigma_X_Xk <- torch::torch_mm(X_c$t(),  Xk_c) / (n - 1L)
  
  # Off-diagonal: 
  diff_cov      <- Sigma_X_Xk - Sigma_X
  diag_mask     <- torch::torch_eye(p, device = X$device)
  off_mask      <- 1L - diag_mask
  off_diag_loss <- torch::torch_sum((diff_cov * off_mask)^2) / (torch::torch_sum(Sigma_X^2) + 1e-8)
  
  # Diagonal: 
  diag_cross <- torch::torch_diag(Sigma_X_Xk)
  var_X      <- torch::torch_diag(Sigma_X)
  var_Xk     <- torch::torch_diag(Sigma_Xk)
  corr_j     <- diag_cross / torch::torch_sqrt(var_X * var_Xk + 1e-8)
  
  # Determine target correlation
  if (!is.null(target_corr)) {
    target <- torch::torch_tensor(
      as.numeric(target_corr),
      dtype  = torch::torch_float(),
      device = X$device
    )
    diag_loss <- torch::torch_mean((corr_j - target)^2)
  } else {
    diag_loss <- torch::torch_mean(corr_j^2)
  }
  
  # Individual weights
  w_diag    <- if (is.null(object$pars$lambda_cov_diag))
    1.0 else object$pars$lambda_cov_diag
  w_offdiag <- if (is.null(object$pars$lambda_cov_offdiag))
    1.0 else object$pars$lambda_cov_offdiag
  
  # Return list with both components
  list(
    diag = w_diag * diag_loss,
    offdiag = w_offdiag * off_diag_loss
  )
}


# ── Loss 4: Marginal second moment loss ───────────────────────────────────────
.km_loss_marginal_sm_single <- function(object, X, Xk) {
  n <- X$shape[1]
  
  # Center X and Xk
  X_c  <- X  - torch::torch_mean(X,  dim = 1L, keepdim = TRUE)
  Xk_c <- Xk - torch::torch_mean(Xk, dim = 1L, keepdim = TRUE)
  
  # Covariance matrices
  Sigma_X  <- torch::torch_mm(X_c$t(),  X_c)  / (n - 1L)
  Sigma_Xk <- torch::torch_mm(Xk_c$t(), Xk_c) / (n - 1L)
  
  # Frobenius norm of difference / Frobenius norm of Sigma_X
  diff_frob_sq <- torch::torch_sum((Sigma_X - Sigma_Xk)^2)
  X_frob_sq    <- torch::torch_sum(Sigma_X^2)
  
  # cat(sprintf("Second moment loss = %f\n", diff_frob_sq / (X_frob_sq + 1e-8)))
  
  diff_frob_sq / (X_frob_sq + 1e-8)
}


# ── Loss 5: Mean matching loss ────────────────────────────────────────────────
.km_loss_mean_single <- function(object, X, Xk) {
  mu_X  <- torch::torch_mean(X,  dim = 1L)
  mu_Xk <- torch::torch_mean(Xk, dim = 1L)
  
  torch::torch_mean((mu_X - mu_Xk)^2)
}

.km_val_loss <- function(object, X_val, Y_val) {
  object$net$eval()
  loss_result <- NULL
  torch::with_no_grad({
    Xk_val     <- .km_generate(object, X_val, Y_val)
    loss_result <- .km_compute_loss(object, X_val, Xk_val, Y_val)
  })
  object$net$train()
  # Return the list of losses
  loss_result
}

# Random column swap permutation (used in second moment swap loss)
.km_swap_cols <- function(p) {
  swap_mask <- as.array(torch::torch_rand(p) > 0.5)
  swap_idx  <- which(swap_mask)
  col_order <- seq_len(2L * p)
  if (length(swap_idx) > 0L) {
    col_order[swap_idx]     <- swap_idx + p
    col_order[swap_idx + p] <- swap_idx
  }
  col_order
}


# Y converter 
.to_y_vec <- function(Y) {
  if (inherits(Y, "torch_tensor")) {
    if (length(Y$shape) > 1L && Y$shape[2] == 1L) {
      Y <- Y$squeeze(2L)
    }
    as.numeric(Y$cpu())
  } else {
    as.numeric(Y)
  }
}
