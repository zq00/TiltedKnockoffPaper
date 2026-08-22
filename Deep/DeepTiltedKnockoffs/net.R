# Deep neural net with FiLM architecture 

knockoff_net <- torch::nn_module(
  classname = "KnockoffNet",
  
  initialize = function(p,
                        noise_dim  = as.integer(p),
                        hidden_dim = 2L * as.integer(p),
                        num_layers = 6L,
                        y_dim      = 1L) {
    
    self$p         <- as.integer(p)
    self$noise_dim <- as.integer(noise_dim)
    self$y_dim     <- as.integer(y_dim)
    hidden_dim     <- as.integer(hidden_dim)
    num_layers     <- as.integer(num_layers)
    
    # Input block: takes [X, noise] only — Y is applied via FiLM 
    input_dim  <- self$p + self$noise_dim
    self$fc_in <- torch::nn_linear(input_dim, hidden_dim)
    self$bn_in <- torch::nn_batch_norm1d(hidden_dim)
    
    # Hidden blocks
    n_hidden <- max(0L, num_layers - 1L)
    self$fc_hidden <- torch::nn_module_list(
      lapply(seq_len(n_hidden), function(i)
        torch::nn_linear(hidden_dim, hidden_dim))
    )
    self$bn_hidden <- torch::nn_module_list(
      lapply(seq_len(n_hidden), function(i)
        torch::nn_batch_norm1d(hidden_dim))
    )
    
    # FiLM generators 
    n_film <- num_layers   # one per block
    self$film_generators <- torch::nn_module_list(
      lapply(seq_len(n_film), function(i)
        torch::nn_linear(self$y_dim, 2L * hidden_dim))  # outputs [gamma | beta]
    )
    
    self$fc_out <- torch::nn_linear(hidden_dim, self$p)
    self$act    <- torch::nn_leaky_relu(negative_slope = 0.1)
  },
  
  forward = function(x, noise = NULL, y = NULL) {
    n <- x$shape[1]
    
    if (is.null(noise)) {
      noise <- torch::torch_randn(n, self$noise_dim, device = x$device)
    }
    
    # Ensure y is (n, y_dim)
    if (!is.null(y) && length(y$shape) == 1L) {
      y <- y$unsqueeze(2L)
    }
    
    # apply FiLM modulation using the i-th FiLM generator
    apply_film <- function(h, film_idx) {
      if (!is.null(y)) {
        film_out <- self$film_generators[[film_idx]](y)      # (n, 2*hidden_dim)
        gamma    <- film_out[, 1:(film_out$shape[2] %/% 2)]  # (n, hidden_dim)
        beta     <- film_out[, ((film_out$shape[2] %/% 2) + 1):film_out$shape[2]]
        h <- gamma * h + beta                                 # scale and shift
      }
      h
    }
    
    # Input block
    h <- torch::torch_cat(list(x, noise), dim = 2L)  
    h <- self$bn_in(self$fc_in(h))
    h <- apply_film(h, 1L)                           
    h <- self$act(h)
    
    # Hidden blocks
    for (i in seq_along(self$fc_hidden)) {
      h <- self$bn_hidden[[i]](self$fc_hidden[[i]](h))
      h <- apply_film(h, i + 1L)                      
      h <- self$act(h)
    }
    
    self$fc_out(h)
  }
)