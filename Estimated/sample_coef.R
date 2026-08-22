# sample beta and gamma 
coef_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/correct/coef/"

p <- 200

beta <- numeric(p)
gamma <- numeric(p)

n_nonnull_gamma <- 40
n_nonnull_beta <- 40

nonnull_loc_beta <- sample(1:p, n_nonnull_beta, replace = F)
nonnull_loc_gamma <- sample(1:p, n_nonnull_gamma, replace = F)

beta[nonnull_loc_beta] <- rnorm(n_nonnull_beta, 0, 0.4)
gamma[nonnull_loc_gamma] <- rnorm(n_nonnull_gamma, 0, 0.4)

write.table(beta, paste0(coef_dir, "beta.txt"), col.names = F, row.names = F)
write.table(gamma, paste0(coef_dir, "gamma.txt"), col.names = F, row.names = F)

# sample coef when there's an interaction 
coef_dir <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/result/correct/coef/"

p <- 200
beta_int <- numeric(p)
gamma_int <- numeric(p)

# sample beta
n_nonnull_beta <- 40
nonnull_loc_beta <- sample(1:p, n_nonnull_beta, replace = F)

beta_int[nonnull_loc_beta] <- rnorm(n_nonnull_beta, 0, 0.4)


# sample gamma
n_nonnull_gamma_mean <- 20
n_nonnull_gamma_int_x <- 15
n_nonnull_gamma_int_y <- 5

non_null_gamma_mean <- sample(1:p, n_nonnull_gamma_mean, replace = F)
non_null_gamma_interaction_x <- matrix(NA, nrow = n_nonnull_gamma_int_x * 2, ncol = 2)
non_null_gamma_interaction_x[,1] <- sample(1:p, n_nonnull_gamma_int_x * 2, replace = F)
non_null_gamma_interaction_x[,2] <- sample(1:p, n_nonnull_gamma_int_x * 2, replace = F)
non_null_gamma_interaction_y <- sample(1:p, n_nonnull_gamma_int_y, replace = F)
non_null_gamma_interaction_x <- unique(non_null_gamma_interaction_x)[1:n_nonnull_gamma_int_x, ]

gamma_mean <- rnorm(n_nonnull_gamma_mean, 0, 0.4)
gamma_interaction_x <- rnorm(n_nonnull_gamma_int_x, 0, 0.4)
gamma_interaction_y <- rnorm(n_nonnull_gamma_int_y, 0, 0.25)

gamma_int <- numeric(p)
gamma_int[non_null_gamma_mean] <- gamma_mean
write.table(beta_int, paste0(coef_dir, "beta_int.txt"), col.names = F, row.names = F)
write.table(gamma_int, paste0(coef_dir, "gamma_int.txt"), col.names = F, row.names = F)
write.table(non_null_gamma_interaction_x, paste0(coef_dir, "loc_gamma_int_x.txt"), col.names = F, row.names = F)
write.table(non_null_gamma_interaction_y, paste0(coef_dir, "loc_gamma_int_y.txt"), col.names = F, row.names = F)
write.table(gamma_interaction_x, paste0(coef_dir, "gamma_int_x.txt"), col.names = F, row.names = F)
write.table(gamma_interaction_y, paste0(coef_dir, "gamma_int_y.txt"), col.names = F, row.names = F)







