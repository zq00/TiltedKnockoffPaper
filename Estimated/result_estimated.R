# Results for the second order tilted knockoff
rm(list = ls())

library(rlang, lib.loc = "/home/qianzhao_umass_edu/R/x86_64-pc-linux-gnu-library/4.5/")
library(tidyverse)
library(purrr)
library(ggpubr)
library(forcats)
library(patchwork)


result_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/correct/"
fig_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/fig/"

# sample size
n <- 1600

# load and combine results from all 4 batches
n_batches <- 4

result_vars <- c("nselect_ko", "nfd_ko", "nselect_ko_case_control", "nfd_ko_case_control",
                 "nselect_so", "nfd_so", "nselect_est", "nfd_est", "nselect_ko_ipw",
                 "nfd_ko_ipw", "nselect_missing", "nfd_missing", "nselect_probit",
                 "nfd_probit", "nselect_xgboost", "nfd_xgboost", "nselect_ko_case_control", "nfd_ko_case_control")

# Load first batch as base
load(paste0(result_dir, "n_", n, "_batch_1_V2.RData"))
combined_results <- lapply(result_vars, function(var) get(var))
names(combined_results) <- result_vars

# Load and combine remaining batches
for(batch in 2:n_batches) {
  load(paste0(result_dir, "n_", n, "_batch_", batch, "_V2.RData"))
  for(var in result_vars) {
    combined_results[[var]] <- rbind(combined_results[[var]], get(var))
  }
}

# Put combined results back in global environment
for(var in result_vars) {
  assign(var, combined_results[[var]], envir = .GlobalEnv)
}

# parameters
fdr_levels <- as.character(q)
nnonnull <- 40

# Function to compute FDP and power
# same as in the intro_result (can separate this file into a src script)
compute_fdr_power <- function(nfd, nselect, nnonull, fdr_levels, method){
  fdp <- as_tibble(nfd / apply(nselect, c(1,2), function(t) max(t, 1)))
  power <- as_tibble((nselect - nfd) / nnonnull)
  
  colnames(fdp) <- fdr_levels
  colnames(power) <- fdr_levels
  
  fdp = pivot_longer(fdp, cols = all_of(fdr_levels),
                     names_to = "q", values_to = "FDP")
  fdp$Method <- method
  
  power <- pivot_longer(power, cols = all_of(fdr_levels),
                        names_to = "q", values_to = "Power")
  power$Method <- method
  
  list(fdp = fdp, power = power)
}

# Compute FDP and power
result_ko <- compute_fdr_power(nfd_ko, nselect_ko, nnonnull, fdr_levels, "No adjustment")
result_tilted <- compute_fdr_power(nfd_so, nselect_so, nnonnull, fdr_levels, "Known")
result_ko_ipw <- compute_fdr_power(nfd_ko_ipw, nselect_ko_ipw, nnonnull, fdr_levels, "Knockoff + IPW")
result_est <- compute_fdr_power(nfd_est, nselect_est, nnonnull, fdr_levels, "Estimated")
result_missing<- compute_fdr_power(nfd_missing, nselect_missing, nnonnull, fdr_levels, "Missing")
result_probit <- compute_fdr_power(nfd_probit, nselect_probit, nnonnull, fdr_levels, "Incorrect")
result_xgboost <- compute_fdr_power(nfd_xgboost, nselect_xgboost, nnonnull, fdr_levels, "XGBoost")
result_case_control <- compute_fdr_power(nfd_ko_case_control, nselect_ko_case_control, nnonnull, fdr_levels, "Case/control")

# color palette
# standard knockoff (no adjustment) #B5AAEA 
# KO + IPW: #CC6677
# tilted knockoff (known/exact) #3C7711
# tilted knockoff: estimated #49B3A1
# missing: #88CCEE
# incorrect: #DDCC77
# xgboost: #AA4499

# combine the FDP 
fdp_combined <- rbind(
  result_ko$fdp,
  result_tilted$fdp,
  result_est$fdp,
  result_ko_ipw$fdp, 
  result_missing$fdp,
  result_probit$fdp,
  result_xgboost$fdp,
  result_case_control$fdp
) %>% 
  filter(q != "0.05") %>% 
  mutate(Method = fct_relevel(Method, "No adjustment", "Knockoff + IPW", "Known", "Estimated", "XGBoost", "Missing", "Incorrect", "Case/control")) 

# Average FDP of each method
fdp_combined %>% 
  group_by(Method, q) %>% 
  summarize(avg_fdp = mean(FDP)) %>% View()

# Power
power_combined <- rbind(
  result_ko$power,
  result_tilted$power,
  result_est$power,
  result_ko_ipw$power,
  result_missing$power,
  result_probit$power,
  result_xgboost$power,
  result_case_control$power
) %>% 
  filter(q != "0.05") %>% 
  mutate(Method = fct_relevel(Method, "No adjustment", "Knockoff + IPW", "Known", "Estimated", "XGBoost", "Missing", "Incorrect", "Case/control")) 

power_combined %>% 
  group_by(Method, q) %>% 
  summarize(avg_fdp = mean(Power)) %>% View()


# Define the method order (as specified in your code)
method_order <- c("No adjustment", "Knockoff + IPW", "Case/control", "Known", "Estimated", "XGBoost", "Missing", "Incorrect")

# Define colorblind-friendly color mapping (Okabe-Ito palette - tested for colorblindness)
color_mapping <- c(
  "No adjustment" = "grey60",      
  "Knockoff + IPW" = "#E69F00",     
  "Known" = "#56B4E9",              
  "Estimated" = "#009E73",          
  "XGBoost" = "#F0E442",            
  "Missing" = "#D55E00",            
  "Incorrect" = "#CC79A7",          
  "Case/control" = "#0072B2"        
)

# Prepare FDP data
fdp_combined <- fdp_combined %>% 
  mutate(Method = fct_relevel(Method, method_order))

# Prepare Power data (filter and reorder as in your code)
power_combined <- power_combined %>% 
  filter(q != "0.05") %>% 
  mutate(Method = fct_relevel(Method, method_order))

# Plot FDP
fig_fdp_gaussian <- ggplot() + 
  geom_boxplot(aes(x = q, y = FDP, fill = Method), data = fdp_combined, width = 0.6) + 
  geom_segment(aes(x = c(0.75, 1.75, 2.75, 3.75, 4.75),
                   xend = c(1.25, 2.25, 3.25, 4.25, 5.25), 
                   y = c( 0.1, 0.15, 0.2, 0.25, 0.3),
                   yend = c( 0.1, 0.15, 0.2, 0.25, 0.3)), color = "red", linetype = "solid") + 
  scale_fill_manual(values = color_mapping) + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

# Plot Power
fig_power_gaussian <- ggplot() + 
  geom_boxplot(aes(x = q, y = Power, fill = Method), data = power_combined, width = 0.6) + 
  ylim(c(0, 1)) + 
  scale_fill_manual(values = color_mapping) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

# Figure 4-6 in the supplement (change n value to get different figures)
combined_fig <- fig_fdp_gaussian + fig_power_gaussian + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

ggsave(combined_fig, filename = paste0(fig_dir, "fdp_power_correct_", n, ".png"), width = 12, height = 5, unit = "in")

# Figure 4 in the paper (uses only second order knockoff with known or estimated selection probability)
fdp_comparison <- fdp_combined %>% 
  filter(Method %in% c("No adjustment", "Known", "Estimated"))

power_comparison <- power_combined %>% 
  filter(Method %in% c("No adjustment", "Known", "Estimated"))

fig_fdp_comparison <- ggplot() + 
  geom_boxplot(aes(x = q, y = FDP, fill = Method), data = fdp_comparison, width = 0.45) + 
  geom_segment(aes(x = c(0.75, 1.75, 2.75, 3.75, 4.75),
                   xend = c(1.25, 2.25, 3.25, 4.25, 5.25), , 
                   y = c( 0.1, 0.15, 0.2, 0.25, 0.3),
                   yend = c( 0.1, 0.15, 0.2, 0.25, 0.3)), color = "red", linetype = "solid") + 
  scale_fill_manual(values = color_mapping) + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

fig_power_comparison <- ggplot() + 
  geom_boxplot(aes(x = q, y = Power, fill = Method), data = power_comparison, width = 0.45) + 
  scale_fill_manual(values = color_mapping) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

combined_fig <- fig_fdp_comparison + fig_power_comparison + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

ggsave(combined_fig, filename = paste0(fig_dir, "fdp_power_correct_small_", n, ".png"), width = 10, height = 4, unit = "in")

