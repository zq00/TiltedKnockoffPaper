# Simulation result in Supplement section 2 and paper section 3.2
rm(list = ls())

library(rlang, lib.loc = "/home/qianzhao_umass_edu/R/x86_64-pc-linux-gnu-library/4.5/")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 

# Load results
result_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/intro/"
fig_dir <-  "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/fig/"

load(paste0(result_dir, "intro_V2.RData"))

# parameters
fdr_levels <- c("0.05", "0.1", "0.15", "0.2", "0.25","0.3")
nnonnull <- 400 * 0.1

# Function to compute FDP and power
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
result_tilted <- compute_fdr_power(nfd_tilted, nselect_tilted, nnonnull, fdr_levels, "Tilted (Exact)")
result_ipw <- compute_fdr_power(nfd_ipw, nselect_ipw, nnonnull, fdr_levels, "IPW")
result_ipw_ko <- compute_fdr_power(nfd_ko_ipw, nselect_ko_ipw, nnonnull, fdr_levels, "KO+IPW")
result_ko_case_control <- compute_fdr_power(nfd_ko_case_control, nselect_ko_case_control, nnonnull, fdr_levels, "Case/Control KO")

# Define the method order
method_order <- c("IPW", "No adjustment", "Case/Control KO", "KO+IPW", "Tilted (Exact)")

# Define color mapping 
color_mapping <- c(
  "No adjustment" = "#E69F00",
  "Tilted (Exact)" = "#56B4E9",
  "IPW" = "#009E73",
  "Case/Control KO" = "#F0E442",
  "KO+IPW" = "#0072B2"
)

# Plot the FDP 
# combine the FDP 
fdp_combined <- rbind(
  result_ko$fdp,
  result_tilted$fdp,
  result_ipw$fdp,
  result_ipw_ko$fdp,
  result_ko_case_control$fdp
) %>% 
  mutate(Method = factor(Method, levels = method_order))

fig_fdp_gaussian <- ggplot() + 
  geom_boxplot(aes(x = q, y = FDP, fill = Method), data = fdp_combined, width = 0.45) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8, 4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
                   yend = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)), color = "red", linetype = "solid") + 
  scale_fill_manual(values = color_mapping) + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

# Plot the power 
# combine the power
power_combined <- rbind( 
  result_ko$power,
  result_tilted$power,
  result_ipw$power,
  result_ipw_ko$power,
  result_ko_case_control$power
) %>% 
  mutate(Method = factor(Method, levels = method_order))


fig_power_gaussian <- ggplot() + 
  geom_boxplot(aes(x = q, y = Power, fill = Method), data = power_combined, width = 0.45) + 
  ylim(c(0, 1)) + 
  scale_fill_manual(values = color_mapping) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

# Supplement Figure 2
combined_fig <- fig_fdp_gaussian + fig_power_gaussian + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

ggsave(combined_fig, filename = paste0(fig_dir, "fdp_power_exact.png"), width = 11, height = 4.5, unit = "in")


# plot the p-value in one simulation 
bin_loc <- seq(0, 1, by = 0.05)

# Supplement Figure 1
fig_pval <- ggplot() + 
  geom_histogram(aes(x = p_val[-nonnull_beta]), 
                 bins = 20, color = "grey50", fill = "grey",
                 breaks = bin_loc) + 
  xlab("Null p-values") + 
  ylab("Counts") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

ggsave(fig_pval, filename = paste0(fig_dir, "p_val.png"), width = 5, height = 4, unit = "in")

# Figure 2 in paper
fdp_comparison <- fdp_combined %>% 
  filter(Method %in% c("No adjustment", "Tilted (Exact)"))

fig_fdp_comparison <- ggplot() + 
  geom_boxplot(aes(x = q, y = FDP, fill = Method), data = fdp_comparison, width = 0.45) + 
  geom_segment(aes(x = c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75),
                   xend = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), , 
                   y = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
                   yend = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)), color = "red", linetype = "solid") + 
  scale_fill_manual(values = c("No adjustment" = "#E69F00", "Tilted (Exact)" = "#56B4E9")) + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

ggsave(fig_fdp_comparison, filename = paste0(fig_dir, "fdp_exact_small.png"), width = 5, height = 4, unit = "in")


