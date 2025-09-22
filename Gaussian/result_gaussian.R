# Fig 2
# Boxplot of the FDR and power when the sample is from a case-control study 
library(ggplot2)
library(dplyr)
library(tidyr)

# Load results
result_dir <- "/Results/Gaussian/"
fig_dir <-  "/Fig/Gaussian/"
  
nselect_ko <- read.table(paste0(result_dir, "nselect_ko.txt"))
nfd_ko <- read.table(paste0(result_dir, "nfd_ko.txt"))
nselect_tilted <- read.table(paste0(result_dir, "nselect_tilted.txt"))
nfd_tilted <- read.table(paste0(result_dir, "nfd_tilted.txt"))

# parameters
fdr_levels <- c("0.05", "0.1", "0.15", "0.2", "0.25","0.3")
nnonnull <- 400 * 0.1

# Function to compute FDP and power
compute_fdr_power <- function(nfd, nselect, nnonull, fdr_levels, method){
  fdp <- tibble(nfd / apply(nselect, c(1,2), function(t) max(t, 1)))
  power <- tibble((nselect - nfd) / nnonnull)
  
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
result_ko <- compute_fdr_power(nfd_ko, nselect_ko, nonnull, fdr_levels, "No adjustment")
result_tilted <- compute_fdr_power(nfd_tilted, nselect_tilted, nonnull, fdr_levels, "Tilted (Exact)")

# Plot the FDP 
# combine the FDP 
fdp_combined <- rbind(
  result_ko$fdp,
  result_tilted$fdp
)

fig_fdp_gaussian <- ggplot() + 
  geom_boxplot(aes(x = q, y = FDP, fill = Method), data = fdp_combined, width = 0.25) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8,4.8,5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2,6.2), 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
      axis.text = element_text(size = 15, color = "black"),
      legend.position="top",
      legend.text=element_text(size=12))

ggsave(fig_fdp_gaussian, filename = paste0(fig_dir, "fdp_exact.png"), width = 5, height = 4, unit = "in")

# Plot the power 
# combine the power
power_combined <- rbind( 
  result_ko$power,
  result_tilted$power
)

fig_power_gaussian <- ggplot() + 
  geom_boxplot(aes(x = q, y = Power, fill = Method), data = power_combined, width = 0.25) + 
  ylim(c(0, 1)) + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

ggsave(fig_power_gaussian, filename = paste0(fig_dir, "power_exact.png"), width = 5, height = 4, unit = "in")



