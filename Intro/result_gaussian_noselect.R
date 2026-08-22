# Fig 3 supplement
# FDR of standard knockoff when there's no selection
library(ggplot2)
library(dplyr)
library(tidyr)

# Load results
result_dir <- "/Results/Gaussian/"
fig_dir <-  "/Fig/Gaussian/"


nselect <- read.table(paste0(result_dir, "nselect_ko_NoSelection.txt"))
nfd <- read.table(paste0(result_dir, "nfd_ko_NoSelection.txt"))

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
result <- compute_fdr_power(nfd, nselect, nnonnull, fdr_levels, "Standard")

# Plot the FDP 
fig_fdp_noselect <- ggplot() + 
  geom_boxplot(aes(x = q, y = FDP), data = result$fdp, width = 0.25, fill = "#E69F00") + 
  geom_segment(aes(x = c(0.9, 1.9, 2.9, 3.9,4.9,5.9) - 0.1,
                   xend = c(1.1, 2.1, 3.1, 4.1, 5.1,6.1) + 0.1, 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

ggsave(fig_fdp_noselect, filename = paste0(fig_dir, "fdp_noselect.png"), width = 5, height = 4, unit = "in")

# Plot the power 
fig_power_noselect <- ggplot() + 
  geom_boxplot(aes(x = q, y = Power), data = result$power, width = 0.25) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"))

ggsave(fig_power_gaussian, filename = paste0(fig_dir, "power_noselect.png"), width = 5, height = 4, unit = "in")



