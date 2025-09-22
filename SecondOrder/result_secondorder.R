# Results for the second order tilted knockoff
library(tidyverse)
library(purrr)
library(ggpubr)

# Load results
result_dir <- "/Results/SecondOrder/"
fig_dir <-  "/Fig/SecondOrder/"

# find all the filenames 
filenames <-  list.files(result_dir)

# compile the number of selections and false discoveries in each simulation
result <- NULL
for(i in 1:length(filenames)){
  new_data <- read.table(paste0(result_dir, filenames[i]))
  new_data$index <- i
  new_data$fdr_level <- rownames(new_data)
  
  result <- rbind(result, new_data)
}

# parameters
fdr_levels <- c("0.05", "0.1", "0.15", "0.2", "0.25","0.3")
nnonnull <- 200 * 0.1

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
result <- result %>% mutate(
  nselect_ko_max1 = map_dbl(nselect_ko, ~max(.x, 1)),
  fdp_ko = ndf_ko / nselect_ko_max1,
  power_ko = (nselect_ko - ndf_ko) / nnonnull, 
  nselect_so_max1 = map_dbl(nselect_so, ~max(.x, 1)),
  fdp_so = nfd_so / nselect_so_max1,
  power_so = (nselect_so - nfd_so) / nnonnull
) 

# Plot the FDP 
fdp <- result %>% 
  mutate(
  Standard = ndf_ko / map_dbl(nselect_ko, ~max(.x, 1)),
  `Tilted (Second Order)` = nfd_so / map_dbl(nselect_so, ~max(.x, 1)),
  ) %>% 
  select(Standard, `Tilted (Second Order)`, fdr_level) %>% 
  rename("No adjustment" = Standard ) %>% 
  pivot_longer(cols = c(`No adjustment`, `Tilted (Second Order)`), names_to = "Method") 

fig_fdp_secondorder <- ggplot() + 
  geom_boxplot(aes(x = fdr_level, y = value, fill = Method), data = fdp, width = 0.25) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8,4.8,5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2,6.2), 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  ylim(c(0, 0.75)) + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

ggsave(fig_fdp_secondorder, filename = paste0(fig_dir, "fdp_secondorder.png"), width = 5, height = 4, unit = "in")



# Plot the power 
power <- result %>% 
  mutate(
    Standard = (nselect_ko - ndf_ko) / nnonnull,
    `Tilted (Second Order)` = (nselect_so - nfd_so) / nnonnull,
  ) %>% 
  select(Standard, `Tilted (Second Order)`, fdr_level) %>% 
  rename("No adjustment" = Standard ) %>% 
  pivot_longer(cols = c(`No adjustment`, `Tilted (Second Order)`), names_to = "Method") 

# color palette at http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
fig_power_secondorder <- ggplot() + 
  geom_boxplot(aes(x = fdr_level, y = value, fill = Method), data = power, width = 0.25) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

ggsave(fig_power_secondorder, filename = paste0(fig_dir, "power_secondorder.png"), width = 5, height = 4, unit = "in")

fig_secondorder <- ggarrange(fig_fdp_secondorder, fig_power_secondorder, ncol = 2, common.legend = T)

ggsave(fig_secondorder, filename = paste0(fig_dir, "fig_secondorder.png"), width = 9, height = 4, unit = "in")


