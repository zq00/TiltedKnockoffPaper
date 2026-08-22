rm(list = ls())

# Results for estimating the population distribution 
library(tidyverse)
library(purrr)
library(forcats)
library(ggpubr)

# Load results
result_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/estimated/"
fig_dir <-  "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/fig/"

# find all the filenames 
filenames <-  list.files(result_dir, pattern = "FDP")

# compile the number of selections and false discoveries in each simulation
result <- NULL
for(i in 1:length(filenames)){
  new_data <- read.table(paste0(result_dir, filenames[i]))
  colnames(new_data) <- c("nselect", "nfd", "q", "method")
  new_data$index <- i
  new_data$fdr_level <- rownames(new_data)
  
  result <- rbind(result, new_data)
}

# parameters
fdr_levels <- c("0.05", "0.1", "0.15", "0.2", "0.25","0.3")
nnonnull <- 40

# compute FDP and power
dat <- result %>% 
  mutate(fdp  = nfd / max(nselect, 1),
         power = (nselect - nfd) / nnonnull) %>% 
  mutate(method = case_when(
    method == "ko" ~ "No adjustment",
    method == "so" ~ "Known",
    method == "gaussian" ~ "Gaussian"
  )
  ) %>% 
  mutate(q = as.factor(q),
         method = fct_relevel(method, "No adjustment", "Known", "Gaussian")) 


# Define new color mapping
color_mapping <- c(
  "No adjustment" = "grey60",    # Red
  "Known" = "#56B4E9",            # Blue
  "Gaussian" = "#27C9B3"          # Green
)


# Plot FDP 
fig_fdp <- ggplot() + 
  geom_boxplot(aes(x = q, y = fdp, fill = method), data = dat, width = 0.6) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8, 4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
                   yend = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)), 
               color = "red", linetype = "solid") +
  scale_fill_manual(values = color_mapping) + 
  ylim(c(0, 0.65)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

# Plot power 
fig_power <- ggplot() + 
  geom_boxplot(aes(x = q, y = power, fill = method), data = dat, width = 0.6) + 
  ylim(c(0, 1)) + 
  scale_fill_manual(values = color_mapping) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 12))

# Combine figures with merged legend
combined_fig <- fig_fdp + fig_power + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

ggsave(combined_fig, filename = paste0(fig_dir, "fdp_power_est.png"), width = 10, height = 4, unit = "in")



