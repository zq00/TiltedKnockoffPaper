rm(list = ls())
source("/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/src_V2/00_config.R")

library(rlang, lib.loc = "/home/qianzhao_umass_edu/R/x86_64-pc-linux-gnu-library/4.5/")
library(tibble)
library(forcats)
library(dplyr)
library(ggplot2)
library(patchwork)

fig_dir <- "/work/pi_qianzhao_umass_edu/TiltedKnockoff/result/fig/"

# ============================================================
# Collect results from all replicates
# ============================================================
all_tilted   <- list()
all_standard <- list()
all_so <- list()

for (b in 1:n_sim) {
  res_file <- file.path(sim_dir(b), "results.RData")
  if (file.exists(res_file)) {
    load(res_file)
    all_tilted[[b]]   <- result_tilted
    all_standard[[b]] <- result_standard
    all_so[[b]] <- result_so
  } else {
    cat("WARNING: missing results for replicate", b, "\n")
  }
}

df_tilted   <- bind_rows(all_tilted) %>% 
  mutate(method = "Deep")
df_standard <- bind_rows(all_standard) %>% 
  mutate(method = "No adjustment")
df_so <- bind_rows(all_so) %>% 
  mutate(method  = "SO")

# ============================================================
# Summarize: mean FDP and mean number of selections across sims
# ============================================================
summary_tilted <- df_tilted %>%
  group_by(q) %>%
  summarise(
    mean_fdp     = mean(fdp, na.rm = TRUE),
    sd_fdp       = sd(fdp, na.rm = TRUE),
    mean_nselect = mean(nselect, na.rm = TRUE),
    n_reps       = n(),
    .groups      = "drop"
  )

summary_standard <- df_standard %>%
  group_by(q) %>%
  summarise(
    mean_fdp     = mean(fdp, na.rm = TRUE),
    sd_fdp       = sd(fdp, na.rm = TRUE),
    mean_nselect = mean(nselect, na.rm = TRUE),
    n_reps       = n(),
    .groups      = "drop"
  )

summary_so <- df_so %>%
  group_by(q) %>%
  summarise(
    mean_fdp     = mean(fdp, na.rm = TRUE),
    sd_fdp       = sd(fdp, na.rm = TRUE),
    mean_nselect = mean(nselect, na.rm = TRUE),
    n_reps       = n(),
    .groups      = "drop", 
  )

cat("=== Deep Tilted Knockoff (averaged over", n_sim, "sims) ===\n")
print(summary_tilted)

cat("\n=== Standard Knockoff (averaged over", n_sim, "sims) ===\n")
print(summary_standard)

cat("\n=== SO Knockoff (averaged over", n_sim, "sims) ===\n")
print(summary_so)


# ============================================================
# Plot results 
# ============================================================

# parameters
nnonnull <- p * 0.2

# compute FDP and power
df_all <- bind_rows(df_tilted,
                    df_standard,
                    df_so) %>% 
  filter(q <= 0.3 )

dat <- df_all %>% 
  mutate( power = (nselect - nfd) / nnonnull) %>% 
  mutate(q = as.factor(q),
         method = fct_relevel(method, "No adjustment", "SO", "Deep")) %>% 
  mutate(method = recode(method,
                         "No adjustment" = "No adjustment",
                         "Deep" = "Deep Tilted",
                         "SO" = "SO (binned)"))


# Define color mapping
color_mapping <- c("No adjustment" = "grey60", 
                   "Deep Tilted" = "#3C7711", 
                   "SO (binned)" = "#FE6100")

# Plot FDP 
fig_fdp <- ggplot() + 
  geom_boxplot(aes(x = q, y = fdp, fill = method), data = dat, width = 0.6) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8, 4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
                   yend = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)), 
               color = "red", linetype = "solid") +
  scale_fill_manual(values = color_mapping) + 
  ylim(c(0, 1)) + 
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

# Combine figures with merged legend (same style as your Gaussian figures)
combined_fig <- fig_fdp + fig_power + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

ggsave(combined_fig, filename = paste0(fig_dir, "fdp_power_deep.png"), width = 10, height = 4, unit = "in")

# ============================================================
# Save aggregated results
# ============================================================
save(df_tilted, df_standard, df_so, summary_tilted, summary_standard, summary_so,
     file = file.path(path_data, "aggregated_results.RData"))

cat("\nAggregated results saved to:", file.path(path_data, "aggregated_results.RData"), "\n")

