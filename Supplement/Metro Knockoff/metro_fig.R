# 5/8/2024
# Figures for the results for the deep knockoff
library(tidyverse)
library(ggpubr)

result_dir <- "/Results/metro/"
fig_dir <-  "/Fig/metro/"

# 1 -- Compute power and FDR
nselect_ko <- as.matrix(read.table(paste0(result_dir, "nselect_ko.txt")))
nselect_so <- as.matrix(read.table(paste0(result_dir, "nselect_so.txt")))
nselect_metro <- as.matrix(read.table(paste0(result_dir, "nselect_metro.txt")))

nfd_ko <- as.matrix(read.table(paste0(result_dir, "nfd_ko.txt")))
nfd_so <- as.matrix(read.table(paste0(result_dir, "nfd_so.txt")))
nfd_metro <- as.matrix(read.table(paste0(result_dir, "nfd_metro.txt")))

dat <- tibble(
  method = rep(c("ko", "so", "metro"), each = 250 * 6),
  q = rep(c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), times = 250 * 3),
  nselect = c(as.vector(t(nselect_ko)), as.vector(t(nselect_so)), 
              as.vector(t(nselect_metro))),
  nfd = c(as.vector(t(nfd_ko)), as.vector(t(nfd_so)), 
          as.vector(t(nfd_metro)))
) %>% 
  mutate(FDP = nfd / map_dbl(nselect, ~max(.x, 1))) %>% 
  mutate(power =  (nselect - nfd) / 15) %>% 
  mutate(q = as.factor(q)) %>% 
  mutate(Method = factor(method, levels = c("ko", "so", "metro"),
                         labels = c("No adjustment", "SO tilted", "Metro")))

dat %>% group_by(method, q) %>% 
  summarize(mean(FDP))
# 2 -- Main texts
# 2.1 -- FDP
g_fdp_metro <- ggplot() + 
  geom_boxplot(aes(x = q , y = FDP, fill = Method), width = 0.3, data = dat) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8,4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  ylim(c(0, 0.75)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#785EF0")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.2 -- Power
g_power_metro <- ggplot() + 
  geom_boxplot(aes(x = q , y = power, fill = Method), width = 0.3, data = dat) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#785EF0")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.3 -- put the two plots together
g_metro <- ggarrange(g_fdp_metro, g_power_metro, ncol = 2, common.legend = T)

ggsave(g_metro, filename = paste0(fig_dir, "metro.png"), width = 10, height = 4, unit = "in")

