# 3/15/2024
# Figures for the results for the deep knockoff
library(tidyverse)
library(ggpubr)

result_dir <- "/Results/deep/"
fig_dir <-  "/Fig/deep/"



# 1 -- Compute power and FDR
nselect_ko <- as.matrix(read.table(paste0(result_dir, "nselect_ko.txt")))
nselect_so <- as.matrix(read.table(paste0(result_dir, "nselect_so.txt")))
nselect_deep <- as.matrix(read.table(paste0(result_dir, "nselect_deep.txt")))

nfd_ko <- as.matrix(read.table(paste0(result_dir, "nfd_ko.txt")))
nfd_so <- as.matrix(read.table(paste0(result_dir, "nfd_so.txt")))
nfd_deep <- as.matrix(read.table(paste0(result_dir, "nfd_deep.txt")))

dat <- tibble(
  method = rep(c("ko", "so", "deep"), each = 500 * 6),
  q = rep(c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), times = 500 * 3),
  nselect = c(as.vector(t(nselect_ko)), as.vector(t(nselect_so)), 
              as.vector(t(nselect_deep))),
  nfd = c(as.vector(t(nfd_ko)), as.vector(t(nfd_so)), 
          as.vector(t(nfd_deep)))
) %>% 
  mutate(FDP = nfd / map_dbl(nselect, ~max(.x, 1))) %>% 
  mutate(power =  (nselect - nfd) / 20) %>% 
  mutate(q = as.factor(q)) %>% 
  mutate(Method = factor(method, levels = c("ko", "so", "deep"),
                         labels = c("No adjustment", "SO tilted", "Deep")))


# 2 -- Main texts
# 2.1 -- FDP
g_fdp_deep <- ggplot() + 
  geom_boxplot(aes(x = q , y = FDP, fill = Method), width = 0.3, data = dat) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8,4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  ylim(c(0, 0.7)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#D81B60")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.2 -- Power
g_power_deep <- ggplot() + 
  geom_boxplot(aes(x = q , y = power, fill = Method), width = 0.3, data = dat) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#D81B60")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.3 -- put the two plots together
g_deep <- ggarrange(g_fdp_deep, g_power_deep, ncol = 2, common.legend = T)

ggsave(g_deep, filename = paste0(fig_dir, "deep.png"), width = 10, height = 4, unit = "in")

