# Oct 11, 2023
# Figures for the results for the estimated selection probabilities

library(tidyverse)
library(ggpubr)

result_dir <- "/Results/estimated/"
fig_dir <-  "/Fig/estimated/"

# 1 -- load FDP and power results
names <- list.files(result_dir)
indices <- sapply(names, function(t) strsplit(t, "_")[[1]][2])

dat <- NULL
for(i in 1:length(indices)){
  new_dat <- read.table(paste0(result_dir, "FDP_", indices[i]))
  new_dat$ind <- i
  print(new_dat[12, ])
  dat <- rbind(dat, new_dat)
}
colnames(dat) <- c("nselect", "nfd", "q", "method", "ind")

dat <- dat %>% 
  mutate(FDP = nfd / map_dbl(nselect, ~max(.x, 1))) %>% 
  mutate(power =  (nselect - nfd) / 40) %>% 
  mutate(q = as.factor(q))


# 2 -- Main texts
# 2.1 -- FDP
dat_main <- dat %>%  
  filter(method %in% c("ko", "so", "logistic")) %>% 
  mutate(Method = factor(method, levels = c("ko", "so", "logistic"),
                         labels = c("No adjustment", "SO tilted (known)", "SO tilted (estimated)")))


g_main_1 <- ggplot() + 
  geom_boxplot(aes(x = q , y = FDP, fill = Method), width = 0.3, data = dat_main) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8,4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  ylim(c(0, 0.8)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.2 -- Power
g_main_2 <- ggplot() + 
  geom_boxplot(aes(x = q , y = power, fill = Method), width = 0.3, data = dat_main) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.3 -- put the two plots together
g_main <- ggarrange(g_main_1, g_main_2, ncol = 2, common.legend = T)

ggsave(g_main, filename = paste0(fig_dir, "estimated_logistic.png"), width = 10, height = 4, unit = "in")

# 3 -- figure in the appendix 
dat_app <- dat %>%  
  filter(method %in% c("ko", "lambda-min", "lambda-ko-0.25")) %>% 
  mutate(Method = factor(method, levels = c("ko", "lambda-min", "lambda-ko-0.25"),
                         labels = c("No adjustment", "SO tilted (L1 lambda.min)", "SO tilted (KO q = 0.25)")))

g_app_1 <- ggplot() + 
  geom_boxplot(aes(x = q , y = FDP, fill = Method), width = 0.3, data = dat_app) + 
  geom_segment(aes(x = c(0.8, 1.8, 2.8, 3.8,4.8, 5.8),
                   xend = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), 
                   y = c(0.05, 0.1, 0.15,0.2,0.25,0.3),
                   yend = c(0.05, 0.1, 0.15,0.2,0.25,0.3)), color = "red", linetype = "solid") + 
  ylim(c(0, 0.8)) + 
  xlab("Target FDR") + 
  ylab("FDP") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.2 -- Power
g_app_2 <- ggplot() + 
  geom_boxplot(aes(x = q , y = power, fill = Method), width = 0.3, data = dat_app) + 
  ylim(c(0, 1)) + 
  xlab("Target FDR") + 
  ylab("Power") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        legend.position="top",
        legend.text=element_text(size=12))

# 2.3 -- put the two plots together
g_app <- ggarrange(g_app_1, g_app_2, ncol = 2, common.legend = T)

ggsave(g_app, filename = paste0(fig_dir, "/estimated_logistic_appendix.png"), width = 10, height = 4, unit = "in")




