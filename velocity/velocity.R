library(tidyverse)
library(readr)

vel_cal_1 <- read_csv("excel/vel_cal_6.csv")

vel_cal_1_ss <- vel_cal_1 %>% 
  group_by(type) %>%
  summarise(mean_vel = mean(velocity),
            min_vel = min(velocity),
            max_vel = max(velocity),
            n_vel = n(),
            sd_vel = sd(velocity),
            se_vel = (sd_vel/sqrt(n_vel)),
            UL = (mean_vel + se_vel),
            LL = (mean_vel - se_vel))

vel_cal_1_ss <- vel_cal_1_ss %>% 
  mutate(mean_vel_1 = (mean_vel/1.62500),
         sd_vel_1 = (sd_vel/1.62500),
         se_vel_1 = (sd_vel_1/sqrt(n_vel)),
         UL_1 = (mean_vel_1 + se_vel_1),
         LL_1 = (mean_vel_1 - se_vel_1))

vel_cal_1_ss %>% 
  ggplot(aes(x = mean_vel_1, 
             y = type,
             color = type)) +
  geom_point(size = 6, alpha = 0.8, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = LL_1, xmax = UL_1, height = 0.2), show.legend = F) +
  labs(x = "Velocity of RNAP relative to HTkA-free conditions", 
       y = "HTkA condition of template") +
  ggtitle("Comparison of RNAP velocities at 60 seconds") +
  scale_y_discrete(limits = c("E_G", "E19K", "G17D", "E34A", "G52K","R11A", 
                              "WT", "T55L", "E3A","R20S","naked"),
                   labels = c("E_G"="E19K/G52K", "E19K"="E19K", "G17D"="G17D", 
                              "E34A"="E34A", "G52K"="G52K","R11A"="R11A",
                              "WT"="WT", "T55L"="T55L", "E3A"="E3A",
                              "R20S"="R20S","naked"="HTkA-free")) +
  scale_colour_manual(breaks = c("naked", "R20S","E3A","T55L", "WT","R11A",
                                 "G52K", "E34A", "G17D", "E19K", "E_G"), 
                      values = c("#004134", "greenyellow", "steelblue3", 
                                 "springgreen4", "grey68", "navy",
                                 "lightslateblue", "lightsteelblue", 
                                 "#FF6633", "purple4","darkviolet")) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 16),
        axis.title = element_text(family = '', face = 'bold', size = 14),
        axis.text = element_text(family = "", size = 11, color = "black"))

vel_60s <- read_csv("excel/vel_cal_60s.csv")

vel_60s_ss <- vel_60s %>% 
  group_by(type, TFS) %>%
  summarise(mean_vel = mean(velocity),
            min_vel = min(velocity),
            max_vel = max(velocity),
            n_vel = n(),
            sd_vel = sd(velocity),
            se_vel = (sd_vel/sqrt(n_vel)),
            UL = (mean_vel + se_vel),
            LL = (mean_vel - se_vel))

vel_60s_ss <- vel_60s_ss %>%
  mutate(mean_vel_1 = (mean_vel/1.625),
         sd_vel_1 = (sd_vel/1.625),
         se_vel_1 = (sd_vel_1/sqrt(n_vel)),
         UL_1 = (mean_vel_1 + se_vel_1),
         LL_1 = (mean_vel_1 - se_vel_1))

vel_60s_ss %>% group_by(TFS) %>% 
  ggplot(aes(x = mean_vel_1, 
             y = type,
             color = type, 
             shape = TFS, 
             linetype = TFS)) +
  geom_point(aes(size = TFS), 
             alpha = 0.8, 
             stroke = 1.5) +
  scale_size_manual(values = c(5,2)) +
  guides(color = "none", size = "none") +
  geom_errorbarh(aes(xmin = LL_1, xmax = UL_1, height = 0.2),
                 show.legend = FALSE) +
  labs(x = "Velocity of RNAP relative to HTkA-free conditions", 
       y = "HTkA condition of template") +
  ggtitle("Comparison of RNAP velocities at 60 seconds") +
  scale_y_discrete(limits = c("E_G", "E19K", "G17D", "E34A", "G52K","R11A", 
                              "WT", "T55L", "E3A","R20S","naked"),
                   labels = c("E_G"="E19K/G52K", "E19K"="E19K", "G17D"="G17D", 
                              "E34A"="E34A", "G52K"="G52K","R11A"="R11A",
                              "WT"="WT", "T55L"="T55L", "E3A"="E3A",
                              "R20S"="R20S","naked"="HTkA-free")) +
  scale_colour_manual(breaks = c("naked", "R20S","E3A","T55L", "WT","R11A",
                                 "G52K", "E34A", "G17D", "E19K", "E_G"), 
                      values = c("#004134", "greenyellow", "steelblue3", 
                                 "springgreen4", "grey68", "navy",
                                 "lightslateblue", "lightsteelblue", 
                                 "#FF6633", "purple4","darkviolet")) +
  scale_shape_manual(values = c(16,6), 
                     labels = c("- TFS", "+ TFS"),
                     guide = guide_legend(override.aes = list(size = c(5, 2)))) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 16),
        axis.title = element_text(family = '', face = 'bold', size = 14),
        axis.text = element_text(family = "", size = 11, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(family = '', face = 'bold', size = 11))


