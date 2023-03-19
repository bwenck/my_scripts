library(tidyverse)

ppp <- read_csv("excel/all_reps_ppp.csv")

ppp_ss <- ppp %>% 
  group_by(Time, Type) %>% 
  summarise(mean_per = mean(Percent),
            sd_per = sd(Percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

ppp_ss %>% ggplot(aes(x = Time, 
                      y = mean_per,
                      color = Type)) +
  geom_line() +
  geom_point(size = 3, alpha = 0.7)+
  geom_errorbar(aes(ymin = LL, ymax = UL), width = 5, alpha = 0.5) +
  labs(x = "Time (sec)", y = "Percent of TECs (70 nt)") +
  ggtitle("Prominent pause position (70 nt)") +
  scale_colour_manual(breaks = c("Naked", "WT","G17D", "E19K","G52K", "E19K/G52K", 
                                 "R20S","T55L","E3A", "R11A","E34A"), 
                      values = c("#004134", "grey68", "#FF6633", 
                                 "purple4", "lightslateblue", "darkviolet", 
                                 "greenyellow", "springgreen4", "steelblue3", 
                                 "navy", "lightsteelblue"),
                      name = "Condition",
                      labels = c("HTkA-free","WT","G17D", "E19K","G52K", "E19K/G52K", 
                                 "R20S","T55L","E3A", "R11A","E34A")) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 16),
                        axis.title = element_text(family = '', 
                                                  face = 'bold',
                                                  size = 12),
                        legend.title = element_text(family = '', face = 'bold', 
                                                    size = 12))

ppp_TFS <- read_csv("TFS/all_TFS_ppp.csv")

ppp_TFS_ss <- ppp_TFS %>% 
  group_by(time, type) %>% 
  summarise(mean_per = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

ppp_TFS_ss %>% ggplot(aes(x = time, 
                      y = mean_per,
                      color = type)) +
  geom_line() +
  geom_point(size = 3, alpha = 0.7, shape = 6, stroke = 1.5) +
  geom_errorbar(aes(ymin = LL, ymax = UL), width = 5, alpha = 0.5) +
  labs(x = "Time (sec)", y = "Percent of TECs (70 nt)") +
  ggtitle("Prominent pause position +TFS (70 nt)") +
  scale_colour_manual(breaks = c("naked", "WT","G17D", "E19K","G52K", "E19K/G52K", 
                                 "R20S","T55L","E3A", "R11A","E34A"), 
                      values = c("#004134", "grey68", "#FF6633", 
                                 "purple4", "lightslateblue", "darkviolet", 
                                 "greenyellow", "springgreen4", "steelblue3", 
                                 "navy", "lightsteelblue"),
                      name = "Condition",
                      labels = c("HTkA-free","WT","G17D", "E19K","G52K", "E19K/G52K", 
                                 "R20S","T55L","E3A", "R11A","E34A")) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 16),
        axis.title = element_text(family = '', face = 'bold', size = 12),
        legend.title = element_text(family = '', face = 'bold', 
                                    size = 12))


