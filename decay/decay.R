library(tidyverse)
library(emmeans)
library(broom)

decay <- read_csv("excel/no_TFS_decay_1.csv")
decay_1 <- read_csv("excel/no_TFS_decay_1.csv")
decay_2 <- read_csv("TFS/TFS_decay.csv")


decay <- decay %>% 
  mutate(TFS = as.factor(TFS),
         type = as.factor(type))

decay_ss <- decay %>% 
  group_by(type, time, TFS) %>% 
  summarize(mean_per = mean(percent),
            min_per = min(percent),
            max_per = max(percent),
            n_per = n(),
            sd_per = sd(percent),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

decay_aov <- aov(mean_per ~ type * TFS * time, data = decay_ss)
summary(decay_aov)
tidy(decay_aov)

LMfit <- lm(mean_per ~ type, data = decay_ss)

anova(LMfit)
emout <- emmeans(LMfit, ~type)
pairs(emout, adjust = "none")

decay_1_ss <- decay_1 %>% 
  group_by(type, time) %>% 
  summarize(mean_per = mean(percent),
            min_per = min(percent),
            max_per = max(percent),
            n_per = n(),
            sd_per = sd(percent),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

decay_2_ss <- decay_2 %>% 
  group_by(type, time) %>% 
  summarize(mean_per = mean(percent),
            min_per = min(percent),
            max_per = max(percent),
            n_per = n(),
            sd_per = sd(percent),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))


decay_1_ss %>% 
  ggplot(aes(x = time, y = mean_per, color = type)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line() +
  geom_errorbar(aes(ymin = LL, ymax = UL, width = 4), alpha = 0.6) +
  labs(x = "Time (sec)", 
       y = "Percentage of complexes") +
  ggtitle("Percentage of decay from 58 nt") +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                  "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue")) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold',),
        legend.title = element_blank(),
        legend.text = element_text(family = '', face = 'bold'))

decay_2_ss %>% 
  ggplot(aes(x = time, y = mean_per, color = type)) +
  geom_point(size = 3, alpha = 0.6, shape = 6, stroke = 1.5) +
  geom_line() +
  geom_errorbar(aes(ymin = LL, ymax = UL, width = 4), alpha = 0.6) +
  labs(x = "Time (sec)", 
       y = "Percentage of complexes") +
  ggtitle("Percentage of decay from 58 nt (+TFS)") +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","#FF6633","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue")) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold',),
        legend.title = element_blank(),
        legend.text = element_text(family = '', face = 'bold'))

decay_ss %>% 
  group_by(TFS, time, type) %>% 
  ggplot(aes(y = type, 
             x = mean_per,
             color = type,
             shape = TFS)) +
  geom_point(alpha = 0.8, 
             stroke = 1.5) +
  scale_size_manual(values = c(2,2)) +
  guides(color = FALSE, size = FALSE) +
  facet_wrap(~ time) +
  geom_errorbarh(aes(xmin = LL, xmax = UL, height = 0.5), 
                alpha = 0.3) +
  labs(x = "Mean percentage of complexes at +58", 
       y = "Condition of template") +
  ggtitle("Comparison of decay from +58 (-/+ TFS)") +
  scale_y_discrete(limits = c("E34A","R11A","E3A","T55L","R20S","E19K/G52K",
                              "G52K","E19K","G17D", "WT","HTkA-free")) +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","violetred1","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue")) +
  scale_shape_manual(values = c(16,6),
                     labels = c("- TFS", "+ TFS"),
                     guide = guide_legend(override.aes = list(size = c(5, 2)))) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold',),
        legend.title = element_blank(),
        legend.text = element_text(family = '', face = 'bold'))

#Statistics
hist(decay_1$percent) 
hist(decay_2$percent)
hist(decay$percent)

hist(decay_1_group$mean)
hist(decay_2_group$mean)

boxplot(decay_1$percent)
boxplot(decay_2$percent)

mean(decay_1$percent)
mean(decay_2$percent)

qqnorm(decay_1$percent)
qqline(decay_1$percent)

qqnorm(decay_2$percent)
qqline(decay_2$percent)

qqnorm(decay$percent)
qqline(decay$percent)

t.test(decay_1$percent, decay_2$percent)
t.test(decay$percent)

decay_1_group <- decay_1_ss %>% 
  group_by(time, type) %>% 
  summarize(mean = mean(mean_per))

decay_2_group <- decay_2_ss %>% 
  group_by(time, type) %>% 
  summarize(mean = mean(mean_per))

t.test(decay_1_group$mean, decay_2_group$mean)

#Welch Two Sample t-test

#data:  decay_1$percent and decay_2$percent
#t = 9.2069, df = 308.06, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  6.621986 10.221832
#sample estimates:
#  mean of x mean of y 
#27.65409  19.23218 

#Welch Two Sample t-test

#data:  decay_1_group$mean and decay_2_group$mean
#t = 2.6373, df = 7.7224, p-value = 0.03076
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1.011734 15.832084
#sample estimates:
#  mean of x mean of y 
#27.65409  19.23218 

decay_means <- read_csv("excel/decay_means.csv")

decay_means %>% 
  group_by(TFS, time) %>%
  ggplot(aes(y = mean,
             x = type,
             fill = TFS)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  facet_wrap(~ time)
  
decay_ss %>% 
  group_by(TFS, time, type) %>% 
  ggplot(aes(x = mean_per, 
             y = type,
             fill = TFS)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbarh(aes(xmin = LL, xmax = UL, height = 0.3), 
                 alpha = 0.3, position = position_dodge(0.9), color = "black") +
  facet_wrap(~ time) +
  labs(x = "Mean percentage of complexes at +58", 
       y = "Condition of template") +
  ggtitle("Comparison of decay from +58 (-/+ TFS)") +
  scale_y_discrete(limits = c("E34A","R11A","E3A","T55L","R20S","E19K/G52K",
                              "G52K","E19K","G17D", "WT","HTkA-free")) +
  scale_fill_manual(values = c("mediumturquoise", "darkmagenta"),
                    labels = c("- TFS", "+ TFS")) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold',),
        legend.title = element_blank(),
        legend.text = element_text(family = '', face = 'bold'))

decay_ss %>% 
  ggplot(aes(x = time, 
             y = mean_per, 
             color = type, 
             shape = TFS,
             linetype = TFS)) +
  geom_point(aes(size = TFS), 
             alpha = 0.6,
             stroke = 1.5) +
  scale_size_manual(values = c(3,2)) +
  guides(size = FALSE, linetype = FALSE) +
  geom_line() +
  geom_errorbar(aes(ymin = LL, ymax = UL, width = 4), alpha = 0.4) +
  labs(x = "Time (sec)", 
       y = "Percentage of complexes") +
  ggtitle("Percentage of decay from 58 nt") +
  scale_colour_manual(breaks = c("HTkA-free", "WT","G17D","E19K","G52K",
                                 "E19K/G52K","R20S","T55L","E3A", "R11A",
                                 "E34A"), 
                      values = c("#004134", "grey68","violetred1","purple4",
                                 "lightslateblue","darkviolet","greenyellow", 
                                 "springgreen4","steelblue3","navy",
                                 "lightsteelblue")) +
  scale_shape_manual(values = c(16,6), 
                     labels = c("- TFS", "+ TFS"),
                     guide = guide_legend(override.aes = list(size = c(5, 2)))) +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold'),
        axis.title = element_text(family = '', face = 'bold',),
        legend.title = element_blank(),
        legend.text = element_text(family = '', face = 'bold'))

