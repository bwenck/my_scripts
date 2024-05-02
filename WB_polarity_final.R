library(tidyverse)
library(emmeans)
library(broom)
library(ggpubr)
library(rstatix)

WB <- read_csv("pol_WB_v5.csv")

WB <- WB %>% 
  mutate(type = as.factor(type),
         type = fct_inorder(type))
  
WB_ss <- WB %>% 
  group_by(type) %>% 
  summarise(mean_quant = mean(avg),
            sd_quant = sd(avg),
            n_quant = n(),
            se_quant = (sd_quant/sqrt(n_quant)),
            UL = (mean_quant + se_quant),
            LL = (mean_quant - se_quant),
            .groups = "drop")

hist(WB$avg)
qqnorm(WB_avg$avg_quant)
qqline(WB_avg$avg_quant)

WB %>% identify_outliers(avg)
###test for normality
WB %>% shapiro_test(avg)
###test for equal variances
WB %>% levene_test(avg ~ type)

### There are no extreme outliers, the data is normally distributed, and there 
### are equal variances; I will use the Student t-test

LMfit <- lm(avg ~ type, data = WB)

summary(LMfit)
par(mfrow=c(2,2))
plot(LMfit)
anova(LMfit)

WB_stest <- WB %>% 
  t_test(avg ~ type, ref.group = "TS559")

WB_ss %>% group_by(type) %>% 
  ggplot(aes(x = type,
             y = mean_quant,
             color = type,
             fill = type)) +
  geom_col(alpha = 0.7) +
  guides(color = "none", fill = "none") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif),
               data = WB_stest, y.position = c(43,47,51),
               inherit.aes = FALSE) +
  xlab("Strain") +
  ylab("Absolute quantity of FttA (ng)") +
  geom_errorbar(aes(ymin = LL, ymax = UL, width = 0.1), 
                color = "black",
                alpha = 0.5) +
  scale_color_manual(breaks = c("TS559", "CM003", "Mut95", "CM005"),
                     labels = c("TS559", "CM003", "Mut95", "CM005"),
                     values = c("steelblue3", "springgreen4", "#FF6633",
                                "lightslateblue")) +
  scale_fill_manual(breaks = c("TS559", "CM003", "Mut95", "CM005"),
                    labels = c("TS559", "CM003", "Mut95", "CM005"),
                    values = c("steelblue3", "springgreen4", "#FF6633",
                               "lightslateblue"),
                    name = "Strain") +
  labs(tag = "Note: *"~italic(p)~"< .5, **"~italic(p)~"< .01, ***"~italic(p)~"< .001;"
       ~italic(p)~"values were determined using the Student's" ~italic(t)~"-test.") +
  theme_minimal() +
  theme(plot.margin = margin(1.0, 0.5, 1.0, 0.5, "cm"),
        plot.tag.position = c(0.5, -0.05),
        plot.tag = element_text(family = '', face = 'bold', size = rel(0.85)),
        axis.title = element_text(family = '', face = 'bold', size = 14),
        axis.text = element_text(family = "", face = 'italic',
                                 size = 12, color = "#263238"))
