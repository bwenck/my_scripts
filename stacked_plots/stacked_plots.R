library(tidyverse)
library(viridis)
library(viridisLite)

trx_all <- read_csv("excel/hist_meas_all.csv")

trx_all <- trx_all %>%
  mutate(time = as.factor(time),
         type = as.factor(type),
         type = fct_inorder(type),
         bin = as.factor(bin),
         bin = fct_rev(bin))

trx_all_ss <- trx_all %>% 
  group_by(type, time, bin) %>% 
  summarise(mean_per = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

trx_all_ss %>% ggplot(aes(x = time,
                          y = mean_per,
                          fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (sec)", y = "Percentage of transcripts") +
  ggtitle("Variant chromatin alters the percentage of full-length RNA products") +
  facet_wrap(~ type, scales = "free_x") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', face = 'bold', size = 11),
        legend.title = element_text(family = '', face = 'bold', size = 9),
        legend.text = element_text(size = 8),
        legend.position = "left",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 9))


TFS_trx_all <- read_csv("TFS/TFS_meas_all.csv")

TFS_trx_all <- TFS_trx_all %>%
  mutate(time = as.factor(time),
         type = as.factor(type),
         type = fct_inorder(type),
         bin = as.factor(bin),
         bin = fct_rev(bin))

TFS_trx_all_ss <- TFS_trx_all %>% 
  group_by(type, time, bin) %>% 
  summarise(mean_per = mean(percent),
            sd_per = sd(percent),
            n_per = n(),
            se_per = (sd_per/sqrt(n_per)),
            UL = (mean_per + se_per),
            LL = (mean_per - se_per))

TFS_trx_all_ss %>% ggplot(aes(x = time,
                              y = mean_per,
                              fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (sec)", y = "Percentage of transcripts") +
  ggtitle("Variant chromatin alters the percentage of full-length RNA products (+TFS)") +
  facet_wrap(~ type, scales = "free_x") +
  theme_minimal() +
  theme(plot.title = element_text(family = '', face = 'bold', size = 12),
        axis.title = element_text(family = '', 
                                  face = 'bold', size = 11),
        legend.title = element_text(family = '', 
                                    face = 'bold', size = 9),
        legend.text = element_text(size = 8),
        legend.position = "left",
        strip.text.x = element_text(family = '', 
                                    face = 'bold', size = 9))