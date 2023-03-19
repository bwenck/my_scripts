library(tidyverse)
library(readr)
library(viridis)

ec <- read_csv("excel/ext_chase_msmts_long.csv")

ec <- ec %>% mutate(time = as.factor(time),
                    mins = as.factor(mins),
                    type = as.factor(type),
                    type = fct_inorder(type),
                    bin = as.factor(bin),
                    bin = fct_rev(bin))

ec %>% ggplot(aes(x = mins,
                  y = percent,
                  fill = bin)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = "RNA \nlength (nt)",
                     alpha = 0.8,
                     labels = c("209-231","179-208","149-178","119-148",
                                "89-118","59-88","58")) +
  labs(x = "Time (min)", y = "Percentage of transcripts") +
  ggtitle("The E19K/G52K variant impedes TEC progression in extended timepoints") +
  facet_wrap(~ type) +
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
