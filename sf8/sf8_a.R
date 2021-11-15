library(tidyverse)
library(pals)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20()[seq(1, 9, 2)], samps)

read_csv('../data/western/summary.csv') %>% 
  mutate(cell = factor(cell, samps)) %>% 
  filter(grepl('CTCF', trg)) %>% 
  na.omit() %>%
  mutate(trg = c(CTCF = 'Whole cell lysate', CTCF_chromatin = 'Chromatin fraction')[trg]) %>%
  ggplot(aes(x = cell, y = ratio)) + 
  geom_line(aes(linetype = trg, group = trg)) +
  geom_point(aes(shape = trg, y = ratio * 100), size = 3) +
  geom_point(aes(shape = trg, color = cell), show.legend = F, size = 3) +
  scale_color_manual(values = clrs) +
  ylab('Normalized\nCTCF signal') +
  coord_cartesian(ylim = c(.4, 1.4)) +
  facet_grid(.~'CTCF quantification') +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey90'),
        legend.title = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black')) -> p
ggsave('sf8_a.pdf', p, height = 1.7, width = 3.58)

