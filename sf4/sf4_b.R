library(tidyverse)
library(pals)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20()[seq(1, 9, 2)], samps)

read_csv('../data/western/summary.csv') %>% 
  mutate(cell = factor(cell, samps)) %>% 
  filter(grepl('K9', trg)) %>%
  na.omit() %>%
  mutate(trg = paste0('H3', trg)) %>%
  ggplot(aes(x = cell, y = ratio)) + 
  geom_line(aes(linetype = trg, group = trg)) +
  geom_point(aes(shape = trg, y = ratio * 100), size = 3) +
  geom_point(aes(shape = trg, color = cell), show.legend = F, size = 3) +
  scale_color_manual(values = clrs) +
  ylab('H3K9me / H3') +
  coord_cartesian(ylim = c(.2, 1.6)) +
  facet_grid(.~'H3K9me quantification') +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey90'),
        legend.title = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,.8),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(size = 11),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black')) -> p
ggsave('sf4_b.pdf', p, height = 2.4, width = 3.41)
