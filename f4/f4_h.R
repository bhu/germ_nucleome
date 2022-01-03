library(tidyverse)
library(pals)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

read_csv('../data/western/summary.csv') %>%
  mutate(samp = factor(cell, samps)) %>%
  na.omit() %>%
  filter(trg == 'Laminb1') %>%
  ggplot(aes(x = samp, y = ratio, color = cell)) +
  geom_point(size = 3) +
  facet_grid(.~'Lamin B1 quantification') +
  ylab(expression(frac('Lamin B1',beta*'-Actin'))) +
  scale_color_manual(values = clrs) +
  scale_y_continuous(breaks = c(1, 1.5), expand = expansion(.1)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey90'),
        strip.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) -> p
  ggsave('f4_h.pdf', height = 1.2, width = 3.1, device = cairo_pdf, bg = 'transparent')

