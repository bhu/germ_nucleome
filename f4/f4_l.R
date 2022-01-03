library(tidyverse)
library(pals)
library(ggdist)
library(gghalves)

samps <- c('EpiLC', 'GSC') 
clrs <- setNames(tableau20()[c(3,9)], samps)

read_csv('../data/image/Major_satellite_count.csv') %>% 
  mutate(Cell = factor(Cell, samps)) %>%
  ggplot(aes(x = Cell, y = 100 * interior_num / all_num)) +
  coord_cartesian(xlim = c(1.1,2)) +
  geom_half_point(aes(color = Cell), range_scale = 0.4, size = 1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = Cell), alpha = .5, scale = .5) +
  stat_pointinterval(aes(color = Cell)) +
  ylab('% of interior pericentromere') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank()) -> p
  ggsave('f4_l.pdf', p, height = 2.2, width = 3.7)
