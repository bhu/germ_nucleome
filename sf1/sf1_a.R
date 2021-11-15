library(tidyverse)
library(pals)
library(ggpubr)
library(ggdist)
library(gghalves)

samps <- c('ESC', 'EpiLC', 'GSC')
clrs <- setNames(tableau20()[c(1,3,9)], samps)

load('../data/image/FISH.res.rda')
d <- FISH.res %>%
  filter(Variable == 'Volume' & chr == 'chr16') %>%
  mutate(samp = factor(cell, samps),
         y = Value)


anns <- compare_means(y ~ samp, data = d) %>%
  arrange(group2, group1) %>%
  mutate(y.position = 30 + (1:n()) * 3) 


ggplot(d, aes(x = samp, y = y, color = samp, fill = samp)) +
  geom_half_point(aes(color = samp), range_scale = 0.4, size = .25, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(side = 'right', alpha = .5, color = NA, scale = .8) +
  stat_pointinterval() +
  stat_pvalue_manual(data = anns, label = 'p.signif', inherit.aes = F) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  ylab(expression('chr1 FISH signal volume ('*mu*'m'^3*')')) +
  scale_y_continuous(breaks = seq(0,30,10)) +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', x)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(l=10),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey80')) -> p

ggsave('sf1_a.pdf', p, height = 3.1, width = 2.8)

