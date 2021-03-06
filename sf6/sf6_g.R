library(tidyverse)
library(pals)
library(ggpubr)
library(ggdist)
library(gghalves)
library(patchwork)

samps <- c('ESC', 'EpiLC','GSC') 
rnm <- function(x) sub('ESC', 'mESC', x)
clrs <- setNames(tableau20()[c(1,3,9)], samps)

load('../data/image/FISH.res.rda')

d <- FISH.res %>%
  filter(chr == 'chrY' & Variable == 'Sphericity') %>%
  mutate(samp = factor(cell, rev(samps)),
         y = Value)

anns <- compare_means(y ~ samp, data = d) %>%
  arrange(group1, group2) %>%
  mutate(y.position = .85 + (1:n()) * .1) 


ggplot(d, aes(x = samp, y = y, color = samp, fill = samp)) +
  geom_half_point(aes(color = samp), range_scale = 0.4, size = .25, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(side = 'right', alpha = .5, color = NA) +
  stat_pointinterval() +
  stat_pvalue_manual(data = anns, label = 'p.signif', inherit.aes = F,
                     coord.flip = T) +
  coord_flip() +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  ylab('chrY FISH sphericity') +
  scale_x_discrete(labels = rnm) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.y = element_line(color = 'grey80')) -> p


d2 <- FISH.res %>%
  filter(chr == 'chrY' & Variable == 'Volume') %>%
  mutate(samp = factor(cell, rev(samps)),
         y = Value) %>%
  na.omit()

anns2 <- compare_means(y ~ samp, data = d) %>%
  arrange(group2, group1) %>%
  mutate(y.position = 40 + (1:n()) * 10) 


ggplot(d2, aes(x = samp, y = y, color = samp, fill = samp)) +
  geom_half_point(aes(color = samp), range_scale = 0.4, size = .25, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(side = 'right', alpha = .5, color = NA) +
  stat_pointinterval() +
  stat_pvalue_manual(data = anns2, label = 'p.signif', inherit.aes = F,
                     coord.flip = T) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  ylab(expression('Signal volume ('*mu*'m'^3*')')) +
  scale_x_discrete(labels = rnm) +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.y = element_line(color = 'grey80'))  -> p2


{wrap_plots(p2,p,nrow=1) &
    theme(plot.background = element_blank()) } %>%
  ggsave('sf6_g.pdf', ., height = 2.2, width = 8)

