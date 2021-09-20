library(tidyverse)
library(pals)
library(ggnewscale)

samps <- c('ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}

load('../data/tads/10.methods.rda')

d %>%
  count(samp, method) %>%
  mutate(samp = factor(samp, samps)) %>%
  na.omit() %>%
  group_by(method) %>%
  mutate(rank = as.character(rank(n))) %>%
  ungroup() %>%
  ggplot(aes(x = rank)) +
  geom_bar(aes(fill = samp), alpha = .7) +
  scale_fill_manual(values = clrs, labels = rnm) +
  scale_x_discrete(breaks = as.character(1:5),
                   labels = c('1\nLeast\nboundaries', '2', '3', '4', '5\nMost\nboundaries')) +
  scale_y_continuous('# of supporting algorithms', breaks = seq(0,10,2),
                     expand = expansion(0)) +
  facet_grid(.~'Insulation ranking') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color=  'grey70', linetype = 'dashed'),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        legend.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face =  'bold')) +
  ggsave('f2_j.pdf', height = 2.9, width = 3.5)

