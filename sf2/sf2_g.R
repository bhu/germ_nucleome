library(data.table)
library(tidyverse)
library(pals)
library(ggdist)
library(ggpubr)
library(ggnewscale)
library(ggforce)


samps <- c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC")
clrs <- setNames(tableau20(10)[seq(1, 10, 2)], samps)

d <- list.files('../data/abc', full.names = T, pattern = 'csv.gz') %>%
  setNames(., sub('.csv.gz', '', basename(.))) %>%
  grep('C.csv', ., invert = T, value = T) %>%
  lapply(fread) %>%
  rbindlist(idcol = 'samp') %>%
  mutate(type = factor(sub('_.*', '', samp), rev(samps))) %>%
  na.omit() %>%
  arrange(type) %>%
  mutate(rep = sub('.*_', '', samp))

p <- compare_means(distance ~ type, d, method = 'wilcox.test') %>%
  filter(group1 == 'd4c7PGCLC' | group2 == 'd4c7PGCLC') %>%
  mutate(y.position = 7 + (1:n()) * .3)


ggplot(d, aes(y = type, x = distance, color = type, fill = type, group = CellType)) +
  geom_boxplot(position = position_dodge(.8), width = .5, notch = T) +
  stat_summary(geom = "crossbar", width = 0.3, fatten = 0, color = "white", position = position_dodge(.8),
               fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_zoom(xlim = c(4.83, 5), zoom.size = 1) +
  scale_y_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGCLC', ' mPGCLC', x))) +
  xlab('Distance between ABC-based E-P pairs') +
  scale_x_continuous(
    labels = function(r) {
      if (diff(log10(r))[1] < 1) {
        c('70', '80', '90', '100kb')
      } else {
        c('1kb', '10kb', '100kb', '1mb')
      }
    },
    breaks = function(r) {
    if (diff(log10(r)) < 1) {
      c(7e4, 8e4, 9e4, 1e5)
    } else {
      c(1e3, 1e4, 1e5, 1e6)
    }
  }, trans = 'log10') + 
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed')) +
  ggsave('sf2_g.pdf', height = 4.6, width = 3.9)

