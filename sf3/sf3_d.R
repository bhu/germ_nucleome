library(tidyverse)
library(rtracklayer)
library(pals)

d <- list.files('../data/peaks/CTCF', pattern = 'narrowPeak', full.names = T) %>%
  split(., sub('_.*', '', basename(.))) %>%
  lapply(function(x) {
    lapply(x, import) %>%
      sapply(length) %>%
      tibble(n = .) %>%
      mutate(rep = 1:n())
  }) %>%
  bind_rows(.id = 'samp')

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC',  sub('PGC', '\nmPGC', x)) 
samps <- rnm(samps)
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)
d %>%
  mutate(samp = factor(rnm(samp), rev(samps))) %>%
  arrange(samp) %>%
  na.omit() %>%
  ggplot(aes(x = samp, y = n, fill = samp)) +
  geom_col(aes(group = rep), position = position_dodge(.5), width = .3, alpha = .8) +
  scale_fill_manual(values = clrs) +
  scale_y_continuous('# of CTCF peaks', expand = expansion(c(0, .05)), breaks = c(0, 5e4, 1e5), labels = c(0, '5e4', '1e5')) +
  facet_grid(.~'Loss of CTCF\nbinding in GSCs', scales = 'free') +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        strip.clip = 'off',
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, family = 'Arial', face = 'bold'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed')) +
  ggsave('sf3_d.pdf', height = 5.85, width = 2, device = cairo_pdf, bg = 'transparent')

