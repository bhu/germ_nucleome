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

samps <- c("d4c7PGCLC","GSC",  'GSCLC')
rnm <- function(x) sub('ESC', 'mESC',  sub('PGC', ' mPGC', x)) 
samps <- rnm(samps)
clrs <- setNames(c('#d62728', '#9467bd', '#499894'), samps)
d %>%
  mutate(samp = factor(rnm(samp), samps)) %>%
  arrange(samp) %>%
  na.omit() %>%
  ggplot(aes(x = samp, y = n, fill = samp)) +
  geom_col(aes(group = rep), position = position_dodge(.5), width = .3, alpha = .8) +
  scale_fill_manual(values = clrs) +
  scale_y_continuous('# of CTCF peaks', expand = expansion(c(0, .05)), breaks = c(0, 5e4, 1e5), labels = c(0, '5e4', '1e5')) +
  facet_grid(.~'Reduced CTCF\nloss in GSCLCs', scales = 'free') +
  #coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        strip.clip = 'off',
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        plot.margin = margin(5, 20, 5, 5),
        strip.text = element_text(color = 'black', size = 13, family = 'Arial', face = 'bold'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed')) -> p
  ggsave('f6_k.pdf', p, height = 2.5, width = 2.2, device = cairo_pdf, bg = 'transparent')

