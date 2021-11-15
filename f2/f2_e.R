library(data.table)
library(tidyverse)
library(pals)
library(scales)

samps <- c('ESC','EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}
clrs <- setNames(tableau20(9)[seq(1, 9, 2)], samps)

d <- list.files('../data/abc', full.names = T, pattern = 'csv.gz') %>%
  setNames(., sub('.csv.gz', '', basename(.))) %>%
  grep('C.csv', ., invert = T, value = T) %>%
  lapply(fread) %>%
  rbindlist(idcol = 'samp') %>%
  count(samp) %>%
  mutate(type = factor(sub('_.*', '', samp), samps)) %>%
  na.omit() %>%
  arrange(type)


ggplot(d, aes(x = type, y = n)) +
  geom_col(aes(fill = type, group = samp), position = position_dodge(.5), width = .3) +
  scale_fill_manual(values = clrs) +
  scale_y_continuous('# of predicted E-P pairs',
                     expand = expansion(c(0, .05)),
                     labels = function(x) sub('0000', 'e4', x)) +
  facet_grid(.~'Activity-by-contact pairs') +
  scale_x_discrete(labels = rnm) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed')) -> p
  ggsave('f2_e.pdf', p, height = 3.1, width = 3.3)

  