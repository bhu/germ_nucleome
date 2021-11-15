library(tidyverse)
library(pals)
library(ggrepel)
library(patchwork)

load('../data/contactdecay/germline.smooth.rda')

clrs <- setNames(c(tableau20()[seq(1,9,2)], brewer.set2(6)[-2]),
                 c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC',
                   'E11.5_PGC', 'E13.5_mPGC', 'E13.5_fPGC', 'ICM', 'ESC_serum'))
rnm <- function(x) {
  sub('_', ' ', x) %>%
    sub('serum', '(serum)', .) %>%
    sub('^ESC', 'mESC', .) %>%
    sub('PGCLC', ' mPGCLC', .) %>%
    sub('ESC$', 'ESC (2i)', .)
}

d <- o %>%
  pivot_longer(-sep, names_to = 'samp', values_to = 'v') %>%
  filter(sep > 1000) %>%
  filter(samp %in% names(clrs)) %>%
  group_by(samp) %>%
  arrange(sep) %>%
  mutate(x = log10(sep),
         y = log10(v),
         x2 = (x + lead(x)) / 2,
         dy = (lead(y) - y) / (lead(x) - x),
         x3 = (x2 + lead(x2)) / 2,
         ddy = (lead(dy) - dy) / (lead(x2) - x2)) %>%
  filter(between(x2, 5, 7)) %>%
  slice_min(ddy, n = 1) %>%
  ungroup() 


d %>%
  arrange(x3) %>%
  mutate(samp = fct_inorder(samp)) %>%
  ggplot(aes(x = samp, y = x3, color = samp)) +
  geom_point(size = 5) +
  scale_color_manual(values = clrs) +
  scale_x_discrete(labels = rnm) +
  coord_flip() +
  facet_wrap(~'Point of fastest decline in contact decay rate') +
  ylab('Genomic separation with most\nnegative acceleration of P(s)') +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color= 'black'),
        axis.text = element_text(size = 11, color = 'black'),
        panel.grid = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.y = element_line(color = 'grey70'),
        strip.clip = 'off',
        axis.title.y = element_blank()) -> p
ggsave('sf6_d.pdf', p, height = 5.77, width = 5, bg = 'transparent')
