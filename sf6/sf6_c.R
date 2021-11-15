library(tidyverse)
library(pals)
library(ggrepel)
library(patchwork)

load('../data/contactdecay/germline.smooth.rda')

clrs <- setNames(brewer.set2(6)[-2],
                 c('E11.5 PGC', 'E13.5 mPGC', 'E13.5 fPGC', 'ICM', 'mESC (serum)'))
rnm <- function(x) {
  sub('^ESC', 'mESC', sub('serum', '(serum)', sub('_', ' ', x) )) 
}

d <- o %>%
  pivot_longer(-sep, names_to = 'samp', values_to = 'v') %>%
  filter(sep > 1000) %>%
  filter(grepl('serum|E1|ICM', samp)) %>%
  group_by(samp) %>%
  arrange(sep) %>%
  mutate(x = log10(sep),
         y = log10(v),
         x2 = (x + lead(x)) / 2,
         dy = (lead(y) - y) / (lead(x) - x),
         x3 = (x2 + lead(x2)) / 2,
         ddy = (lead(dy) - dy) / (lead(x2) - x2)) %>%
  ungroup()

p1 <- d %>%
  filter(grepl('serum|E1|ICM', samp)) %>%
  filter(between(x2, 4.5, 7.5)) %>%
  split(., .$samp) %>% 
  lapply(function(x) {
    x %>%
      {smooth.spline(.$x2, .$dy)} %>%
      predict(tibble(x = seq(4.5, 7.5, .01))) %>% 
      {tibble(x = .$x[[1]], y = .$y[[1]])}
  }) %>%
  bind_rows(.id = 'samp') %>%
  mutate(samp = rnm(samp) %>% 
           factor(names(clrs))) %>%
  ggplot(aes(x = x, y = y, color = samp)) +
  geom_line() +
  scale_x_continuous(breaks = 5:7,
                     labels = c('100kb', '1mb', '10mb')) +
  coord_cartesian(xlim = c(5, 7)) +
  scale_color_manual(values = clrs) +
  facet_grid(. ~ 'Contact frequency decay rate') +
  labs(x = 'Genomic separation, s', y = 'Derivative of P(s)') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11))
  
p2 <- d %>%
  filter(between(x2, 5, 7)) %>%
  group_by(samp) %>%
  slice_min(ddy, n = 1) %>%
  ungroup() %>%
  mutate(samp = rnm(samp)) %>%
  ggplot(aes(x = 10^x3, color = samp)) +
  geom_linerange(aes(ymin = 0, ymax = 1)) +
  geom_text_repel(aes(label = samp, x = 10^x3), y = 0, direction = 'x', force_pull = 10,
                  min.segment.length = 0, angle = 90, ylim = c(-1,0)) +
  coord_cartesian(clip = 'off', ylim = c(0, 1), xlim = c(1e5,1e7)) + 
  scale_y_continuous(expand = expansion(0)) +
  scale_x_continuous('Point of fastest decline in\ncontact frequency decay rate',
                     sec.axis = dup_axis(breaks = c(1e5, 1e6, 1e7),
                                         labels = c('100kb', '1mb', '10mb')), trans = 'log10') +
  scale_color_manual(values = clrs) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x.bottom  = element_blank(),
        axis.ticks.x.bottom  = element_blank(),
        axis.title.x.top  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x =  element_text(color = 'black', size = 13, face = 'bold',
                                     margin = margin(t = 80)),
        legend.position = 'none',
        strip.placement = 'outside',
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x.top = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11))

{wrap_plots(p1, p2, ncol=1, heights = c(2,1)) &
    theme(plot.background = element_blank())} %>%
ggsave('sf6_c.pdf', ., height = 6, width = 7.5, bg = 'transparent')

