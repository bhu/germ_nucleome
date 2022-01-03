library(tidyverse)
library(pals)
library(patchwork)
load('../data/contactdecay.rda')

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)

dat <- o %>%
  Map(function(x, k) {
    d <- bind_rows(x, .id = 's')
    if (k == 'der') select(d, sep = s_bp, y = slope, s)
    else select(d, sep = s_bp, y = balanced.avg, s)
  }, ., names(.)) %>%
  bind_rows(.id = 'k') %>%
  separate(s, c('study', 'x'), '/') %>%
  filter(study %in% 'Nagano' & x %in% samps) %>%
  mutate(x = factor(x, samps))

fg <- range(dat[dat$k == 'log',]$y)
bg <- range(dat[dat$k == 'der',]$y)
b <- diff(fg)/diff(bg)
a <- fg[1] - b * bg[1]


p1 <- dat[dat$k == 'log',] %>%
  ggplot(aes(x = sep, y = y)) +
  geom_line(aes(color = x, linetype = x), alpha = .8) +
  scale_color_manual(values = clrs) +
  scale_y_log10() +
  facet_grid(. ~ 'Contact probability decay') +
  labs(x = 'Genomic separation, s', y = 'P(s)') +
  coord_cartesian(xlim = c(1e4, 1e8),
                  ylim = c(1e-7,1e-3)) +
  scale_x_log10(breaks = 10^(4:8),
                labels = c('10kb', '100kb', '1mb', '10mb', '100mb')) +
  theme(legend.position = c(0.01,0.01),
        legend.justification = c(0,0),
        legend.background = element_rect(fill = '#ffffff66', color = NA),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major = element_line(color = 'grey85', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) 

p2 <- dat[dat$k == 'der',] %>%
  ggplot(aes(x = sep, y = y)) +
  geom_line(aes(color = x, linetype = x), alpha = .8) +
  scale_color_manual(values = clrs) +
  labs(x = 'Genomic separation, s', y = 'Derivative of P(s)') +
  coord_cartesian(xlim = c(1e4, 1e8),
                  ylim = c(-2.7, .1)) +
  scale_x_log10(breaks = 10^(4:8),
                labels = c('10kb', '100kb', '1mb', '10mb', '100mb')) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        panel.grid.major = element_line(color = 'grey85', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) 

{ wrap_plots(p1, p2, nrow = 2) &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf7_b.pdf', ., height = 3.8, width = 3.3)
