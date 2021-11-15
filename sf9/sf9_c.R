library(tidyverse)
library(pals)
library(mgcv)
library(gratia)
library(ggnewscale)
library(patchwork)
library(ggdist)

load('../data/tracks/inpnorm.50kb.rda')

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
samps <- rnm(samps)
clrs <- setNames(tableau20(12)[seq(1, 10, 2)],  samps)
clrs2 <- c(setNames(tableau20(12)[seq(1, 10, 2)], paste('p', samps)),
           setNames(tableau20(12)[seq(2, 10, 2)], paste('t', samps)))

pd <- d50 %>% 
  select(chr, start, samp, Laminb1) %>%
  mutate(samp = rnm(samp)) %>%
  pivot_wider(names_from = 'samp', values_from = 'Laminb1') %>%
  #filter(chr %in% c('chr1', 'chr2')) %>%
  group_by(chr) %>%
  mutate(x = (1:n())/n()) %>%
  ungroup() %>%
  select(all_of(samps), x, chr) %>%
  pivot_longer(-c(x, chr), names_to = 'samp', values_to = 'y') %>%
  split(., .$samp) %>%
  lapply(function(d) {
    m <- gam(y ~ s(x, k = 10), data = d, method = 'REML')
    confint(m, parm = 's(x)', type = 'confidence',  level = .9999,
            newdata = tibble(x = seq(0, 1, .001))) %>%
      mutate(x = x$x)
  }) %>%
  bind_rows(.id = 'samp') %>%
  mutate(samp = factor(samp, samps))

ggplot(pd, aes(x = x, y = est, color = samp, fill = samp)) + 
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5, color = NA) +
  scale_fill_manual(values = clrs) +
  scale_color_manual(values = clrs) +
  labs(x = 'Position along chromosomes',
       y = 'log2(ChIP / input)') +
  scale_x_continuous(breaks = c(0, .5, 1),
                     labels = c('pter', '', 'qter'),
                     expand = expansion(0)) +
  facet_grid(.~'Distribution of lamin B1 across chromosomes') +
  theme(legend.position = 'bottom',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.margin = margin(5,15,5,5),
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) -> p

ggsave('sf7_c.pdf', p, width = 6, height = 3)

