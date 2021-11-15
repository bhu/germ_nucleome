library(tidyverse)
library(pals)
library(mgcv)
library(gratia)

load('../data/tracks/inpnorm.50kb.rda')

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
samps <- rnm(samps)
clrs <- setNames(tableau20(12)[seq(1, 10, 2)],  samps)
clrs2 <- c(setNames(tableau20(12)[seq(1, 10, 2)], paste('p', samps)),
           setNames(tableau20(12)[seq(2, 10, 2)], paste('t', samps)))

pd <- c('Laminb1', 'PC1') %>%
  setNames(., .) %>%
  lapply(function(m) {
    d50 %>% 
      select(chr, start, samp, all_of(m)) %>%
      mutate(samp = rnm(samp)) %>%
      pivot_wider(names_from = 'samp', values_from = all_of(m)) %>%
      mutate(across(-c(chr, start), ~ GSC - .x)) %>%
      group_by(chr) %>%
      mutate(x = (1:n())/n()) %>%
      ungroup() %>%
      select(all_of(samps), x, chr) %>%
      select(-GSC) %>%
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
  }) %>%
  bind_rows(.id = 'mark')

p <- pd %>% 
  filter(samp == 'EpiLC') %>%
  arrange(mark) %>%
  mutate(ttl = sub('Laminb1', 'lamin B1', mark) %>%
           sub('PC1', 'comp. score', .) %>%
           sprintf('Chromosome-wide \u0394%s', .) %>%
           fct_inorder()) %>%
  group_by(mark) %>%
  mutate(ym = max(abs(est))) %>%
  ungroup() %>%
  mutate(y = est * min(ym) / ym) %>%
  ggplot(aes(x = x, y = est)) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, color = NA,
              fill = 'grey50') +
  geom_line(aes(color = y)) +
  scale_color_gradientn(colors = clrs[c('EpiLC', 'GSC')], limits = c(-.05, .05),
                        oob = scales::squish) +
  facet_grid(ttl ~ 'In vitro', scale = 'free_y') +
  scale_x_continuous(breaks = c(0, 1),
                     labels = c('pter\n(Centromere)',  'qter\n(Telomere)'),
                     expand = expansion(0)) +
  scale_y_continuous('GSC - EpiLC', breaks = scales::pretty_breaks(3)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.clip = 'off',
        legend.text = element_text(size = 11),
        plot.margin = margin(5,35,5,5),
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) 

  
ggsave('f4_h.pdf', p, width = 4, height = 6.77, bg = 'transparent', device = cairo_pdf)

