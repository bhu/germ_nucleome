library(tidyverse)
library(pals)
library(mgcv)
library(gratia)
library(rtracklayer)
library(ggh4x)

bands <- getTable(ucscTableQuery("rheMac10", table = "cytoBandIdeo")) %>%
  filter(gieStain == 'acen') %>%
  group_by(chrom) %>% 
  summarise(cstart = min(chromStart),
            cend = max(chromEnd)) %>%
  ungroup()


load('../data/compscore/vivo.100kb.rda')

pd <- e %>%
  Map(function(x, nm) {
    o <- x[[1]][,1:3] %>%
      mutate(score = x[[grep('sp.*a', names(x), ignore.case = T)]]$E1 -
               x[[grep('fib', names(x), ignore.case = T)]]$E1) %>%
      na.omit()
    
    o <- if (nm == 'Vara2019') {
      o %>%
        group_by(chrom) %>%
        arrange(start) %>%
        mutate(x = (1:n())/n()) %>%
        ungroup() %>%
        select(chrom, x, y = score) %>%
        gam(y ~ s(x, k = 10), data = ., method = 'REML') %>%
        confint(parm = 's(x)', type = 'confidence',  level = .9999,
                newdata = tibble(x = seq(0, 1, .001))) %>%
        mutate(x = x$x)
    } else {
      o %>%
        merge(bands) %>%
        mutate(arm = case_when(end < cstart ~ 'p',
                               start > cend ~ 'q')) %>%
        na.omit() %>%
        group_by(chrom, arm) %>%
        arrange(start) %>%
        mutate(x = (1:n())/n()/2,
               x = case_when(arm == 'q' ~ x + .5, T ~ x)) %>%
        ungroup() %>%
        select(chrom, x, y = score) %>%
        gam(y ~ s(x, k = 10), data = ., method = 'REML') %>%
        confint(parm = 's(x)', type = 'confidence',  level = .9999,
                newdata = tibble(x = seq(0, 1, .001))) %>%
        mutate(x = x$x)
    }
    o
  }, ., names(.)) %>%
  bind_rows(.id = 'study') %>%
  mutate(species = c(Wang2019 = 'Rhesus macaque', Vara2019 = 'House mouse')[study] %>%
           factor(c('Rhesus macaque', 'House mouse'))) 


p <- ggplot(pd, aes(x = x, y = est)) +
  geom_hline(aes(yintercept = y), data = mutate(distinct(pd, species), y = 0)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, color = NA,
              fill = 'grey50') +
  geom_line() + 
  scale_y_continuous('Spermatogonia - fibroblast', breaks = scales::pretty_breaks(4)) +
  facet_nested('\u0394Compartment score' + species ~ 'In vivo', scales = 'free', independent = 'x') +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = c(0, .5, 1),
                       labels = c('pter\n(Telomere)', 'Centromere',  'qter\n(Telomere)'),
                       expand = expansion(0)),
      scale_x_continuous(breaks = c(0, 1),
                         labels = c('pter\n(Centromere)',  'qter\n(Telomere)'),
                         expand = expansion(0))
  )) +
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
        strip.text = element_text(color = 'black', size = 13, face = 'bold'))

ggsave('sf5_b.pdf', p, width = 4, height = 6.77, bg = 'transparent', device = cairo_pdf)

