library(tidyverse)
library(pals)
library(cowplot)
library(patchwork)
library(scales)
library(ggrepel)
load('../data/contactdecay/differentiation.rda')

samps <- c('ESC (2i)','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
           'ESC (serum)', 'Neural progenitors', 'Cortical neurons',
           'Day 6','Day 4','Day 2','B\u03B1',
           'Cardiac mesoderm', 'Cardiac progenitors',
           'Primitive cardiomyocytes', 'Ventricular cardiomyocytes')

rnm <- function(x) {
  sub('ESC', 'mESC', sub('PGCLC', ' mPGCLC', x))
}

samps <- rnm(samps)

clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(16)[-4]), samps)

dat <- o %>%
  Map(function(x, k) {
    d <- bind_rows(x, .id = 's')
    if (k == 'der') select(d, sep = s_bp, y = slope, s)
    else select(d, sep = s_bp, y = balanced.avg, s)
  }, ., names(.)) %>%
  bind_rows(.id = 'k') %>%
  separate(s, c('study', 'x'), '/') %>%
  mutate(x = case_when(x == 'Ba' ~ 'B\u03B1',
                       x == 'CM' ~ 'Cardiac mesoderm',
                       x == 'CPC' ~ 'Cardiac progenitors',
                       x == 'PCM' ~ 'Primitive cardiomyocytes',
                       x == 'VCM' ~ 'Ventricular cardiomyocytes',
                       x == 'ESC' & study == 'Nagano' ~ 'ESC (2i)',
                       x == 'ESC' & study == 'Bonev2017' ~ 'ESC (serum)',
                       x == 'NPC' ~ 'Neural progenitors',
                       x == 'CN' ~ 'Cortical neurons',
                       T ~ x) %>%
           sub('^D', 'Day ', .) %>%
           rnm() %>%
           factor(samps),
         study = c(Nagano = 'Germline', Bonev2017 = 'Neural',
                   Stadhouders2018 = 'Reprogram',
                   Zhang2019 = 'Cardiac')[study] %>%
           factor(c('Germline', 'Neural', 'Reprogram', 'Cardiac'))) %>%
  na.omit()


dat[dat$k == 'der',] %>%
  ggplot(aes(x = sep, y = y, color = x)) +
  geom_line(alpha = .8) +
  geom_label_repel(aes(label = x), data = ~subset(., x == 'd4c7 mPGCLC' & between(sep, 5e5,5.5e5)),
                   nudge_x = -.5, nudge_y = -.5) +
  scale_color_manual(values = clrs) +
  scale_x_log10(breaks = 10^(4:8),
                labels = c('10kb', '100kb', '1mb', '10mb', '100mb')) +
  facet_wrap(~'Rate of contact probability decay') +
  labs(x = 'Genomic separation, s', y = 'Derivative of P(s)') +
  coord_cartesian(xlim = c(1e4, 1e7),
                  ylim = c(-2.7, .1)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1,),
        panel.grid = element_blank()) -> p

ggsave('f2_h.pdf', height = 2.93, width = 4)
