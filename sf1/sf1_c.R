library(tidyverse)
library(pals)
library(cowplot)
library(patchwork)
library(scales)
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


p1 <- dat[dat$k == 'log',] %>%
  ggplot(aes(x = sep, y = y)) +
  geom_line(aes(color = x), alpha = .8) +
  scale_color_manual(values = clrs) +
  scale_x_log10(breaks = 10^(4:8),
                labels = c('10kb', '100kb', '1mb', '10mb', '100mb')) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(.~ study) +
  labs(x = 'Genomic separation, s', y = 'P(s)') +
  coord_cartesian(xlim = c(1e4, 1e8),
                  ylim = c(1e-7,1e-3)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) 

p2 <- dat[dat$k == 'log',] %>%
  ggplot(aes(x = sep, y = y)) +
  geom_line(aes(color = x), alpha = .8) +
  scale_color_manual(values = clrs) +
  scale_x_log10(breaks = 10^(4:8),
                labels = c('10kb', '100kb', '1mb', '10mb', '100mb')) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(. ~ study) +
  labs(x = 'Genomic separation, s', y = 'P(s)') +
  coord_cartesian(xlim = c(1e7, 1e8),
                  ylim = c(1e-7,1e-5)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        panel.grid = element_blank()) 

{wrap_plots(p1,p2,nrow=2) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf1_c.pdf', ., height = 4.3, width = 9)


leg <- ggplot(dat, aes(x = sep, y = y, color = x)) +
  geom_line() +
  scale_color_manual(values = clrs) +
  guides(color = guide_legend(ncol = 2)) +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.justification = 'top',
        legend.text = element_text(size = 11)) 

ggsave('sf1_c_leg.pdf', get_legend(leg), height = 4.3, width = 4, device = cairo_pdf, bg = 'transparent')

