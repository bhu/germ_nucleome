library(data.table)
library(tidyverse)
library(pals)

samps <- c('ESC (2i)','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
           'ESC (serum)', 'Neural progenitors (NPC)', 'Cortical neurons (CN)',
           sprintf('Day %d with OSKM', c(6,4,2)),'B cell with C/EBP\u03B1 (B\u03B1)',
           'Cardiac mesoderm (CM)', 'Cardiac progenitors (CPC)',
           'Primitive cardiomyocytes (PCM)', 'Ventricular cardiomyocytes (VCM)')
rnm <- function(x) {
  #x %>%
  #  sub('C ', 'Cs ', .) %>%
  #  sub('C$', 'Cs', .) %>%
  #  sub("^(.*[A-Z])\\)", '\\1s)', .)
  x <- sub('ESC', 'mESC', x)
  x <- case_when(x == 'ESC' ~ 'mESC',
                 x == 'd2PGCLC' ~ 'd2 mPGCLC',
                 x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
                 T ~ x)
}
samps <- rnm(samps)

clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(12)[-4]), samps)

load('../data/compscore/100kb.rda')
dat.ratio <- c(mats$mm10, mats$hg38) %>%
  lapply(function(m) {
    na.omit(m) %>%
      pivot_longer(everything(), names_to = 'x', values_to = 'y') %>% 
      dplyr::filter(!grepl('_', x)) 
  }) %>%
  rbindlist(idcol = 'study') %>%
  mutate(x = case_when(x == 'Ba' ~ 'B cell with C/EBP\u03B1 (B\u03B1)',
                       x == 'CM' ~ 'Cardiac mesoderm (CM)',
                       x == 'CPC' ~ 'Cardiac progenitors (CPC)',
                       x == 'PCM' ~ 'Primitive cardiomyocytes (PCM)',
                       x == 'VCM' ~ 'Ventricular cardiomyocytes (VCM)',
                       x == 'ESC' & study == 'Nagano' ~ 'ESC (2i)',
                       x == 'ESC' & study == 'Bonev2017' ~ 'ESC (serum)',
                       x == 'NPC' ~ 'Neural progenitors (NPC)',
                       x == 'CN' ~ 'Cortical neurons (CN)',
                       T ~ x) %>%
           sub('^D([0-9])', 'Day \\1 with OSKM', .) %>%
           rnm() %>%
           factor(samps),
         study = c(Nagano = 'Germline', Bonev2017 = 'Neural',
                   Stadhouders2018 = 'Reprogram',
                   Zhang2019 = 'Cardiac')[study] %>%
           factor(c('Germline', 'Neural', 'Reprogram', 'Cardiac'))) %>%
  na.omit() %>%
  group_by(x, study) %>%
  summarise(A = 100 * sum(y > 0) / n(),
            B = 100 * sum(y < 0) / n()) %>%
  ungroup() 




dat.ratio <- dat.score %>%
  group_by(x, study) %>%
  summarise(y = 100 * sum(y > 0) / n()) %>%
  ungroup()

m <- 0.2
b <- -10

sclr <- 'gray60'
ggplot(dat.score, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_violin(aes(fill = x), alpha = .7, color = NA) +
  geom_line(aes(y = m * y + b, group = study), data = dat.ratio, 
            size = 1.5, color = sclr, alpha = .7) +
  geom_point(aes(y = m * y + b, fill = x), data = dat.ratio, pch = 21, 
             color = alpha(sclr, .7), size = 4, stroke = 2) +
  facet_grid(.~study, scales = 'free_x', space = 'free_x') +
  coord_cartesian(ylim = c(-2.5, 1.5)) +
  scale_y_continuous('Compartment score',
                     sec.axis = sec_axis(~ (. - b) / m, '% of bins in compartment A')) +
  scale_fill_manual(values = clrs) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y.right = element_text(color = sclr),
        axis.text.y.right = element_text(color = sclr),
        panel.grid = element_blank()) -> p



ggsave('f1_g.pdf', p, width = 5.8, height = 4, device = cairo_pdf, bg = 'transparent')
