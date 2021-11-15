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



p <- ggplot(dat.ratio, aes(x = x, y = A, color = x)) +
  geom_hline(yintercept = 50, color = 'grey60') +
  geom_line(aes(group = study), color = 'grey50', size = 2, alpha = .5) +
  geom_point(size = 4) +
  facet_grid(.~study, scales = 'free_x', space = 'free_x')  +
  scale_color_manual(values = clrs) +
  scale_y_continuous('% of genome in compartment A') +
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
        panel.grid = element_blank())
ggsave('f1_g.pdf', p, width = 5.5, height = 4, device = cairo_pdf, bg = 'transparent')
