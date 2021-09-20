library(data.table)
library(tidyverse)
library(pals)

tot <- fread('../data/resources/mm10.chrom.sizes') %>%
  filter(V1 %in% paste0('chr', 1:19)) %>%
  pull(V2) %>% 
  sum()

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

list.files('../data/peaks/LADs', full.names = T)  %>% 
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(fread) %>% 
  bind_rows(.id = 'samp') %>% 
  mutate(w = V3 - V2) %>% 
  filter(V1 %in% paste0('chr', 1:19)) %>% 
  group_by(samp) %>%
  summarise(w = sum(w)) %>%
  mutate(f = 100 * w / tot) %>%
  mutate(study = case_when(grepl('Poleshko', samp) ~ 'Poleshko 2017',
                           grepl('BDF121|AAG', samp) ~ 'This study',
                           grepl('OPC', samp) ~ 'Yattah 2020',
                           grepl('diffe', samp) ~ 'Robson 2016',
                           T ~ 'Peric-Hupkes 2010') %>%
           factor(rev(c('This study', 'Yattah 2020', 'Poleshko 2017', 'Robson 2016', 'Peric-Hupkes 2010'))),
         samp = case_when(samp == 'NIH3T3' ~ 'MEF',
                          samp == 'OPC_PF' ~ 'OPC',
                          samp == 'OPC_T3' ~ 'OL',
                          samp == 'C2C12_differentiated' ~ 'Myotube',
                          samp == 'C2C12_undifferentiated' ~ 'Myoblast',
                          T ~ sub('_.*', '', samp)) %>%
           factor(rev(c(samps, 'NPC', 'Astrocyte', 'MEF', 'CM', 'OPC', 'OL', 'Myoblast', 'Myotube')))) %>%
  na.omit() %>%
  arrange(study, samp) %>%
  mutate(study = sub('\\ [0-9]*$', '', study) ,
         study = case_when(study == 'This study' ~ study,
                           study == 'Peric-Hupkes' ~ 'P-H',
                           T ~ substr(study, 1, 1)) %>%
           fct_inorder(),
         clr = case_when(study == 'This study' ~ clrs[as.character(samp)],
                         T ~ 'grey30'))  %>%
  ggplot(aes(y = samp, x = f, fill = clr)) +
  geom_col(alpha = .8) +
  facet_grid(study ~., space = 'free_y', scale = 'free_y') +
  scale_x_continuous('% of genome in LADs', expand = expansion(c(0, .05)),
                     breaks = c(0,25,50)) +
  scale_fill_identity() +
  scale_y_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.line.y = element_line(color = 'black'),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed')) +
  ggsave('f4_d.pdf', height = 4.6, width = 2.4)
