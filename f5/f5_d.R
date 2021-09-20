library(tidyverse)
library(ppcor)
library(pals)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')

load('../data/tracks/inpnorm.50kb.rda')
d50 %>%
  dplyr::select(-c(chr, start, end, PC1)) %>%
  split(., .$samp) %>%
  lapply(function(x) {
    x %>%
      dplyr::select(-samp) %>%
      {tibble(mark = names(.),
              partial = pcor(., method = 'pearson')$estimate[which(names(.) == 'Laminb1'),],
              full = cor(., method = 'pearson')['Laminb1',])}
  }) %>%
  bind_rows(.id = 'samp') %>%
  mutate(samp = factor(samp, samps)) %>%
  na.omit() %>%
  filter(mark %in% c('K9me3', 'K9me2')) %>%
  mutate(mark = paste0('H3', mark)) %>%
  ggplot(aes(x = samp, y = mark, fill = full)) +
  geom_tile() +
  geom_text(aes(label = round(full, 2))) +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  scale_fill_gradientn(expression('Spearman\'s'~rho),
                       colors = coolwarm(25), limits = c(-1,1),
                       breaks = c(-1,0,1)) +
  #guides(fill = guide_colorbar(barwidth = 4.7, barheight = .5, title.position = 'top',
  #                             frame.linewidth = 1, ticks = F)) +
  guides(fill = guide_colorbar(barheight = 2, barwidth = .5, title.position = 'top', ticks = F)) +
  coord_cartesian(expand = F) +
  facet_grid(.~'Lamin B1 correlation') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11),
        #legend.position = c(1,1),
        #legend.direction = 'horizontal',
        #legend.justification = c(1,-.3),
        #plot.margin = margin(50,2,2,2),
        strip.background = element_rect(fill = NA),
        legend.title = element_text(angle = 270, vjust = 0),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave('f5_d.pdf', width = 3.5, height = 2.5)
