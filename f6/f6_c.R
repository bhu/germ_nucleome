library(tidyverse)
library(pals)

load('../data/tracks/inpnorm.50kb.rda')
smps <- c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC')
mrks <- c('ATAC','CTCF','Rad21','Stag1','Stag2',
          'Ring1b','H2Aub','K27me3', 'K27ac', 'K4me1','K4me3',
          'K36me2','K36me3','K9me2','K9me3', 'Laminb1', 'PC1')

d50 %>%
  split(., .$samp) %>%
  lapply(function(x) {
    x[,-21:-18] %>%
      cor(method = 'spearman') %>%
      {.['K36me2',]} %>%
      data.frame(dist = .) %>%
      rownames_to_column('mark')
  }) %>%
  bind_rows(.id = 'samp') %>%
  filter(mark %in% c('PC1', 'Laminb1')) %>%
  mutate(samp = factor(samp, rev(smps)),
         mark = factor(mark, rev(mrks))) %>%
  na.omit() %>%
  arrange(samp, mark) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           sub('PC1', 'Compartment score', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x = mark, y = samp, fill = dist)) + 
  geom_tile() +
  geom_text(aes(label = round(dist, 2))) +
  scale_fill_gradientn(expression('Spearman\'s' ~ rho), breaks = c(-1, 0, 1),
                       colors = coolwarm(25), limits = c(-1, 1)) +
  coord_cartesian(expand = F, clip = F) +
  scale_y_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  scale_x_discrete(labels = function(x) sub('b1', ' B1', x)) +
  facet_grid(.~'H3K36me2 corr.') +
  theme(plot.background = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.position = c(1,0),
        legend.direction = 'horizontal',
        strip.clip = 'off',
        legend.justification = c(1,2),
        plot.margin = margin(5,5,30,5),
        legend.key = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 4.7, barheight = 0.5, 
                               frame.linewidth = 1, ticks = F, title.position = 'top')) +
  ggsave(file = 'f6_c.pdf', height = 4.5, width = 2.7)

