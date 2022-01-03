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
           sub('b1', ' B1', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x = samp, y = dist, color = mark, group = mark)) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_line(aes(linetype = mark), size = 1, alpha = .5) +
  geom_point(aes(shape = mark), size = 3) +
  coord_flip() +
  scale_color_manual(values = c('firebrick', 'forestgreen')) +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  ylab('Corr. vs H3K36me2\nin specified cell type') +
  guides(color = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2)) +
  theme(plot.background = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        legend.direction = 'horizontal',
        strip.clip = 'off',
        legend.text = element_text(size = 11),
        legend.key = element_blank()) -> p
ggsave(file = 'sf6_b.pdf', p, height = 4.5, width = 2.5)

