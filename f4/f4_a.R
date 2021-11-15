library(tidyverse)
library(pals)

load('../data/tracks/inpnorm.50kb.rda')
smps <- c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC')
mrks <- c('ATAC','CTCF','Rad21',
          'Ring1b','H2Aub','K27me3', 'K27ac', 'K4me1','K4me3',
          'K36me2','K36me3','K9me2','K9me3', 'Laminb1')

d50 %>%
  split(., .$samp) %>%
  lapply(function(x) {
    x[,-21:-18] %>%
      cor(method = 'spearman') %>%
      {.['PC1',]} %>%
      data.frame(dist = .) %>%
      rownames_to_column('mark')
  }) %>%
  bind_rows(.id = 'samp') %>%
  filter(mark != 'PC1') %>%
  mutate(samp = factor(samp, smps),
         mark = factor(mark, mrks)) %>%
  na.omit() %>%
  arrange(samp, mark) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           fct_inorder()) %>%
  ggplot(aes(x = samp, y = mark, fill = dist)) + 
  geom_tile() +
  #geom_text(aes(label = round(dist, 2))) +
  scale_fill_gradientn(expression('Spearman\'s' ~ rho), breaks = c(-1, 0, 1),
                       colors = coolwarm(25), limits = c(-1, 1)) +
  coord_cartesian(expand = F) +
  facet_grid(.~'Comp. score corr.') +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.justification = c(1,2.03),
        plot.margin = margin(5,5,45,5),
        legend.direction = 'horizontal',
        legend.position = c(1,0)) +
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, title.hjust = .5,
                               title.position = 'top',
                               frame.linewidth = 1, ticks = F)) -> p
ggsave(file = 'f4_a.pdf', p, height = 4.55, width = 2.6)
