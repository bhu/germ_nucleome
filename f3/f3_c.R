library(tidyverse)
library(scattermore)
library(pals)

load('../data/clust/promoter.rda')

cclrs <- kelly(22)[c(8,10,11,12,15,14,20)]

d %>%
  mutate(lab = c('Repressed', 'Polycomb', 'Active', 'Bivalent', 'CTCF', 'Readthrough', 'Enhancer')[cl]) %>%
  arrange(cl) %>%
  mutate(lab = fct_inorder(lab)) %>%
  ggplot(aes(x, y, color = lab)) +
  geom_scattermore() +
  scale_color_manual('Cluster', values = cclrs, guide = guide_legend(ncol = 1)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_wrap(~'Epigenetic embedding of promoters') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank()) -> p

ggsave('f3_c.pdf', p, height = 3, width = 5.5)
