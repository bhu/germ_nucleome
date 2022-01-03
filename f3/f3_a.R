library(tidyverse)
library(scattermore)

load('../data/clust/open.rda')

cclrs <- c('#377eb8', '#66a61e', '#984ea3','#00d2d5', 
           '#ff7f00', '#af8d00', '#7f80cd', '#b3e900', 
           '#c42e60', '#a65628', '#f781bf', '#8dd3c7', 
           '#bebada', '#fb8072', '#80b1d3', '#fdb462', 
           '#fccde5', '#bc80bd', '#ffed6f', '#bbbbbb')
odr <- c(16, 17, 19, 1, 3, 4, 5, 6, 15, 18,
         8, 9, 10, 13, 2, 11, 12, 7, 14, 'NA')

d %>%
  mutate(lab = case_when(cl > -1 ~ cl + 1,
                         T ~ cl) %>%
           as.character() %>%
           factor(odr))  %>%
  ggplot(aes(x, y, color = lab)) +
  geom_scattermore() +
  scale_color_manual('Cluster', values = cclrs, guide = guide_legend(ncol = 2)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_wrap(~'Epigenetic embedding of open sites') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.clip = 'off',
        legend.key = element_blank()) -> p

ggsave('f3_a.pdf', p, height = 3, width = 3.5)

