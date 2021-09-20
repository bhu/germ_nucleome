library(data.table)
library(tidyverse)
library(ggtext)
library(rasterly)
library(cowplot)
library(patchwork)
library(pals)
library(ggnewscale)


cclrs <- kelly(22)[c(8,10,11,12,15,14,20)]
odr <- c('Expression', 'CTCF', 'Rad21',  'ATAC', 'K4me1', 'K4me3', 'K27ac',
         'K27me3', 'Ring1b', 'H2Aub', 'K36me2', 'K36me3', 'K9me2', 'K9me3') %>%
  sub('^K', 'H3K', .)

load('../data/clust/promoter.rda')

pd <- d %>%
  filter(cl != -1 & type != 'GSCLC') %>%
  {.[.$gene %in% Reduce(intersect, split(.$gene, .$type)),]} %>%
  select(clu = cl, Expression = SC3, CTCF, Rad21,  ATAC, K4me1, K4me3, K27ac,
         K27me3, Ring1b, H2Aub, K36me2, K36me3, K9me2, K9me3) %>%
  pivot_longer(-clu, names_to = 'mark', values_to = 'v') %>%
  mutate(mark = case_when(grepl('^K', mark)~ paste0('H3', mark), T ~ mark) %>%
           fct_inorder()) %>%
  group_by(mark, clu) %>%
  summarise(v = mean(v, na.rm = T)) %>%
  ungroup() %>%
  mutate(clu = as.character(clu),
         mark = factor(mark, odr),
         y = as.numeric(mark)) 
ggplot(pd, aes(x = clu, y = y)) +
  geom_tile(aes(fill = v), data = ~subset(., mark=='Expression')) +
  scale_fill_gradientn('log2(RPM)', colors = brewer.blues(25), breaks = c(6, 7),
                       guide = guide_colorbar(barheight = .5, barwidth = 4, title.vjust = 1, order = 1)) +
  new_scale_fill() +
  geom_tile(aes(fill = v), data = ~subset(., mark !='Expression')) +
  scale_fill_gradientn('Enrichment', colors = rev(brewer.brbg(25)), 
                       limits = c(-3,3), breaks = c(-3,0,3), oob = scales::squish,
                       guide = guide_colorbar(barheight = .5, barwidth = 4, title.vjust = 1, order = 2)) +
  geom_hline(yintercept = 1.5) +
  scale_x_discrete(labels = sprintf('<span style=\'color:%s\'>%s</span>', cclrs, 1:7)) +
  scale_y_continuous(breaks = distinct(pd, y, mark)$y, labels = distinct(pd, y, mark)$mark) +
  labs(x = 'Promoter cluster', y = 'Mark') +
  coord_flip(expand = F) +
  facet_grid(.~ 'Promoter clusters\' characteristics') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(color = 'black', angle = 30, hjust = 1),
        legend.position = 'bottom',
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        #legend.justification = c(1,2.3),
        legend.background = element_blank(),
        legend.direction = 'horizontal',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text.y = element_markdown(),
        axis.title.x = element_blank(),
        plot.margin = margin(2,2,10,15)) +
  ggsave('f3_e.pdf', height = 3, width = 4.25)

