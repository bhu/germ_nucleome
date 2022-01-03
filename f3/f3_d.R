library(data.table)
library(tidyverse)
library(rasterly)
library(pals)
library(cowplot)
library(patchwork)

load('../data/clust/promoter.rda')

cclrs <- kelly(22)[c(8,10,11,12,15,14,20)]

leg <- tibble(clr = cclrs) %>% 
  mutate(x = as.character(1:7) %>%
           fct_inorder()) %>% 
  {ggplot(., aes( x, x, color = x)) + 
  geom_point() + 
  scale_color_manual(values = cclrs) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_blank())} %>%
  get_legend()

pts <- d %>%
  mutate(lab = cl %>%
           as.character() %>%
           fct_inseq()) %>%
  arrange(lab, x, y) %>%
  as.data.table() %>%
  .[,c('lab','x','y')] %>%
  rasterly(mapping = aes(x = x, y = y, color = lab), color = cclrs) %>%
  rasterly_points()
pts['background']$rasterly_env <- '#ffffff00'
pts['background']$rasterlyPoints1 <- '#ffffff00'
pts <- rasterlyGrob(pts, xlab = 'UMAP 1', ylab = 'UMAP 2', legend = F, axes = F)

wrap_plots(pts, leg, nrow = 1) &
  theme(plot.background = element_blank()) &
  ggsave('f3_d.pdf', height = 4, width = 7.5)

