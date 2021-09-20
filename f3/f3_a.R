library(data.table)
library(tidyverse)
library(rasterly)
library(pals)
library(cowplot)
library(patchwork)

load('../data/clust/open.rda')

cclrs <- c('#bbbbbb', '#377eb8', '#66a61e', '#984ea3',
           '#00d2d5', '#ff7f00', '#af8d00', '#7f80cd', 
           '#b3e900', '#c42e60', '#a65628', '#f781bf',
           '#8dd3c7', '#bebada', '#fb8072', '#80b1d3', 
           '#fdb462', '#fccde5', '#bc80bd', '#ffed6f')

leg <- tibble(clr = cclrs) %>% 
  mutate(x = {c(0, (1:n())[-1]) - 1} %>%
           as.character() %>%
           fct_inorder()) %>% 
  {ggplot(., aes( x, x, color = x)) + 
  geom_point() + 
  scale_color_manual(values = cclrs, guide = guide_legend(nrow = 10)) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key = element_blank())} %>%
  get_legend()

pts <- d %>%
  mutate(lab = case_when(cl > -1 ~ cl + 1,
                         T ~ cl) %>%
           as.character() %>%
           fct_inseq(),
         x = -x) %>%
  arrange(lab, x, y) %>%
  as.data.table() %>%
  .[,c('lab','x','y')] %>%
  rasterly(mapping = aes(x = y, y = x, color = lab), color = cclrs[1:20]) %>%
  rasterly_points()
pts['background']$rasterly_env <- '#ffffff00'
pts['background']$rasterlyPoints1 <- '#ffffff00'
pts <- rasterlyGrob(pts, xlab = 'UMAP 1', ylab = 'UMAP 2', legend = F, axes = F)

wrap_plots(pts, leg, nrow = 1) &
  theme(plot.background = element_blank()) &
  ggsave('f3_a.pdf', height = 4, width = 7.5)

