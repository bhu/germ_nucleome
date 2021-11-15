library(tidyverse)
library(matrixStats)
library(pals)
library(dendextend)
library(ComplexHeatmap)

load("../data/peaks/ATAC/public.fpkm.rda")
samps <- c("ESC", "EpiLC", "d2PGCLC", "d4PGCLC", "d4c7PGCLC", "GSC", "MEF")

mat <- atac.my.pub.pks %>% 
  select(matches('BDF121|MEF_Gia|GSC_AAG')) %>%
  slice_max(rowVars(as.matrix(.)), n = 2e3) %>%
  mutate(idx = 1:n()) %>%
  pivot_longer(-idx, names_to = 'samp', values_to = 'v') %>%
  mutate(samp = sub('_[0-9]$', '', samp)) %>%
  group_by(idx, samp) %>%
  summarise(v = mean(v), .groups = 'drop') %>%
  pivot_wider(names_from = 'samp', values_from = 'v') %>%
  select(-idx) %>%
  {names(.) <- sub('_.*', '', names(.)) ; .} %>%
  .[,samps] %>%
  `colnames<-`(paste0(samps, 's')) %>%
  t()

dend <- dist(t(mat)) %>% 
  hclust() %>%
  as.dendrogram() %>%
  color_branches(7, col = kelly(8)[-1])

pdf('sf5_c.pdf', width = 12, height = 2.65, bg = 'transparent')
ht <- Heatmap(mat, show_column_names = F, cluster_rows = F, col = viridis(100),
              use_raster = T, raster_quality = 20, cluster_columns = dend, column_split = 7,
              heatmap_legend_param = list(title = 'log2(FPKM + 1)', direction = "horizontal",
                                          title_position = "topcenter", title_gp = gpar(fontface = 'plain')))
draw(ht, background = "transparent", padding = unit(c(2, 2, 2, 10), "mm"))
dev.off()
