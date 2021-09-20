library(tidyverse)
library(gplots)
library(pals)
library(dendextend)

samps <- c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC", "GSCLC")
clrs <- setNames(c(tableau20()[seq(1, 9, 2)], '#499894'), samps)

load('../data/hicrep.rda')

pdf('sf1_b.pdf', bg = 'transparent', width = 3.4, height = 2.7)
m <- d %>%
  pivot_wider(names_from = 's1', values_from = 'v') %>%
  column_to_rownames('s2') %>%
  as.matrix()
colnames(m) <- rownames(m) <- sub('_', ' ', colnames(m)) %>%
  sub('ESC', 'mESC', .) %>%
  sub('PGCLC', ' mPGCLC', .)
h <- sqrt(1 - m) %>% 
  as.dist() %>% 
  hclust()
dend <- h %>%
  as.dendrogram() %>%
  color_branches(6, col = clrs[unique(sub('[\\ ]*m', '', sub('\\ [12]$', '', h$labels[h$order])))])
heatmap.2(m, trace = "none", col = viridis(10), keysize = 2,
          RowSideColors = clrs[sub('[\\ ]*m', '', sub('\\ [12]$', '', colnames(m)))],
          ColSideColors = clrs[sub('[\\ ]*m', '', sub('\\ [12]$', '', colnames(m)))], key=F,
           Rowv = dend, Colv = dend, labCol=F,lhei = c(1,4),lwid = c(1,4),
          margins = c(3, 8),  key.xlab = "SCC",key.title = NA)
dev.off()

pdf('sf1_b_key.pdf', bg = 'transparent', width = 4, height = 2.7)
heatmap.2(m, trace = "none", col = viridis(10), keysize = 3,
          RowSideColors = clrs[sub('_.*', '', rownames(m))],
          ColSideColors = clrs[sub('_.*', '', colnames(m))], 
          Rowv = dend, Colv = dend, labCol=F,
          margins = c(3, 8),  key.xlab = "SCC",key.title = NA)
dev.off()

