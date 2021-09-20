library(tidyverse)
library(readxl)
library(pals)
library(MASS)
library(patchwork)
library(gplots)
library(dendextend)
library(plotly)

anchors <- c('ESC', 'EpiLC','d4c7PGCLC', 'GSC', 'GSCLC')
types <- c('ESC', 'EpiLC', 'd2PGCLC','d4c7PGCLC', 'GSC', 'MEF')

rnm <- function(x) {
  sub('PGCLC', ' mPGCLLC', sub('ESC', 'mESC', x))
}
types <- rnm(types)

clrs <- setNames(tableau20(20)[c(1, 3, 5, 7, 9, 17)], types)
mods <- c("K4me1", "K4me3", "K9me2", "K9me3", "K27me3", "K27ac", "K36me2", "K36me3")

d <- c('old', 'new') %>%
  setNames(., .) %>%
  lapply(function(s) {
    read_excel('../data/histone_ratios.xlsx', sheet = s) %>%
      {.[,colSums(is.na(.)) != nrow(.)]} %>% 
      na.omit() %>% 
      {.[,c(1, which(.[1,] == 'Area'))]} %>%
      tail(-1) %>%
      rename(p = 1) %>%
      `names<-`(., sub('\\..*', '', names(.))) %>%
      separate(p, c('pep', 'mod'), '\\ ') %>%
      mutate(pep = sub('H33', 'H3', pep)) %>%
      group_by(pep) %>%
      mutate(across(-mod, function(x) as.numeric(x) / sum(as.numeric(x)))) %>%
      ungroup() %>% 
      filter(grepl('^H[34]', pep) & mod != 'unmod') %>%
      mutate(mod = strsplit(mod, '(?<=.)(?=K)', perl = T)) %>% 
      unnest(mod) %>%
      mutate(mod = paste(substr(pep, 1, 2), mod)) %>%
      dplyr::select(-pep) %>%
      group_by(mod) %>%
      summarise(across(everything(), sum), .groups = 'drop') %>%
      pivot_longer(-mod, names_to = 'samp', values_to = 'a') %>%
      separate(samp, c('type', 'rep'), '_', F) 
  })

m <- lapply(d, function(x) {
  x %>%
    filter(type %in% anchors) %>%
    group_by(mod, type) %>%
    summarise(a = median(a), .groups = 'drop')
}) %>%
  bind_rows(.id = 'batch') %>%
  pivot_wider(names_from = 'batch', values_from = 'a') %>%
  split(., .$mod) %>% 
  lapply(function(x) {
    fit <- tryCatch(
      rlm(old ~ new, x)$coefficients,
      error = function(msg) {
        lm(old ~ new, x)$coefficients
      },
      warning = function(msg) {
        lm(old ~ new, x)$coefficients
      }
    )
    d$new %>%
      filter(mod == x$mod[1]) %>%
      mutate(a = fit[1] + a * fit[2])
  }) %>%
  bind_rows() %>%
  mutate(batch = 'new') %>%
  rbind(mutate(d$old, batch = 'old')) %>%
  separate(mod, c('his', 'mod'), '\\ ') %>%
  mutate(type = factor(rnm(type), types)) %>%
  na.omit() %>%
  filter(his == 'H3') %>%
  arrange(mod, type) %>%
  filter(type != 'unmod') %>%
  group_by(mod) %>% 
  mutate(z = (a - mean(a)) / sd(a)) %>%
  ungroup() %>% 
  dplyr::select(mod, samp, z) %>%
  mutate(samp = sub('pool', '4', samp)) %>%
  pivot_wider(names_from = 'samp', values_from = 'z') %>%
  na.omit() %>%
  column_to_rownames('mod') %>%
  as.matrix()
colnames(m) <- rnm(colnames(m)) %>%
  sub('_', ' ', .)
rownames(m) <- paste0('H3', rownames(m))

h <- cor(m) %>%
  {sqrt(1 - .)} %>%
  as.dist() %>%
  hclust(method = 'ward.D2')
dend <- h %>%
  as.dendrogram() %>%
  color_branches(6, col = clrs[unique(sub('\\ [0-9]$*', '', h$labels[h$order]))])

pdf('f2_b_hm_leg.pdf', bg = 'transparent', width = 5.5, height = 4)
heatmap.2(t(m), trace = "none", col = rev(brewer.puor(100)), Rowv = dend, 
          key.xlab = 'Z-score', key.title = NA, srtCol = 45,
          RowSideColors = clrs[sub('\\ [0-9]$*', '', colnames(m))],
          margins = c(4, 7))
dev.off()

pdf('f2_b_hm2.pdf', bg = 'transparent', width = 4.3, height = 2.8)
heatmap.2(t(m), trace = "none", col = rev(brewer.puor(100)), Rowv = dend, 
          key = F, lwid = c(1,10), lhei = c(1,10),
          RowSideColors = clrs[sub('\\ [0-9]$*', '', colnames(m))],
          margins = c(4, 7))
dev.off()


