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

samps <- c('ESC', 'EpiLC', 'd2PGCLC','d4c7PGCLC', 'GSC', 'MEF') %>% rnm
clrss <- setNames(tableau20(20)[c(1, 3, 5, 7, 9, 17)], samps)

res <- m %>%
  as.data.frame() %>%
  rownames_to_column('mod') %>%
  pivot_longer(-mod, names_to = 'samp', values_to = 'z') %>%
  mutate(samp = sub('\\ [0-9]$', '', samp)) %>%
  group_by(samp, mod) %>%
  summarise(z = mean(z), .groups = 'drop') %>%
  pivot_wider(names_from = 'samp', values_from = 'z') %>%
  column_to_rownames('mod') %>%
  t() %>%
  prcomp(center = F)

comps <- summary(res)$importance[2,] * 100

axs <- paste0('PC', 1:3)


pts <- res$x %>%
  as.data.frame() %>%
  rownames_to_column('samp') %>%
  mutate(samp = factor(samp, samps)) %>%
  arrange(samp)


scn <- setNames(axs, c('xaxis', 'yaxis', 'zaxis')) %>%
  lapply(function(ax) list(title = sprintf('%s (%d%%)', ax, round(comps[ax])),
                           showticklabels = F)) %>%
  c(list(camera = list(eye = list(x = 1.786, y = 1.938, z = 1.069),
                       up = list(x = 0, y = 0, z = 1),
                       center = list(x = 0.174, y = -.028, z = -.517))))

p <- plot_ly(source = 'src') %>%
  add_markers(type = 'scatter3d', x = pts[[axs[1]]],
              y = pts[[axs[2]]], z = pts[[axs[3]]],
              color = pts$samp, colors = clrss,
              marker = list(size = 20))

lns <- lapply(1:(nrow(pts) - 1), function(i) {
  pts[c(i, i + 1), axs] %>%
    lapply(function(x) {seq(x[1], x[2], length.out = 100)}) %>%
    bind_cols() %>%
    `colnames<-`(c('x', 'y', 'z')) %>%
    mutate(idx = 1:n())
})
mz <- min(pts[[axs[3]]])*1.05
cns <- lapply(seq_along(pts), function(i) {
  pts[i, axs] %>%
    `colnames<-`(c('x', 'y', 'z')) %>%
    {rbind(., mutate(., z = mz))} %>%
    mutate(clr = clrs[pts$samp[i]])
})

for (i in 1:4) {
  p <- p %>%
    add_trace(mode = 'lines', type = 'scatter3d', data = lns[[i]],
              x = ~x, y = ~y, z = ~z,
              line = list(width = 10, color = lns[[i]]$idx,
                          colorscale = list(c(0, clrs[i]),
                                            c(1, clrs[i + 1]))),
              showlegend = FALSE)
}
for (i in seq_along(cns)) {
  p <- p %>%
    add_trace(mode = 'lines', type = 'scatter3d', data = cns[[i]],
              x = ~x, y = ~y, z = ~z, opacity = .5,
              line = list(color = cns[[i]]$clr, width = 12),
              showlegend = FALSE)
}

p <- p %>%
  layout(scene = scn,
         font = list(family = 'Arial', size = 24, color = 'black'),
         paper_bgcolor = 'transparent',
         plot_bgcolor = 'transparent',
         legend = list(orientation="h"),
         showlegend = T)

orca(p, 'f2_c.pdf', height = 3000, width = 5000)

