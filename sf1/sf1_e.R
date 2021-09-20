library(data.table)
library(tidyverse)
library(pals)
library(plotly)

samps <- c('ESC (2i)','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
           'ESC (serum)', 'Neural progenitors', 'Cortical neurons',
           'Day 6','Day 4','Day 2','B\u03B1')
rnm <- function(x) {
  sub('ESC', 'mESC', sub('PGCLC', ' mPGCLC', x))
}

samps <- rnm(samps)

clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(8)[-4]), samps)

load('../data/compscore/100kb.rda')
res <- c('Nagano', 'Stadhouders2018', 'Bonev2017') %>%
  lapply(function(s) {
    m <- mats$mm10[[s]]
    names(m) <- paste(s, names(m))
    m
  }) %>%
  bind_cols() %>%
  na.omit() %>%
  select(-contains('_')) %>%
  mutate(idx = 1:n()) %>%
  pivot_longer(-idx, names_to = 'x', values_to = 'y') %>%
  separate(x, c('study', 'x'), '\\ ') %>%
  mutate(x = case_when(x == 'Ba' ~ 'B\u03B1',
                       x == 'CM' ~ 'Cardiac mesoderm',
                       x == 'CPC' ~ 'Cardiac progenitors',
                       x == 'PCM' ~ 'Primitive cardiomyocytes',
                       x == 'VCM' ~ 'Ventricular cardiomyocytes',
                       x == 'ESC' & study == 'Nagano' ~ 'ESC (2i)',
                       x == 'ESC' & study == 'Bonev2017' ~ 'ESC (serum)',
                       x == 'NPC' ~ 'Neural progenitors',
                       x == 'CN' ~ 'Cortical neurons',
                       T ~ x) %>%
           sub('^D', 'Day ', .) %>%
           rnm() %>%
           factor(samps),
         study = c(Nagano = 'Germline', Bonev2017 = 'Neural',
                   Stadhouders2018 = 'Reprogramming')[study] %>%
           factor(c('Germline', 'Neural', 'Reprogramming'))) %>%
  na.omit() %>%
  arrange(idx, study, x) %>%
  select(-study) %>%
  pivot_wider(names_from = 'x', values_from = 'y') %>%
  select(-idx) %>%
  as.matrix() %>%
  t() %>%
  prcomp(scale. = T, center = T)

comps <- summary(res)$importance[2,] * 100

axs <- paste0('PC', 1:3)

pts <- res$x %>%
  data.frame() %>%
  rownames_to_column("samp") %>%
  mutate(samp = fct_inorder(samp))

lns <- lapply(1:(nrow(pts) - 1), function(i) {
  pts[c(i, i + 1), axs] %>%
    lapply(function(x) {seq(x[1], x[2], length.out = 100)}) %>%
    bind_cols() %>%
    `colnames<-`(c('x', 'y', 'z')) %>%
    mutate(idx = 1:n())
})

mz <- min(pts[[axs[3]]]) - 1
cns <- lapply(seq_along(pts), function(i) {
  pts[i, axs] %>%
    `colnames<-`(c('x', 'y', 'z')) %>%
    {rbind(., mutate(., z = mz))} %>%
    mutate(clr = clrs[pts$samp[i]])
})

scn <- lapply(axs, function(ax) {
  list(title = sprintf('%s (%d%%)', ax, round(comps[ax])),
       showticklabels = F, titlefont = list(size = 20))
}) %>% setNames(c('xaxis', 'yaxis', 'zaxis')) %>%
 c(list(#aspectratio= list(x= 1, y= 1, z= 0.8),
         camera = list(eye = list(x = 1.728, y = 1.143, z = 0.660),
                       center = list(x = .215, y = -.086, z = -.282),
                       up = list(x = 0, y = 0, z = 1))))

p <- plot_ly(source = 'src') %>%
  add_markers(type = 'scatter3d', x = pts[[axs[1]]],
              y = pts[[axs[2]]], z = pts[[axs[3]]],
              color = pts$samp, colors = clrs,
              marker = list(size = 14))

brks <- c(5, 8)

for (i in seq_along(lns)) {
  if (i %in% brks) next
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
              x = ~x, y = ~y, z = ~z,
              line = list(color = cns[[i]]$clr, dash = 'dash', width = 6),
              showlegend = FALSE)
}

p <- p %>%
  layout(scene = scn,
         paper_bgcolor='transparent',
         plot_bgcolor='transparent',
         showlegend=T,
         height = 600,
         font=list(family = "Arial", size = 20, color = 'black'))

orca(p, 'sf1_e.pdf', width = 1400,height=600, scale=1)
