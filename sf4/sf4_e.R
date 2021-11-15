library(tidyverse)
library(rtracklayer)
library(pals)
library(patchwork)
library(readxl)

load('../data/tracks/K9me2.rda')

smps <- c('EpiLC', 'd4c7PGCLC')
nms <- c('EpiLC', 'd4c7 mPGCLC')
clrs <- setNames(tableau20(8)[c(3,7)], nms)

fac <- read_excel('../data/histone_ratios.xlsx') %>%
  {.[,colSums(is.na(.)) != nrow(.)]} %>% 
  na.omit() %>% 
  {.[,c(1, which(.[1,] == 'Ratio'))]} %>%
  tail(-1) %>%
  dplyr::rename(mark = 1) %>%
  filter(grepl('K9me2$', mark)) %>%
  pivot_longer(-mark, names_to = 'samp') %>%
  mutate(type = sub('_.*', '', samp),
         value = as.numeric(value)) %>%
  filter(type %in% smps)  %>%
  group_by(type) %>% 
  summarise(v = mean(value)) %>% 
  mutate(v = v / min(v)) %>% 
  deframe()

k <- !overlapsAny(makeGRangesFromDataFrame(bins), import.bed('../data/resources/blacklist.bed'))
pclr <- 'steelblue3'
p2 <- ggplot(dat[k,], aes_string(x = smps[2], y = smps[1])) +
  geom_abline(slope = 1) +
  geom_point(alpha = .1, size = .1, color = pclr) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = nms[2], y = nms[1]) +
  coord_cartesian(xlim = c(10, 150), ylim = c(10, 150)) +
  annotation_logticks() +
  facet_grid(.~'Original') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.clip = 'off',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size = 11))

p4 <- dat[k,] %>% 
  dplyr::rename(y := !!smps[1],
                x := !!smps[2]) %>%
  mutate(y = y * fac[smps[1]], 
         x = x * fac[smps[2]]) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1) +
  geom_point(alpha = .1, size = .1, color = pclr) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = nms[2], y = nms[1]) +
  coord_cartesian(xlim = c(10, 200), ylim = c(10, 200)) +
  annotation_logticks() +
  facet_grid(. ~ 'Rescaled') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        axis.title.y = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size = 11))




p1 <- dat[k & bins$chr == 'chr6',] %>%
  mutate(idx = 1e5 * (1:n() - .5 )) %>%
  pivot_longer(-idx, names_to = 'samp', values_to = 'v') %>%
  mutate(samp = setNames(nms, smps)[samp] %>%
           factor(nms)) %>%
  ggplot(aes(x = idx, y = v, color = samp, fill = samp)) +
  geom_line() +
  geom_area(alpha = .5) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_x_continuous(breaks = c(0, 5e7, 1e8, 1.5e8),
                     labels = c('0b', '50mb', '100mb', '150mb'),
                     expand = expansion(0)) +
  facet_grid(samp ~ 'Original') +
  scale_y_continuous(expand = expansion(c(0,.05)),
                     breaks = c(25,50,75)) +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = 'chr6', y = 'CPM') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        legend.position = 'none',
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size = 11))

p3 <- dat[k & bins$chr == 'chr6',] %>%
  dplyr::rename(y := !!smps[1],
                x := !!smps[2]) %>%
  mutate(y = y * fac[smps[1]], 
         x = x * fac[smps[2]],
         idx = 1e5 * (1:n() - .5 )) %>%
  pivot_longer(-idx, names_to = 'samp', values_to = 'v') %>%
  mutate(samp = setNames(nms, c('y', 'x'))[samp] %>%
           factor(nms)) %>%
  ggplot(aes(x = idx, y = v, color = samp, fill = samp)) +
  geom_line() +
  geom_area(alpha = .5) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_x_continuous(breaks = c(0, 5e7, 1e8, 1.5e8),
                     labels = c('0b', '50mb', '100mb', '150mb'),
                     expand = expansion(0)) +
  facet_grid(samp ~ 'Rescaled') +
  scale_y_continuous(expand = expansion(c(0,.05)),
                     breaks = c(50,100,150)) +
  coord_cartesian(ylim = c(0,200)) +
  labs(x = 'chr6', y = 'Normalized CPM') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        legend.position = 'none',
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size = 11))

g <- ggplot_gtable(ggplot_build(p3))
strip <- which(grepl('strip-r', g$layout$name))
fills <- clrs
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}

{wrap_plots(p1, wrap_ggplot_grob(g), nrow = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf4_e_left.pdf', ., height = 3.35, width = 5)

{wrap_plots(p2, p4, nrow = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf4_e_right.pdf', ., height = 3.35, width = 5)

