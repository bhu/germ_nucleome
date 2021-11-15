library(tidyverse)
library(pals)
library(ggh4x)

odr <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(12)[seq(1, 10, 2)], odr)
cclrs <- kelly(22)[c(8,10,11,12,15,14,20)]
rnm <- function(x) case_when(x == 'ESC' ~ 'mESC',
                             x == 'd2PGCLC' ~ 'd2 mPGCLC',
                             x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
                             T ~ x)

load('../data/clust/pileup.rda')

pd <- lapply(cpup, function(y) {
    y %>% 
      mutate(y = 1:n()) %>% 
      pivot_longer(-y, names_to = 'x', values_to = 'z') %>% 
      mutate(x = as.integer(sub('^V', '', x))) 
}) %>%
  bind_rows(.id = 'f') %>%
  separate(f, c('samp', 'cl')) %>%
  filter(cl != '7') %>%
  mutate(cl = as.integer(cl),
         clu = c('Repressed', 'Polycomb', 'Active', 'Bivalent',
                 'CTCF', 'Readthrough', 'Enhancer')[cl]) %>%
  arrange(cl) %>%
  mutate(clu = fct_inorder(clu),
         samp = factor(rnm(samp), rnm(odr))) %>%
  na.omit() 


lab <- pd %>% 
  filter(between(x, 20, 22) & between(y, 20, 22)) %>%
  group_by(samp, clu) %>% 
  summarize(z = round(mean(z), 2), .groups = 'drop')

p <- ggplot(pd, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_gradientn('O/E', colors = rev(brewer.rdylbu(25)), trans = 'log2',
                       limits = c(1/1.6, 1.6), breaks = c(1/1.5, 1, 1.5),
                       oob = scales::squish,
                       labels = c('.67', '1.0', '1.5')) +
  geom_label(aes(x = -Inf, y = -Inf, hjust = 0, vjust = 1, label = z),
             data = lab, inherit.aes = F, alpha = .25, label.size = NA) +
  scale_y_reverse(breaks = 21) +
  coord_cartesian(expand = F) +
  facet_grid(samp ~ clu, labeller = label_wrap_gen(width=8)) +
  scale_x_continuous(breaks = 21) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        #axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.text = element_text(size = 11),
        plot.margin = margin(5,5,30,80),
        legend.position = c(0,1),
        legend.justification = c(1,.72)) +
  guides(fill = guide_colorbar(barheight = 3.2, barwidth = .5,
                               title.hjust = 1,
                               label.position = 'left',
                               title.vjust = 1.05))

g <- ggplot_gtable(ggplot_build(p))
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
strip <- which(grepl('strip-t', g$layout$name))
fills <- cclrs
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}
ggsave('f3_e.pdf', g, height = 6.2, width = 8.1, bg = 'transparent')

