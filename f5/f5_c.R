library(tidyverse)
library(pals)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
samps <- rnm(samps)
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)
load('../data/peaks/K9me3/pileup.rda')
z1 <- d %>%
  filter(reg == 'inter' & between(x, 33, 66) & between(y, 33, 66)) %>%
  group_by(samp, reg) %>%
  summarise(z = mean(log2(z)))
z <- d %>%
  filter(reg == 'intra') %>%
  mutate(kind = case_when(between(x, 33, 66) & between(y, 33, 66) ~ 'intra',
                          between(x, 1, 33) & between(y, 33, 66) ~ 'inter',
                          between(x, 33, 66) & between(y, 1, 33) ~ 'inter',
                          T ~ NA_character_)) %>%
  filter(abs(x - y) > 3) %>%
  na.omit() %>%
  group_by(samp, kind) %>%
  summarise(z = mean(log2(z))) %>%
  pivot_wider(names_from = 'kind', values_from = 'z') %>%
  mutate(z = intra - inter,
         reg = 'intra') %>% 
  ungroup() %>%
  select(samp, reg, z) %>% 
  rbind(z1) %>%
  mutate(x = -Inf, y = -Inf,
         samp = factor(rnm(samp), samps),
         z = sprintf('%.3f', z),
         reg = c('intra' = 'Intra- domain', inter = 'Inter- domain')[reg]) %>%
  na.omit()

p <- d %>%
  mutate(samp = factor(rnm(samp), samps),
         v = log2(z),
         v = case_when(z == 0 ~ -.25,
                       T ~ v),
         v = case_when(reg == 'intra' & abs(x - y) < 3 ~ NA_real_,
                       T ~ v),
         reg = c('intra' = 'Intra- domain', inter = 'Inter- domain')[reg]) %>%
  na.omit() %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = v)) +
  geom_label(aes(label = z), data = z, hjust = 0, vjust = 1, label.size = NA, alpha = .25) +
  facet_grid(reg ~ samp, labeller = label_wrap_gen(7)) +
  scale_y_reverse() +
  scale_fill_gradientn('log2(O/E)', colors = rev(pals::brewer.rdylbu(10)), 
                       limits = c(-.25, .25), breaks = c(-.2, 0 , .2), oob = scales::squish) +
  guides(fill = guide_colorbar(title.vjust = 1, barheight = .6, barwidth = 4.25)) +
  coord_cartesian(expand = F) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.position = 'bottom',
        legend.margin = margin(-8,0,0,0),
        legend.background = element_blank(),
        legend.justification = 'right')

g <- ggplot_gtable(ggplot_build(p))
strip <- which(grepl('strip-t', g$layout$name))
fills <- clrs
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}

ggsave(filename = 'f5_c.pdf', g, height = 2.7, width = 5.1)


