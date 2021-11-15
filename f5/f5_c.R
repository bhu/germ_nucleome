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
  summarise(z = mean(z))
z <- d %>%
  filter(reg == 'intra') %>%
  mutate(kind = case_when(between(x, 33, 66) & between(y, 33, 66) ~ 'intra',
                          between(x, 1, 33) & between(y, 33, 66) ~ 'inter',
                          between(x, 33, 66) & between(y, 1, 33) ~ 'inter',
                          T ~ NA_character_)) %>%
  filter(abs(x - y) > 3) %>%
  na.omit() %>%
  group_by(samp, kind) %>%
  summarise(z = mean(z)) %>%
  pivot_wider(names_from = 'kind', values_from = 'z') %>%
  mutate(z = intra/inter,
         reg = 'intra') %>% 
  ungroup() %>%
  select(samp, reg, z) %>% 
  rbind(z1) %>%
  mutate(samp = factor(rnm(samp), samps),
         reg = c('intra' = 'Intra- domain', inter = 'Inter- domain')[reg]) %>%
  na.omit()

p <- d %>%
  mutate(samp = factor(rnm(samp), samps),
         z = case_when(reg == 'intra' & abs(x - y) < 3 ~ NA_real_, T ~ z),
         reg = c('intra' = 'Intra- domain', inter = 'Inter- domain')[reg]) %>%
  na.omit() %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = z)) +
  facet_grid(reg ~ samp, labeller = label_wrap_gen(7)) +
  scale_y_reverse() +
  scale_fill_gradientn('O/E', colors = rev(brewer.rdylbu(25)), 
                       limits = c(1/1.2,1.2), breaks = c(.84, 1, 1.2),
                       labels = c('0.83', '1', '1.2'),
                       oob = scales::squish, trans = 'log2') +
  guides(fill = guide_colorbar(title.vjust = 1, barheight = 4, barwidth = .5,
                               label.position = 'left', title.hjust = 1)) +
  coord_cartesian(expand = F) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.position = 'left',
        legend.background = element_blank(),
        legend.justification = 'top')

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

ggsave(filename = 'f5_c.pdf', g, height = 2.3, width = 5.5, bg = 'transparent')

z %>% 
  ggplot(aes(x = samp, y = z)) +
  geom_line(aes(linetype = reg, group = reg)) +
  geom_point(aes(shape = reg, y = z * 100), size = 3) +
  geom_point(aes(shape = reg, color = samp), show.legend = F, size = 3) +
  scale_color_manual(values = clrs) +
  coord_cartesian(ylim = c(1.03,1.17)) +
  scale_x_discrete(expand = expansion(add = .45)) +
  ylab('O/E interaction') +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey90'),
        legend.title = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,.8),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(size = 11),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black')) -> p

ggsave('f5_c_bot.pdf', height = 2.2, width = 4.63)

