library(data.table)
library(tidyverse)
library(ggdist)
library(ggpubr)
library(gghalves)
library(pals)
library(patchwork)

samps <- c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC", "GSCLC")
clrs <- dark <- setNames(tableau20(12)[seq(1, 12, 2)], samps)
light <- setNames(tableau20(12)[seq(2, 12, 2)], samps)

rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}

thm <- theme(legend.position = "none",
             plot.background = element_blank(),
             panel.background = element_blank(),
             strip.background = element_rect(fill = NA),
             strip.text = element_text(color = 'black', size = 13, face = 'bold'),
             axis.title.y = element_blank(),
             panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
             axis.line.x = element_line(color = 'black'),
             axis.ticks.y = element_blank(),
             axis.text = element_text(color = 'black', size = 11),
             panel.grid = element_blank())

load("../data/image/dense.dat.rda")
load("../data/image/DAPI_variance_all.rda")

dense.dat <- dense.dat %>% 
  filter(between(area, 500, 3000)) %>% 
  dplyr::filter(cell %in% samps) %>%
  mutate(cell = factor(cell, rev(samps)),
         dist = dist * 0.035,
         area = area * (0.035^2)) %>%
  filter(cell != "GSCLC")

stat.dist <- compare_means(dist ~ cell, data = dense.dat, method = 'wilcox.test') %>%
  mutate_at(c('group1','group2'), factor, samps) %>%
  dplyr::filter(((as.numeric(group1) - as.numeric(group2)) == 1) |
                  (group1 == 'd4c7PGCLC' & group2 == 'ESC')) %>%
  mutate(y.position = 6.5 + (1:n()) * .4)
 
p1 <- ggplot(dense.dat, aes(x = cell, y = dist)) +
  geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = cell), alpha = .7, scale = .9) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.dist, coord.flip = T, label = "p.signif",
                     tip.length = .01) +
  scale_fill_manual(values = light) +
  scale_color_manual(values = dark) +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  ylab(expression("Distance of DAPI dense regions to nuclear periphery (\u03BCm)")) +
  facet_grid('Radial positioning' ~ .) +
  coord_flip() +
  scale_x_discrete(labels = rnm) +
  thm
  
stat.area <- compare_means(area ~ cell, data = dense.dat, method = 'wilcox.test') %>%
  mutate_at(c('group1','group2'), factor, samps) %>%
  dplyr::filter(((as.numeric(group1) - as.numeric(group2)) == 1) |
                  (group1 == 'd4c7PGCLC' & group2 == 'ESC')) %>%
  mutate(y.position = 3.8 + (1:n()) * .18)

p2 <- dense.dat %>%
  ggplot(aes(x = cell, y = area)) +
  geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = cell), alpha = .7, scale = .9) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.area, coord.flip = T, label = "p.signif",
                     tip.length = .01) +
  scale_fill_manual(values = light) +
  scale_color_manual(values = dark) +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  ylab(expression("Area of DAPI dense regions (um"^2*")")) +
  facet_grid('Size' ~ .) +
  coord_flip() +
  scale_x_discrete(labels = rnm) +
  thm



var.dat <- dapi.varinace_all %>%
  filter(celltype %in% samps) %>%
  mutate(celltype = factor(celltype, rev(samps))) %>%
  filter(celltype != "GSCLC") %>%
  rename(cell = celltype)

stat.var <- compare_means(Var2 ~ cell, data = var.dat, method = 'wilcox.test') %>%
  mutate_at(c('group1','group2'), factor, samps) %>%
  dplyr::filter(((as.numeric(group1) - as.numeric(group2)) == 1) |
                  (group1 == 'd4c7PGCLC' & group2 == 'ESC')) %>%
  mutate(y.position = .2 + (1:n()) * .035)

p3 <- var.dat %>%
  ggplot(aes(x = cell, y = Var2)) +
  geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = cell), alpha = .7, scale = .9) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.var, coord.flip = T, label = "p.signif",
                     tip.length = .01) +
  scale_fill_manual(values = light) +
  scale_color_manual(values = dark) +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  ylab(expression("Variance of DAPI signal per slice")) +
  facet_grid('Uniformity' ~ .) +
  coord_flip() +
  scale_x_discrete(labels = rnm) +
  thm

{wrap_plots(p2, p1, p3, ncol = 1) & 
  theme(plot.background = element_blank())} %>%
  ggsave('f1_c.pdf', ., height = 8, width = 5.17, device = cairo_pdf, bg = 'transparent')

