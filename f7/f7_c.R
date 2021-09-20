library(data.table)
library(tidyverse)
library(ggdist)
library(ggpubr)
library(gghalves)
library(pals)
library(patchwork)
library(gghalves)

samps <- c('d4c7PGCLC',"GSC",  'GSCLC') 
clrs <- dark<-  setNames(c("#D62728", '#9467bd', '#499894'), samps)
light <- setNames(c("#FF9896", "#C5B0D5", "#98c5c2"), samps)

thm <- theme(legend.position = "none",
             plot.background = element_blank(),
             panel.background = element_rect(fill = NA, color = 'black', size = 1),
             strip.background = element_rect(fill = NA),
             strip.text = element_text(color = 'black', size = 13, face = 'bold'),
             axis.title.x = element_blank(),
             panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
             #axis.line.x = element_line(color = 'black'),
             axis.ticks.y = element_blank(),
             strip.clip = 'off',
             axis.text.x = element_text(angle = 30, hjust = 1),
             axis.text = element_text(color = 'black', size = 11),
             panel.grid = element_blank())

load("../data/image/dense.dat.rda")

dense.dat <- dense.dat %>% 
  filter(between(area, 500, 3000)) %>% 
  dplyr::filter(cell %in% samps) %>%
  mutate(cell = factor(cell, samps),
         dist = dist * 0.035,
         area = area * (0.035^2)) 

stat.dist <- compare_means(dist ~ cell, data = dense.dat, method = 'wilcox.test') %>%
  mutate_at(c('group1','group2'), factor, samps) %>%
  mutate(y.position = 6.5 + (1:n()) * .4)

p1 <- ggplot(dense.dat, aes(x = cell, y = dist)) +
  #geom_violin(aes(fill = cell), color = NA, alpha = .5) +
  #geom_boxplot(aes(color = cell, fill = cell), outlier.color = NA, width = .15, notch = T) +
  #stat_summary(geom = "crossbar", width = 0.07, fatten = 0, color = "white", 
  #             fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = cell), alpha = .7, scale = .9) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.dist, coord.flip = F, label = "p.signif",
                     tip.length = .01) +
  scale_fill_manual(values = light) +
  scale_color_manual(values = dark) +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  ylab(expression("Distance of DAPI dense regions\n   from nuclear periphery (\u03BCm)")) +
  facet_grid(.~'Radial positioning') +
  scale_x_discrete(labels = function(x) sub('PGC', ' mPGC', x)) +
  thm +
  theme(plot.margin = margin(l = 15))

stat.area <- compare_means(area ~ cell, data = dense.dat, method = 'wilcox.test') %>%
  mutate_at(c('group1','group2'), factor, samps) %>%
  mutate(y.position = 3.6 + (1:n()) * .2)

p2 <- dense.dat %>%
  ggplot(aes(x = cell, y = area)) +
  #geom_violin(aes(fill = cell), color = NA, alpha = .5) +
  #geom_boxplot(aes(color = cell, fill = cell), outlier.color = NA, width = .15, notch = T) +
  #stat_summary(geom = "crossbar", width = 0.07, fatten = 0, color = "white", 
  #             fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = cell), alpha = .7, scale = .9) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.area, coord.flip = F, label = "p.signif",
                     tip.length = .01) +
  scale_fill_manual(values = light) +
  scale_color_manual(values = dark) +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  ylab(expression("Area of DAPI dense regions (um"^2*")")) +
  facet_grid(. ~ 'Size') +
  scale_x_discrete(labels = function(x) sub('PGC', ' mPGC', x)) +
  thm

load("../data/image/DAPI_variance_all.rda")

var.dat <- dapi.varinace_all %>%
  filter(celltype %in% samps) %>%
  mutate(celltype = factor(celltype, samps)) %>%
  rename(cell = celltype) 

stat.var <- compare_means(Var2 ~ cell, data = var.dat, method = 'wilcox.test') %>%
  mutate_at(c('group1','group2'), factor, samps) %>%
  mutate(y.position = .2 + (1:n()) * .035)

p3 <- var.dat %>%
  ggplot(aes(x = cell, y = Var2)) +
  geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.1)) +
  stat_slab(aes(fill = cell), alpha = .7, scale = .9) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.var, coord.flip = F, label = "p.signif",
                     tip.length = .01) +
  scale_fill_manual(values = light) +
  scale_color_manual(values = dark) +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  ylab(expression("Variance of DAPI signal per slice")) +
  facet_grid(. ~ 'Uniformity') +
  scale_x_discrete(labels = function(x) sub('PGC', ' mPGC', x)) +
  thm

wrap_plots(p2, p1, p3, nrow = 1) & 
    theme(plot.background = element_blank()) &
  ggsave('f7_c.pdf', width = 6.2, height = 3.58, device = cairo_pdf, bg = 'transparent')

