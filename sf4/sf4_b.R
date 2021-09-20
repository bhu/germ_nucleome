library(tidyverse)
library(gprofiler2)
library(pals)
library(patchwork)
library(ggrepel)
library(ggtext)
library(ggrastr)
library(rrvgo)

load('../data/tss.1kb.MSnorm.rda')

pd <- bind_rows(tss.1kb.norm.mass[c('d2PGCLC', 'EpiLC')], .id = 'samp') %>%
  select(K9me3, samp, name) %>%
  pivot_wider(names_from = 'samp', values_from = 'K9me3') %>%
  mutate(kind = case_when(d2PGCLC >1 & (d2PGCLC-EpiLC) >1 ~ 'dn',
                          EpiLC >1 & (EpiLC-d2PGCLC) >1 ~ 'up',
                          T ~ 'ns'))

l <- c(-3, 5)
m <- mean(l)
anns <- pd %>% count(kind) %>% filter(kind != 'ns') %>% 
  mutate(x = ifelse(kind == 'dn', m, -Inf), y = ifelse(kind == 'dn', -Inf, m),
         hjust = .5, vjust = ifelse(kind == 'dn', 0, 1))

gs <- c("Mael", "Sycp3", "Ddx4", "Dazl", "Sycp2", "Sohlh2", "Sycp1", "Zp3", "Dmrtc2", "Terf1")

ggplot(pd, aes(x = d2PGCLC, y = EpiLC, color = kind)) +
  rasterize(geom_point(alpha = .2, size = .5, shape = 16, data = ~subset(., kind == 'ns')),dpi=600) +
  geom_point(data = ~subset(., kind != 'ns'), size = .5, alpha = .3, shape = 16) +
  geom_abline(slope = 1, intercept = 0, color = 'chartreuse3', alpha = .5) +
  geom_label(aes(x = x, y = y, label = n, hjust = hjust, vjust = vjust), data = anns) +
  geom_label_repel(aes(label = name, color = kind), data = ~subset(.,  name %in% gs),
                   direction = "y", nudge_x = 10) +
  scale_x_continuous(limits = l, expand = expansion(.01)) +
  scale_y_continuous(limits = l, expand = expansion(.01)) +
  scale_color_manual(values = c(dn = 'indianred', up = 'steelblue', ns = 'grey50')) +
  coord_cartesian( clip = F) +
  facet_grid(.~'Promoters depleted of H3K9me3') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed')) +
  ggsave('sf4_b.pdf', height = 6, width = 6)
