library(tidyverse)
library(gprofiler2)
library(pals)
library(patchwork)
library(ggrepel)
library(gprofiler2)

load('../data/tss.1kb.MSnorm.rda')

pd <- bind_rows(tss.1kb.norm.mass[c('d2PGCLC', 'd4c7PGCLC')], .id = 'samp') %>%
  select(K9me3, samp, name) %>%
  pivot_wider(names_from = 'samp', values_from = 'K9me3') %>%
  mutate(kind = case_when(d2PGCLC >1 & (d2PGCLC-d4c7PGCLC) >1 ~ 'dn',
                          d4c7PGCLC >1 & (d4c7PGCLC-d2PGCLC) >1 ~ 'up',
                          T ~ 'ns'))

l <- c(-3, 5)
m <- mean(l)
anns <- pd %>% count(kind) %>% filter(kind != 'ns') %>% 
  mutate(x = ifelse(kind == 'dn', m, -Inf), y = ifelse(kind == 'dn', -Inf, m),
         hjust = .5, vjust = ifelse(kind == 'dn', 0, 1))

gs <- c("Mael", "Sycp3", "Ddx4", "Dazl", "Sycp2", "Sohlh2", "Sycp1", "Zp3", "Dmrtc2", "Terf1")

p1 <- ggplot(pd, aes(x = d2PGCLC, y = d4c7PGCLC, color = kind)) +
  geom_point(alpha = .2, size = .5, data = ~subset(., kind == 'ns')) +
  geom_point(data = ~subset(., kind != 'ns'), size = .5, alpha = .3) +
  geom_abline(slope = 1, intercept = 0, color = 'chartreuse3', alpha = .5) +
  geom_label(aes(x = x, y = y, label = n, hjust = hjust, vjust = vjust), data = anns) +
  geom_label_repel(aes(label = name, color = kind), data = ~subset(.,  name %in% gs),
                   direction = "y", nudge_x = 10) +
  scale_x_continuous(limits = l, expand = expansion(.01)) +
  scale_y_continuous(limits = l, expand = expansion(.01)) +
  scale_color_manual(values = c(dn = 'indianred', up = 'steelblue', ns = 'grey50')) +
  coord_cartesian( clip = F) +
  facet_grid(.~'Promoters depleted of H3K9me3') +
  labs(x = 'd2PGCLCs', y = 'd4c7PGCLCs') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed')) 


res <- gost(pd$name[pd$kind == 'dn'], 'mmusculus') %>%
  .$result %>%
  filter(source == 'GO:BP')

p2 <- res %>%
  arrange(precision) %>%
  mutate(term_name = fct_inorder(term_name)) %>%
  ggplot(aes(y = term_name, x = precision * 100, color = -log10(p_value))) +
  geom_linerange(aes(xmin = 0, xmax = precision * 100)) +
  geom_point(aes(size = intersection_size)) +
  labs(x = '% intersect') +
  scale_size('|Intersect|', breaks= c(10, 20, 30)) +
  scale_color_viridis_c(expression('-log'[10]~'p'['adj']),
                       breaks = c(2,4,6)) +
  scale_x_continuous(limits = c(0, 15), expand = expansion(0)) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  facet_grid(.~'Enriched GO:BP terms') +
  guides(color = guide_colorbar(barheight = 3.1, barwidth = .6)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        legend.position = c(0,1),
        legend.justification = c(6.5,1),
        plot.margin = margin(l = 70),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'))

wrap_plots(p1, p2, nrow = 1, widths = c(1.5,1)) &
  theme(plot.background = element_blank())  &
  ggsave(file = 'f5_j.pdf', height = 4.1, width = 14, device = cairo_pdf, bg = 'transparent')

