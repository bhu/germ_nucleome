library(tidyverse)
library(pals)
library(ggh4x)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
samps <- rnm(samps)
clrs <- setNames(tableau20()[seq(1, 9, 2)], samps)

scls <- lapply(1:3, function(x) {
  scale_x_continuous(breaks = 1:5, labels = case_when(x == 1 ~samps, T ~ rep('', 5)))
})
read_csv('../data/western/summary.csv') %>% 
  mutate(cell = factor(rnm(cell), samps),
         x = as.numeric(cell)) %>% 
  filter(grepl('GLP|G9a|Setdb1', trg)) %>% 
  na.omit() %>%
  mutate(prot = sub('_.*', '', trg),
         bnd = sub('.*_', '', trg) %>%
           sub('L', 'G9a-S', .) %>%
           sub('H', 'G9a-L', .) %>%
           factor(c('G9a-L', 'G9a-S'))) %>%
  ggplot(aes(x = x, y = ratio)) + 
  geom_line(aes(group = 1), data = ~subset(., prot != 'G9a')) +
  geom_line(aes(linetype = bnd, group = bnd), data = ~subset(., prot == 'G9a')) +
  geom_point(aes(shape = bnd, color = cell), show.legend = F, data = ~subset(., prot == 'G9a')) +
  geom_point(aes(shape = bnd, y = ratio * 100), data = ~subset(., prot == 'G9a')) +
  geom_point(aes(color = cell), show.legend = F, data = ~subset(., prot != 'G9a')) +
  scale_color_manual(values = clrs) +
  ylab(expression('Target /'~alpha*'-Tubulin')) +
  facet_grid(.~prot, scales = 'free_x') +
  facetted_pos_scales(x = scls) +
  coord_cartesian(ylim = c(.2,2.2)) + 
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey90'),
        legend.title = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = '#ffffffaa', color = NA),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black')) -> p
ggsave('sf10_b.pdf',p, height = 2.55, width = 3.7)

