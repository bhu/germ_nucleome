library(data.table)
library(tidyverse)
library(tidygraph)
library(igraph)
library(GenomicRanges)
library(pals)

load('../data/cliques/sig.ints.rda')
load('../data/tracks/inpnorm.50kb.rda')

cs <- o$germ$d4c7PGCLC %>%
  mutate(i = sprintf('%s:%d-%d', V1, V2, V3),
         j = sprintf('%s:%d-%d', V4, V5, V6)) %>%
  select(i, j) %>%
  as.matrix() %>%
  graph_from_edgelist() %>%
  max_cliques() %>%
  lapply(function(y) {
    tibble(node = names(y)) %>%
      add_count()
  }) %>%
  bind_rows(.id = 'idx') %>%
  group_by(node) %>%
  slice_max(n, n = 1, with_ties = F) %>%
  ungroup() %>%
  separate(node, c('chr', 'start', 'end'), '[:-]') %>%
  mutate(start = as.numeric(start) + 1) 

d <- d50 %>%
  filter(samp == 'd4c7PGCLC') %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

cs$K9me3 <- makeGRangesFromDataFrame(cs) %>%
  findOverlaps(d) %>%
  as("List") %>%
  extractList(d$K9me3, .) %>%
  mean()

cs %>%
  filter(n > 2) %>%
  mutate(x = factor(n)) %>%
  ggplot(aes(x = x, y = K9me3)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_violin(aes(fill = n), color = NA, alpha = .5) +
  geom_boxplot(aes(color = n, fill = n), outlier.color = NA, width = .25) +
  stat_summary(geom = "crossbar", width = 0.2, fatten = 0, color = "white", 
               fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  scale_color_viridis_c(option = 'E', end = .85) +
  scale_fill_viridis_c(option = 'E', end = .85) +
  scale_y_continuous(breaks = -1:1) +
  labs(x = 'Size of max TAD clique',
       y = 'H3K9me3 in TAD') +
  facet_grid(.~'d4c7PGCLC TADs') +
  theme(legend.position = 'none',
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#ffffff66", color = NA),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13,face = 'bold'),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank()) +
  ggsave(file = 'f6_h.pdf', width = 2.4, height = 1.72)
