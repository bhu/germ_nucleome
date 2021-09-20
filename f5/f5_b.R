library(data.table)
library(tidyverse)
library(plotly)
library(ggalluvial)

state <- function(a, b) {
  case_when(a == b ~ a,
            a == 'S' ~ b,
            b == 'S' ~ a,
            a == 'E' & b == 'C' ~ 'M',
            a == 'C' & b == 'E' ~ 'M')
}

clrs <- c(Missing = 'grey70',
          Contraction = 'indianred',
          Expansion = 'steelblue',
          Steady = 'olivedrab3')

fread('../data/peaks/K9me3/ChromTime.bed.gz', sep = '\t', select = 1:5) %>%
  separate(V5, c('u1', 'u2', 'u3', 'u4',
                 'd1', 'd2', 'd3', 'd4'), '[/-]') %>%
  mutate(s1 = state(u1, d1),
         s2 = state(u2, u2),
         s3 = state(u3, u3),
         s4 = state(u4, u4),
         s = sprintf('%s-%s-%s-%s', s1, s2, s3, s4)) %>%
  dplyr::filter(!grepl('M', s) & (V3 - V2) >= 1e4) %>%
  group_by(s) %>%
  summarise(num = n()) %>%
  separate(s, c('s1', 's2', 's3', 's4'), '-', remove = T) %>%
  mutate_at(paste0('s', 1:4), function(x) {
    factor(c('x' = 'Missing', C = 'Contraction', S = 'Steady', E = 'Expansion')[x],
           c('Missing', 'Contraction', 'Steady', 'Expansion'))
  }) %>%
  to_lodes_form(key = 'segment', value = 'value', id = 'id', axes = 1:4) %>%
  ggplot(aes(x = segment, stratum = value, alluvium = id, y = num)) +
  geom_flow(aes(fill = value), width = 1/12, aes.flow = "backward",
            alpha = 0.5) +
  geom_stratum(aes(fill = value, color = value), width = 1/12) +
  scale_x_discrete(labels = c('mESC\u2192EpiLC', 'EpiLC\u2192d2 mPGCLC', 'd2\u2192d4c7 mPGCLC', 'd4c7 mPGCLC\u2192GSC'),
                   expand = c(.05, .05, .05, .05)) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  labs(x = 'Transition', y = '# of H3K9me3 domains') +
  scale_y_continuous(expand = expansion(0), breaks = c(0,1e4,2e4)) +
  facet_wrap(~'H3K9me3 domain kinetics') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 13, face =  'bold'),
        axis.text = element_text(color = 'black', size = 11, family = 'Arial'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80',
                                          linetype = 'dashed'),
        axis.ticks.y = element_blank()) +
  ggsave('f5_b.pdf', height = 2.8, width = 3.5, device = cairo_pdf, bg = 'transparent')

