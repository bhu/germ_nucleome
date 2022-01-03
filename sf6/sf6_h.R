library(tidyverse)
library(rtracklayer)

g <- read_tsv('../data/resources/mm10.chrom.sizes', col_names = F, col_types = cols()) %>%
  {Seqinfo(.$X1, .$X2)} %>%
  import.bed('../data/resources/Ns.100k.bed', seqinfo = .) %>%
  {setdiff(as(seqinfo(.), 'GRanges'), .)}

import.bed('../data/PMD/GSC.bed') %>%
  intersect(g) %>%
  as_tibble() %>%
  mutate(grp = 'PMD') %>%
  rbind(mutate(as_tibble(g), grp = 'GW')) %>%
  mutate(chr = case_when(seqnames == 'chrY' ~ 'Y',
                         seqnames == 'chrX' ~ 'X',
                         T ~ 'autosome')) %>%
  group_by(chr, grp) %>%
  summarise(w = sum(width), .groups = 'drop') %>%
  pivot_wider(names_from = 'grp', values_from = 'w') %>%
  mutate(PMD = PMD / GW,
         nonPMD = 1 - PMD) %>%
  select(-GW) %>%
  pivot_longer(-chr, names_to = 'grp', values_to = 'y') %>%
  ggplot(aes(x = chr, y = y * 100, fill = grp)) +
  geom_col() +
  scale_y_continuous(expand = expansion(0)) +
  scale_fill_manual(values = c('grey60', 'darkgoldenrod2')) +
  ylab('Percentage') +
  facet_grid(.~'PMD proportion') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.text = element_text(size = 11)) -> p
ggsave('sf6_h.pdf', height = 3.1, width = 4)
