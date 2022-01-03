library(data.table)
library(tidyverse)
library(rtracklayer)
library(ggpubr)
library(gghalves)
library(pals)

pks <- list.files('../data/peaks/K9me3', pattern = 'bed$', full.names = T) %>% 
  setNames(.,sub('.bed', '', basename(.))) %>%
  .[c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC')] %>%
  lapply(function(x) {
    fread(x, col.names = c('chr','start','end')) %>%
      mutate(start = start + 1, end = end + 1) %>% 
      makeGRangesFromDataFrame()
  }) 

rmsk <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/rmsk.txt.gz",
              col.names = c('bin', 'swScore', 'milliDiv', 'milliDel',
                            'milliIns', 'chr', 'start', 'end', 'left',
                            'strand',  'name', 'class', 'family', 
                            'repStart', 'repEnd', 'repLeft', 'id'))

rr <- makeGRangesFromDataFrame(rmsk)
i <- Reduce(intersect, pks)
u <- Reduce(union, pks)

d <- rmsk %>%
  mutate(kind = case_when(overlapsAny(rr, i) ~ 'Constitutive',
                          overlapsAny(rr, u) ~ 'Dynamic',
                          T ~ 'None'),
         age = milliDiv / 4.5e-9 / 1e3) %>%
  filter(!(class %in% c('Simple_repeat', 'Low_complexity')))

anns <- compare_means(age ~ kind, d) %>%
  mutate(y.position = 1e8 + (1:n() - 1) * 1e7)

clrs <- okabe()[1:3]

d %>%
  ggplot(aes(x = kind, y = age, color = kind, fill = kind)) +
  geom_boxplot(position = position_nudge(x = .2), width = .3) +
  stat_summary(geom = "crossbar", width = 0.2, fatten = 0, color = "white",
               position = position_nudge(x = .2),
               fun.data = function(x) c(y = median(x), ymin = median(x), ymax = median(x))) +
  geom_half_violin(color = NA , alpha = .75) +
  stat_pvalue_manual(anns, label = 'p.signif', inherit.aes = F) +
  scale_y_continuous(breaks = c(0, 5e7, 1e8),
                     limits = c(0, 1.3e8),
                     expand = expansion(0),
                     labels = c('0', '50', '100')) +
  labs(x = 'H3K9me3 enrichment across cell types', y = 'Estimated age (Million years)') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_wrap(~'Age of H3K9me3-marked transposons') +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed')) -> p
ggsave('sf5_c.pdf', p, height = 6.6, width = 5)
