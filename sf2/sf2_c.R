library(tidyverse)
library(readxl)
library(pals)
library(ggnewscale)
library(ggpubr)

wb <- read_tsv('../data/western/replicate.tsv', col_names = c('cell', 'mark', 'rep', 'value')) %>%
  group_by(mark, cell) %>% 
  summarise(cv = sd(value) / mean(value), .groups = 'drop') %>%
  mutate(cat = 'Western blot') %>%
  select(cv, cat)

pd <- read_excel('../data/histone_ratios.xlsx') %>%
  {.[,colSums(is.na(.)) != nrow(.)]} %>% 
  na.omit() %>% 
  {.[,c(1, which(.[1,] == 'Ratio'))]} %>%
  tail(-1) %>%
  rename(mark = 1) %>%
  filter(grepl('K9me[23]$|H3_.*K27me3$', mark)) %>% 
  pivot_longer(-mark, names_to = 'cell') %>%
  mutate(cell = sub('_.*', '', cell),
         value = as.numeric(value)) %>%
  group_by(mark, cell) %>% 
  summarise(cv = sd(value)/mean(value), .groups = 'drop') %>%
  mutate(cat = 'Mass spectrometry') %>%
  select(cv, cat) %>%
  rbind(wb)

light <- brewer.paired(4)[c(1,3)]
dark <- brewer.paired(4)[c(2,4)]

anns <- compare_means(cv ~ cat, pd) %>%
  mutate(y.position = .8)

ggplot(pd, aes(x = cat, y = cv, fill = cat, color = cat)) +
  geom_boxplot(outlier.colour = NA) +
  stat_summary(geom = "crossbar", width = 0.65, fatten = 0, color = "white", 
               fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  stat_pvalue_manual(data = anns, label = 'p.signif', inherit.aes = F, coord.flip = T, tip.length = 0) +
  scale_color_manual(values = light) +
  scale_fill_manual(values = light) +
  new_scale_color() +
  geom_jitter(aes(color = cat), width = .1, size = .2) +
  scale_color_manual(values = dark) +
  scale_y_continuous(breaks = c(0,.5,1)) +
  labs(x = "experiment",  y = "Coefficient of variation") +
  coord_flip() +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed')) +
  ggsave('sf2_c.pdf', height = 1.5, width = 4.3)
