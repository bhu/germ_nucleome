library(data.table)
library(tidyverse)
library(pals)
library(ggbeeswarm2)

load("../data/tracks/coverage.perchr.rda")

dat %>% 
  separate(samp, c('samp', NA, 'mark'), '_') %>%
  group_by(samp, mark) %>%
  mutate(ratio = cov / sum(cov)) %>%
  ungroup() %>%
  {merge(filter(., mark != 'Input'), filter(., mark == 'Input'),
         by = c('chr', 'samp'), suffixes = c('', '.input'))} %>%
  mutate(r = ratio / ratio.input,
         kind = case_when(chr == 'chrX' ~ 'X',
                          chr == 'chrY' ~ 'Y',
                          T ~ 'autosome')) %>%
  filter(mark %in% c("K9me3", "Laminb1") & samp != 'GSCLC') %>%
  mutate(samp = factor(samp, c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC")),
         mark = sub('^K', 'H3K', mark) %>%
           sub('b1', ' B1', .)) %>%
  ggplot(aes(x = samp, y = r, color = kind)) +
  geom_beeswarm(aes(shape = kind), spacing = .7, size = 2) +
  facet_wrap(~mark, scales = "free") +
  scale_color_manual(values = c(autosome = '#808080', X = '#01BFC4', Y = '#F8766D')) +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = c(.7,1),
        legend.justification = c(.5,1),
        #legend.background = element_rect(fill = '#ffffff66', color = NA),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_rect(fill = NA),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank()) +
  labs(y = "Coverage ratio (IP/Input)") +
  ggsave('f6_i.pdf', height = 3.8, width = 3.4)

