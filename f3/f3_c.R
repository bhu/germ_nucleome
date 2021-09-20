library(data.table)
library(tidyverse)
library(ggalluvial)
library(ggrepel)
library(pals)
odr <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
odrr <- c('repressed', 'CTCF', 'enhancer',  'bivalent', 'promoterCTCF', 'promoter', 'Missing')
cclrs <- setNames(c(kelly(22)[2:7], 'grey80'), odrr)
load('../data/clust/open.rda')
d %>%
  filter(cl != -1 & type %in% odr) %>%
  mutate(V1 = cl,
         cl = case_when(V1 %in% c(1,10,11) ~ 'promoterCTCF',
                        V1 %in% c(7,8,9,12) ~ 'bivalent',
                        V1 %in% c(6,13) ~ 'promoter',
                        V1 %in% c(4,5,14,17) ~ 'enhancer',
                        V1 %in% c(15,16,18) ~ 'repressed',
                        V1 %in% c(0,2,3) ~ 'CTCF',
                        T ~ NA_character_)) %>%
  select(samp = type, cl, crd) %>%
  pivot_wider(names_from = 'samp', values_from = 'cl', values_fill = 'Missing') %>%
  count(across(all_of(odr))) %>% 
  unite('grp', all_of(odr), remove = F) %>%
  to_lodes_form(axes = odr) %>%
  mutate(stratum = factor(stratum, odrr)) %>%
  ggplot(aes(x = x, stratum = stratum, alluvium = grp, 
             y = n, fill = stratum, label = stratum)) + 
  geom_flow() +
  geom_stratum(alpha = .7, size = .25) +
  #geom_text(stat = 'stratum') +
  scale_x_discrete(expand = expansion(c(.4,0.05)), labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  geom_text_repel(aes(label = ifelse(as.numeric(x) == 1, as.character(stratum), NA), color = stratum),
                  stat = "stratum", size = 4, direction = "y", nudge_x = -1, nudge_y = -10) +
  scale_fill_manual(values = cclrs) +
  scale_color_manual(values = cclrs) +
  facet_grid(.~ 'Open site chromatin state dynamics') +
  scale_y_continuous('Number of open sites', expand = expansion(c(0, 0.01)),
                     breaks = c(0, (1:3)*1e5), labels = c('0', '100k', '200k', '300k')) +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank()) +
  ggsave('f3_c.pdf', width = 5.1, height = 3.2)
