library(tidyverse)
library(pals)
library(ggnewscale)
library(ggalluvial)

samps <- c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC") 

pal <- c("#ff0029", "#377eb8", "#66a61e", "#984ea3", "#00d2d5", "#ff7f00", "#af8d00", "#7f80cd", "#b3e900", "#c42e60", "#a65628", "#f781bf", "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#fccde5", "#bc80bd", "#ffed6f", "#c4eaff", "#cf8c00", "#1b9e77", "#d95f02", "#e7298a", "#e6ab02", "#a6761d", "#0097ff", "#00d067", "#f43600", "#4ba93b", "#5779bb", "#927acc", "#97ee3f", "#bf3947", "#9f5b00", "#f48758", "#8caed6", "#f2b94f", "#eff26e", "#e43872", "#d9b100", "#9d7a00", "#698cff", "#d9d9d9", "#00d27e", "#d06800", "#009f82", "#c49200", "#cbe8ff", "#fecddf", "#c27eb6", "#8cd2ce", "#c4b8d9", "#f883b0", "#a49100", "#f48800", "#27d0df", "#a04a9b")
load('../data/compscore/100kb.rda')
dat <- mats$mm10$Nagano[,samps] %>% 
  na.omit() %>%
  mutate_all(function(x) ifelse(x > 0, 'A', 'B')) %>%
  count(across(all_of(samps))) %>%
  ungroup() %>% 
  unite('grp', all_of(samps), remove = F) %>%
  mutate(n = 100 * n / sum(n)) %>%
  to_lodes_form(axes = samps) %>%
  mutate(x = as.character(x) %>% sub('PGCLC', ' mPGCLC', .) %>% sub('ESC','mESC', .) %>% fct_inorder())

ggplot(dat, aes(x = x, y = n, alluvium = alluvium, stratum = stratum)) +
  geom_alluvium(aes(fill = grp), alpha = .8) +
  scale_fill_manual(values = pal) +
  new_scale_fill() +
  geom_stratum(aes(fill = stratum), color ='black', alpha = .3) +
  scale_fill_manual(values = c(A = 'white', B = 'black')) +
  geom_text(stat = "stratum", aes(label = stratum)) +
  ylab('% of bins') +
  scale_y_continuous(expand = expansion(.01)) +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) -> p
ggsave('sf1_e.pdf', p, height = 3.7, width = 4.5)
