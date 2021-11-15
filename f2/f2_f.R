library(tidyverse)

samps <- c('ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}

load('../data/tads/10.methods.rda')

d %>%
  count(samp, method) %>%
  mutate(samp = factor(samp, samps)) %>%
  na.omit() %>%
  group_by(method) %>%
  mutate(idx = as.character(length(samps) - rank(n) + 1)) %>%
  ungroup() %>%
  arrange(idx) %>%
  mutate(rank = c('1' = 'Most boundaries',
                  '2' = '2nd most',
                  '3' = '3rd most')[idx] %>%
           fct_inorder(),
         samp = factor(samp, samps)) %>%
  na.omit() %>%
  group_by(samp) %>%
  arrange(idx) %>%
  mutate(y = 1:n()) %>%
  ggplot(aes(x = samp, y = y)) +
  #geom_label(aes(color = rank, label = idx), label.r = unit(.5, "lines"), 
  #           fontface = 'bold', label.size = .6) +
  geom_point(aes(color = rank), size = 3) +
  scale_color_manual(values = c('#FAD649', '#C0C0C0', '#C49A6D')) +
  scale_y_continuous('# of TAD callers', breaks = c(2,4,6,8,10)) +
  facet_grid(.~'Insulation ranking') +
  scale_x_discrete(labels = rnm) +
  coord_cartesian(clip = 'off') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color=  'grey70', linetype = 'dashed'),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        legend.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        axis.text.x = element_text(angle = 30, hjust = 1),
        #legend.position = c(1,1),
        #legend.justification = c(1,1),
        #legend.position = 'bottom',
        strip.text = element_text(color = 'black', size = 13, face =  'bold')) -> p
  ggsave('f2_f.pdf', p, height = 3.1, width = 5)

