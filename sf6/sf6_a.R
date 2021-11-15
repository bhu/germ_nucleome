library(tidyverse)
library(pals)
library(ggnewscale)

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
  dplyr::count(samp, method) %>%
  mutate(samp = factor(samp, samps)) %>%
  na.omit() %>%
  arrange(samp) %>%
  ggplot(aes(x = samp, y = n)) +
  geom_boxplot(outlier.color = NA, show.legend = F) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  new_scale_color() +
  geom_point(aes(color = method), alpha = .5) +
  geom_line(aes(color = method, group = method), alpha = .5) +
  scale_color_manual(values = kelly(11)[-1]) +
  ylab('# of TAD boundaries') +
  facet_grid(.~'Insulating loci' ) +
  scale_x_discrete(labels = rnm) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(color = 'black', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed')) +
  #guides(color = guide_legend(nrow = 5)) +
  ggsave('sf6_a.pdf', height = 5, width = 7)

  
