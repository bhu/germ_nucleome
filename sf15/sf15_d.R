library(tidyverse)
library(pals)
library(matrixStats)
library(mgcv)
library(gratia)

load('../data/aggregate/TSS.K27me3.rda')

samps <- c("d4c7PGCLC","GSC",  'GSCLC') 
clrs <- setNames(c("#D62728", '#9467bd', '#499894'), samps)

a %>%
  separate(track, c('samp', 'mark'), '_') %>%
  filter(between(x, 21, 130)) %>%
  mutate(se = sd / sqrt(num)) %>%
  ggplot(aes(x = x, y = mu, color = samp, fill = samp)) +
  geom_vline(xintercept = c(50, 100), color = 'grey70', linetype = 'dashed') +
  annotate('rect', xmin = 50, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'tan',
           color = NA, alpha = .2) +
  #annotate('rect', xmin = 37, xmax = 47, ymin = -Inf, ymax = Inf, fill = 'royalblue',
  #         color = NA, alpha = .2) +
  geom_ribbon(aes(ymin = mu - se, ymax = mu + se), color = NA, alpha = .3) +
  geom_line() +
  scale_x_continuous(expand = expansion(0), 
                     breaks = c(21, 50, 100, 130),
                     labels = c('-15kb', 'TSS', 'TES', '+15kb')) +
  labs(x = 'Gene body', y = 'H3K27me3 (CPM)') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(.~'Male gamete generation: leading edge genes') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = c(.5,1),
        legend.justification = c(.5,1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.title = element_blank(),
        legend.background = element_rect(color = NA, fill = '#ffffff66'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        legend.key = element_blank(),
        strip.clip = 'off',
        legend.text = element_text(size = 11),
        plot.margin = margin(5,15,5,5)) -> p
ggsave('sf15_d.pdf', height = 4.93, width = 4.5)

