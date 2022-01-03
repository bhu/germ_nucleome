library(tidyverse)
library(pals)

samps <- c("d4c7PGCLC","GSC",  'GSCLC')

clrs <- setNames(c('#d62728', '#9467bd', '#499894'), samps)

load("../data/compscore/100kb.rda")

dat.score <- mats$mm10$Nagano[,samps] %>% 
  na.omit()  %>% 
  pivot_longer(everything(), names_to = 'x', values_to = 'y') %>%
  mutate(x = factor(x, rev(samps)))

dat.ratio <- dat.score %>%
  group_by(x) %>%
  summarise(y = 100 * sum(y > 0) / n(), .groups = 'drop')
m <- 0.2
b <- -10

sclr <- 'peru'
ggplot(dat.score, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_violin(aes(fill = x), alpha = .7, color = NA) +
  #geom_line(aes(y = m * y + b, group = 1), data = dat.ratio, color = sclr) +
  geom_point(aes(y = m * y + b, fill = x), data = dat.ratio, pch = 21, 
             color = sclr, size = 4, stroke = .7) +
  #coord_cartesian(ylim = c(-2, 2)) +
  scale_y_continuous('Compartment score',
                     sec.axis = sec_axis(~ (. - b) / m, '% of bins in compartment A')) +
  scale_fill_manual(values = clrs) +
  coord_flip(ylim = c(-2,2)) +
  scale_x_discrete(labels = function(x) sub('PGC', ' mPGC', x)) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.x.top = element_line(color = sclr),
        axis.ticks.x.top = element_line(color = sclr),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x.top = element_text(color = sclr),
        axis.title.x.top = element_text(color = sclr, vjust = 1.5),
        panel.grid = element_blank()) -> p
  ggsave('f6_e.pdf', ., height = 2, width = 5.03)

