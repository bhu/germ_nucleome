library(tidyverse)
library(pals)
library(mgcv)
library(gratia)
library(ggnewscale)
library(patchwork)
library(ggdist)

load('../data/tracks/inpnorm.50kb.rda')

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
samps <- rnm(samps)
clrs <- setNames(tableau20(12)[seq(1, 10, 2)],  samps)
clrs2 <- c(setNames(tableau20(12)[seq(1, 10, 2)], paste('p', samps)),
           setNames(tableau20(12)[seq(2, 10, 2)], paste('t', samps)))

ld <- rep(c('pter (first 300kb)', 'qter (last 300kb)'), 10) %>%
  tibble(x = ., y = rnorm(20, sd = 1e-5))

scl <- .4
shf <- .1

d50 %>% 
  select(chr, start, samp, Laminb1) %>%
  mutate(samp = rnm(samp)) %>%
  pivot_wider(names_from = 'samp', values_from = 'Laminb1') %>%
  group_by(chr) %>%
  arrange(start) %>%
  mutate(i = 1:n(),
         j = n() - i + 1) %>%
  ungroup() %>%
  mutate(reg = case_when(i <= 6 ~ 'p',
                         j <= 6 ~ 't',
                         T ~ NA_character_)) %>%
  na.omit() %>%
  select(all_of(samps), reg) %>%
  pivot_longer(-reg, names_to = 'samp', values_to = 'y') %>%
  mutate(samp = factor(samp, samps),
         grp = paste(reg, samp)) %>%
  ggplot(aes(x = samp, y = y, color = grp, fill = grp)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  stat_slab(side = 'left', data = ~ subset(., reg == 'p'), color = NA, scale = scl,
            position = position_nudge(x = -shf), alpha = .3, show.legend = F) +
  stat_pointinterval(data = ~ subset(., reg == 'p'), show.legend = F,
                     position = position_nudge(x = -shf)) +
  stat_slab(side = 'right', data = ~ subset(., reg == 't'), color = NA, scale = scl,
            position = position_nudge(x = shf), alpha = .3, show.legend = F) +
  stat_pointinterval(data = ~ subset(., reg == 't'), show.legend = F,
                     position = position_nudge(x = shf)) +
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  new_scale_color() +
  new_scale_fill() +
  stat_slab(aes(x = x, y = y, fill = x), data = ld, color = NA,
            position = position_nudge(x = 10), alpha = .3, inherit.aes = F) +
  stat_pointinterval(aes(x = x, y = y, color = x), data = ld,
                     position = position_nudge(x = 10), inherit.aes = F) +
  scale_color_manual(values = c('grey10', 'grey70')) +
  scale_fill_manual(values = c('grey10', 'grey70')) +
  ylab('log2(ChIP / input)') +
  facet_grid(.~'Lamin B1 signal at opposing chromosome ends') +
  coord_cartesian(xlim = c(1,5)) +
  scale_x_discrete(limits = samps) +
  theme(legend.position = 'bottom',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 15, hjust = 1),
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) -> p
  ggsave('f4_j.pdf', p, height = 3.2, width = 5.3)

