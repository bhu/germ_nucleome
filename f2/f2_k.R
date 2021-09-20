library(tidyverse)
library(pals)

load('../data/compscore/25kb.rep.rda')

samps <- c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC") 
clrs <- setNames(tableau20(10)[seq(1, 10, 2)], samps)
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}

bsz <- 25000

dm <- o %>%
  mutate(idx = as.integer(start / bsz)) %>%
  separate(samp, c('samp', 'rep'), '_') %>%
  arrange(samp, chr, idx) 

dl <- distinct(dm, chr, idx) %>%
  split(., .$chr) %>%
  lapply(function(x) {
    group_by(x, gw = cumsum(c(1, diff(idx) != 1))) %>%
      summarise(n = n(), i = idx[1], j = idx[n]) %>%
      ungroup()
  }) %>%
  bind_rows(.id = 'chr') %>%
  group_by(chr) %>%
  top_n(1, n) %>%
  ungroup() %>%
  mutate(l = n * bsz / 1e6)

ac <- split(dm, dm$samp) %>%
  lapply(function(x) {
    split(x, x$chr) %>%
      lapply(function(y) {
        y %>%
          dplyr::filter(y$idx > dl[dl$chr == y$chr[1],]$i &
                          y$idx < dl[dl$chr == y$chr[1],]$j) %>%
          {split(.$score, .$rep)} %>%
          lapply(function(z) {
            acf(z, na.action = na.pass, plot = F,
                lag.max = dl[dl$chr == y$chr[1],]$n) %>%
              {tibble(lag = .$lag[,,1], c = .$acf[,,1])} %>%
              dplyr::filter(lag >= 0) %>%
              mutate(lag = lag * bsz,
                     c = abs(c))
          }) %>%
          bind_rows(.id = 'rep')
      }) %>%
      bind_rows(.id = 'chr') %>%
      dplyr::filter(!(chr %in% c('chrY'))) %>%
      group_by(lag, rep) %>%
      summarise(s = sd(c, na.rm = T) / sqrt(20),
                c = mean(c, na.rm = T),
                lo = c - s, hi = c + s) 
  }) %>%
  bind_rows(.id = 'samp') %>%
  na.omit() %>%
  dplyr::filter(lag != 0) %>%
  mutate(samp = factor(samp, samps)) %>%
  arrange(samp) %>%
  na.omit()

ggplot(ac, aes(x = lag, y = c, color = samp)) + 
  geom_line(aes(color = samp, group = interaction(samp, rep))) +
  coord_cartesian(xlim = c(bsz, 5e6)) +
  scale_color_manual(values = clrs, labels = rnm) +
  labs(x = 'Shift', y = 'Compartment score\nPearson auto-correlation') +
  facet_grid(.~'Compartment broadness') +
  scale_x_continuous(breaks = c(0,1e6,2e6,3e6,4e6,5e6),
                     labels = c('0', '1mb', '2mb', '3mb', '4mb', '5mb')) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        legend.title = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.ticks.y = element_blank()) +
  ggsave('f2_k.pdf', height = 2.8, width = 3.26)

