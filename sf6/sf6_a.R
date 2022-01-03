library(tidyverse)
library(ggh4x)
library(pals)
library(rtracklayer)

load('../data/tracks/inpnorm.50kb.rda')


bsz <- 5e5
rnm <- function(x) sub('PGC', ' mPGC', x)

pd <- d50 %>% 
  filter(samp == 'd4c7PGCLC' & chr == 'chr4') %>%
  select(chr, start, Laminb1, PC1) %>%
  mutate(start = round(start / bsz) * bsz) %>%
  group_by(chr, start) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(end = start + bsz,
         start = start + 1) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

reg <- reduce(pd, min.gapwidth = 1e9)

bws <- list.files('../data/tracks/K36me2', pattern = 'bw$', full.names = T) %>%
  setNames(., sub('_.*', '', basename(.))) %>%
  lapply(function(x) {
    import.bw(x, selection = BigWigSelection(reg))
  })

bins <- granges(pd)
seqlevels(bins) <- seqlevels(bws[[1]])
xs <- mid(bins)

scs <- lapply(bws, function(x) {
  mcolAsRleList(x, 'score') %>%
    binnedAverage(bins, ., 'score') %>%
    score() %>%
    tibble(score = .) %>%
    mutate(idx = xs)
}) %>% 
  bind_rows(.id = 'samp') %>%
  mutate(mark = 'K36me2') %>%
  select(idx, mark, y = score, samp) 

p <- as_tibble(pd) %>%
  mutate(idx = (start + end - 1) / 2) %>% 
  select(Laminb1, PC1, idx) %>% 
  pivot_longer(-idx, names_to = 'mark', values_to = 'y') %>%
  mutate(samp = 'd4c7PGCLC') %>%
  rbind(scs) %>%
  dplyr::rename(x = idx) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           sub('b1', ' B1', .) %>%
           factor(c('H3K36me2', 'PC1', 'Lamin B1')),
         samp = factor(rnm(samp), rnm(c('EpiLC', 'd2PGCLC', 'd4c7PGCLC'))),
         pos = ifelse(y > 0, y, 0),
         neg = ifelse(y < 0, y, 0),
         yint = 0) %>%
  na.omit() %>%
  ggplot(aes(x = x, y = y)) +
  geom_hline(aes(yintercept = yint), color = 'grey70', data = ~distinct(subset(., mark != 'H3K36me2'))) +
  #geom_line(aes(alpha = y > 0, group = 1), data = ~subset(., mark != 'PC1'), size = .25) +
  #geom_line(alpha = .5, data = ~subset(., mark == 'PC1'), size = .25) +
  scale_alpha_manual(values = c(.1,.5)) +
  scale_x_continuous('chr4', breaks = c(5e7, 1e8, 1.5e8),
                     limits = c(3e6, 156225000),
                     labels = c('50mb', '100mb', '150mb'),
                     expand = expansion(0)) +
  facet_nested(mark + samp~., scales = 'free_y', labeller = label_wrap_gen(width = 9)) +
  
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'firebrick', alpha = .5,
              data = ~ subset(., mark == 'PC1')) +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'forestgreen', alpha = .5,
              data = ~ subset(., mark == 'Lamin B1')) +
  geom_ribbon(aes(ymin = -1, ymax = pos), fill = 'gray10', alpha = .5,
              data = ~ subset(., mark == 'H3K36me2')) +
  
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'grey70', alpha = .1,
              data = ~ subset(., mark == 'Lamin B1')) +
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'steelblue', alpha = .5,
              data = ~ subset(., mark == 'PC1')) +

  geom_rug(aes(color = ifelse(y > 0, "Lamin B1", "bg")), sides = "b",
           data = ~ subset(., mark == 'Lamin B1')) +
  geom_rug(aes(color = ifelse(y > 0, "A", "B")), sides = "b",
           data = ~ subset(., mark == 'PC1')) +
  
  geom_vline(xintercept = Inf, color = 'black', size = 1) +
  scale_color_manual(values = c('Lamin B1' = 'forestgreen', H3K9me3 = 'darkorange3', 
                                H3K9me2 = 'darkorchid3', bg = 'grey70',
                                H3K36me2 = 'gray10', H3K27me3 = 'bisque',
                                A = 'firebrick', B = 'steelblue')) +
  ylab('Normalized coverage') +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, .15), breaks = c( 0, .1)),
    scale_y_continuous(limits = c(0, .15), breaks = c( 0, .1)),
    scale_y_continuous(limits = c(0, .12), breaks = c( 0, .1)),
    scale_y_continuous(limits = c(-7.5, 7.5), breaks = c(-5, 0, 5)),
    scale_y_continuous(limits = c(-.1, .2), breaks = c(0, .2))
  )) +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        strip.clip = 'off',
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) 

g <- ggplot_gtable(ggplot_build(p))
strip <- which(grepl('strip-r', g$layout$name))[c(1,3,4, 5,7)]
fills <- tableau20(7)[c(3,5,7,7,7)]
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}
ggsave('sf6_a.pdf', g, height = 4.5, width = 3.5)

