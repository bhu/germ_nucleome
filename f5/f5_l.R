library(tidyverse)
library(ggh4x)
library(rtracklayer)
library(pals)
library(patchwork)
library(ggridges)
library(ggh4x)

samps <- c('ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC') 
clrs <- setNames(tableau20(9)[seq(1, 9, 2)], samps)

g <- import.gff3('../data/resources/gencode.vM25.annotation.gff3.gz')
tss <- g[g$type == 'gene' & g$gene_name %in% c('Dazl', 'Ddx4')] %>%
  promoters(5e3, 5e3) %>%
  split(.,.$gene_name)

bws <- list.files('../data/tracks/K9me3', full.names = T, patter = 'bw$') %>%
  setNames(., sub('_.*', '', basename(.))) %>%
  .[samps] %>%
  lapply(function(x) {
    lapply(tss, function(reg) {
      import.bw(x, selection = BigWigSelection(reg))
    })
  })

scs <- lapply(tss, function(reg) {
  bins <- tile(reg, 500)[[1]]
  seqlevels(bins) <- seqlevels(bws[[1]][[1]])
  xs <- mid(bins)
  lapply(bws, function(x) {
    mcolAsRleList(x[[reg$gene_name]], 'score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }) %>% 
    bind_rows(.id = 'samp')
}) %>%
  bind_rows(.id = 'gene') %>%
  mutate(samp = factor(samp, samps),
         gene = factor(gene, c('Ddx4', 'Dazl'))) %>%
  na.omit()

anns <- distinct(scs, samp) %>%
  mutate(gene = factor('Ddx4', levels(scs$gene)),
         nm = sub('ESC', 'mESC', sub('PGC', ' mPGC', samp)),
         x = -Inf, y = Inf)

scs %>%
  mutate(x = sapply(tss,mid)[as.character(gene)] - idx) %>%
  ggplot(aes(x = x, y = score)) +
  geom_line(aes(color = samp)) +
  geom_area(aes(fill = samp), alpha = .5) +
  geom_text(aes(y = y, label = nm, color = samp), data = anns, inherit.aes = F,
            hjust = 0, vjust = 1.15, fontface = 'bold', size = 4.4,
            x = -4.9e3, y = Inf) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(samp ~ gene) +
  scale_y_continuous('Normalized H3K9me3 coverage', expand = expansion(c(0, .05)), breaks = c(0, .5)) +
  scale_x_continuous(breaks = c(-4e3,0,4e3),
                     labels = c('-4kb', 'TSS', '+4kb'),
                     expand = expansion(0)) +
  coord_cartesian(clip = 'off') +
  geom_hline(yintercept = 0, color = 'black') +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = .5),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x.bottom = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.line = element_blank()) +
  ggsave('f5_l.pdf', height = 3.5, width = 4.2)
  
