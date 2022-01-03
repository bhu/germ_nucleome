library(tidyverse)
library(patchwork)
library(pals)
library(rtracklayer)
library(ggpattern)
library(scales)
library(gggenes)
library(scales)
sess <- browserSession("UCSC")
genome(sess) <- "mm10"

reg <- GRanges("chr4", IRanges(23090000, 23110000))

rmsk.tbl <- getTable(ucscTableQuery(sess, track = "rmsk", range = reg, table = "rmsk"))

bws <- list.files('../data/tracks/c19', pattern = 'bw$', full.names = T) %>% 
  setNames(., sub('ESC_(.*).bw', '\\1', basename(.))) %>% 
  lapply(function(x) {
    import.bw(x, selection = BigWigSelection(reg))
  })

bins <- tile(reg, 500)[[1]]
seqlevels(bins) <- seqlevels(bws[[1]])
xs <- mid(bins)

scs <- lapply(bws, function(x) {
  mcolAsRleList(x, 'score') %>%
    binnedAverage(bins, ., 'score') %>%
    score() %>%
    tibble(score = .) %>%
    mutate(idx = xs)
})

clrs <- kelly(8)[-2:-1]
ps <- seq_along(scs) %>%
  lapply(function(i) {
    p <- ggplot(scs[[i]], aes(x = idx, y = score)) + 
      geom_line(color = clrs[i]) +
      geom_area(fill = clrs[i], alpha = .5) +
      annotate('text', x = -Inf, y = Inf, label = names(scs)[i],
               hjust = -0.1, vjust = 1.5,) +
      scale_y_continuous(expand = expansion(c(0, .05)), breaks = pretty_breaks(2)) +
      scale_x_continuous(breaks = c(23090000, 23100000, 23110000)) +
      labs(x = seqnames(reg), y = 'CPM') +
      theme(legend.position = 'none',
            plot.background = element_blank(),
            panel.background = element_rect(fill = NA, color = 'black', size = 1),
            panel.grid = element_blank(),
            axis.text = element_text(color = 'black', size = 11),
            strip.background = element_rect(fill = NA),
            strip.text = element_text(color = 'black', size = 13, face = 'bold'),
            panel.grid.major = element_blank(),
            strip.background.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(t = 0, b = 2, r = 10),
            axis.title.y = element_blank(),
            axis.ticks.y = element_line(color = 'black'))
    if (i == 1) {
      p <- p +
        facet_grid(.~'Example c19 site in mESC') +
        theme(axis.title.y = element_text(color = 'black'))
    }
    p
  }) 


r <- rmsk.tbl %>% 
  filter(repName == 'L1Md_T') %>%
  dplyr::select(chr = genoName, start = genoStart, end = genoEnd, strand, name = repName, kind = repClass) %>%
  mutate(type = 'exon', utr = 'not', forward = strand == '+', 
         position = case_when(strand == '-' ~ end, T ~ start)) %>%
  ggplot(aes(xmin = start, xmax = end, y = kind, forward = forward)) +
  geom_gene_arrow(fill = 'black', color = NA) +
  scale_x_continuous(as.character(seqnames(reg)), limits = range(xs), breaks =  c(23090000, 23100000, 23110000),
                     labels = c('23.09', '23.1', '23.11mb')) +
  geom_text(aes(x = position, y = kind, label = name), nudge_x = 1e3, hjust = 0) +
  geom_feature(aes(x = position, y = kind, forward = forward)) +
  xlab(unique(reg$chr)) +
  coord_cartesian(clip = 'off') +
  scale_y_discrete(expand = expansion(0)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(t = 0, b = 2, r = 10),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank())



{c(ps, list(r)) %>%
  wrap_plots(ncol = 1, heights = c(1,1,1,1,1,1,.7)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf4_b.pdf', ., height = 3.7, width = 3)
 
