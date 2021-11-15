library(tidyverse)
library(pals)
library(rtracklayer)
library(patchwork)


samps <- c('EpiLC', 'GSC')
clrs <- setNames(tableau20()[c(3,9)], samps)


reg <- GRanges('chrY', IRanges(1, 9e7))
bws <- c('K9me3', 'Laminb1') %>%
  setNames(., .) %>%
  lapply(function(x) {
    samps %>%
      setNames(., .) %>%
      lapply(function(y) {
        import.bw(sprintf('../data/tracks/chrY/%s_%s.inpnorm.western.bw', y, x),
                  selection = BigWigSelection(reg))
      })
  })

bins <- tile(reg, width = 1e5)[[1]]
seqlevels(bins) <- seqlevels(bws[[1]][[1]])
xs <- mid(bins)

scs <- lapply(bws, function(m) {
  lapply(m, function(s) {
    mcolAsRleList(s, 'score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }) %>% 
    bind_rows(.id = 'samp')
}) %>%
  bind_rows(.id = 'mark')

load('../data/tracks/mCG.10kb.rda')
hits <- o %>% 
  mutate(start = start + 1) %>% 
  makeGRangesFromDataFrame() %>%
  findOverlaps(bins, .) %>%
  as("List")

dat <- samps %>%
  lapply(function(x) {
    tibble(mark = 'DNAme', samp = x, score = mean(extractList(o[[x]], hits), na.rm = T), idx = xs)
  }) %>%
  bind_rows() %>%
  rbind(scs) 



p <- o %>%
  mutate(chrom = case_when(!(chr %in% c('chrX', 'chrY')) ~ 'Autosome',
                           T ~ chr)) %>% 
  split(., .$chrom) %>%
  lapply(function(x) {
    p <- ggplot(x, aes(x = EpiLC, y = GSC)) +
      geom_density_2d_filled(aes(color = ..level..), bins = 20, show.legend = F) +
      geom_abline(slope = 1, intercept = 0, color = 'white', alpha = .5) +
      facet_grid(~chrom) + 
      scale_x_continuous('EpiLC mCG / CG per 10kb', breaks = c(0, .5, 1),
                         labels = c('0', '0.5', '1')) +
      scale_y_continuous('GSC mCG / CG per 10kb', breaks = c(0, .5, 1),
                         labels = c('0', '0.5', '1')) +
      coord_cartesian(xlim = 0:1, ylim = 0:1, expand = F) +
      theme(plot.background = element_blank(),
            panel.background = element_rect(fill = 'black', color = 'black', size = 1),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 13, face = 'bold', color = 'black'),
            axis.text = element_text(size = 11, color = 'black')) +
      scale_fill_viridis_d(option = 'A') +
      scale_color_viridis_d(option = 'A')
    
    if (x$chrom[1] == 'Autosome') {
      p + 
        theme(axis.title.x = element_blank())
    } else if (x$chrom[1] == 'chrX') {
      p +
        theme(axis.title.y = element_blank())
    } else {
      p + theme(axis.title = element_blank())
    }
  }) %>%
  wrap_plots(nrow = 1) 

ggsave('sf13_e.pdf', {p  & theme(plot.background = element_blank())}, height = 3.33, width = 8)  

