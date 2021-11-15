library(tidyverse)
library(rtracklayer)
library(pals)
library(patchwork)
library(DiffBind)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) {
  sub('ESC','mESC', sub('PGCLC', ' mPGCLC', x))
}
samps <- rnm(samps)
clrs <- setNames(tableau20(9)[seq(1, 9, 2)], samps)

thm <- theme(legend.position = 'none',
             plot.background = element_blank(),
             panel.background = element_blank(),
             axis.ticks.y = element_blank(),
             panel.grid = element_blank(),
             axis.text = element_text(color = 'black', size = 11),
             panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
             strip.background = element_rect(fill = NA),
             strip.text = element_text(color = 'black', size = 13, face = 'bold'),
             strip.clip = 'off',
             axis.line.x = element_line(color = 'black'))

reg <- GRanges(seqnames = c('chr9'), IRanges(start = 20883000, end = 20957000))
bws <- list.files('../data/tracks/f2', full.names = T) %>%
  setNames(., basename(.)) %>%
  lapply(function(x) {
    list.files(x, full.names = T) %>%
      setNames(., sub('.bw', '', basename(.))) %>%
      lapply(function(y) {
        import.bw(y, selection = BigWigSelection(reg))
      })
  })

bins <- tile(reg, 500)[[1]]
seqlevels(bins) <- seqlevels(bws[[1]][[1]])
xs <- mid(bins)

scs <- lapply(bws, function(x) {
  lapply(x, function(y) {
    mcolAsRleList(y, 'score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }) %>% 
    bind_rows(.id = 'samp') %>%
    mutate(samp = factor(rnm(samp), samps)) %>%
    na.omit()
})

pks <- list.files('../data/peaks/ATAC', pattern = 'bed', full.names = T) %>%
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(function(x) {
    import.bed(x) %>%
      subsetByOverlaps(reg) %>%
      as_tibble() %>%
      select(1:3)
  }) %>%
  bind_rows(.id = 'samp') %>%
  mutate(samp = factor(rnm(samp), samps)) %>%
  na.omit()

p1 <- ggplot(scs$ATAC, aes(x = idx, y = score)) +
  geom_line(aes(color = samp)) +
  geom_area(aes(fill = samp), alpha = .5) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            data = pks, inherit.aes = F, fill = 'black', alpha = .1) +
  facet_grid(samp ~ 'ATAC-seq') +
  scale_y_continuous(expand = expansion(c(0, .05)), breaks = c(0, 3, 6)) +
  scale_x_continuous(breaks = c(20.9e6, 20.94e6)) +
  labs(x = as.character(seqnames(reg)), y = 'Normalized coverage') +
  thm +
  theme(panel.grid.major = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.ticks.y = element_line(color = 'black'),
        axis.line.y = element_line(color = 'black'))

p3 <- ggplot(scs$K4me1, aes(x = idx, y = score)) +
  geom_line(aes(color = samp)) +
  geom_area(aes(fill = samp), alpha = .5) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(samp ~ 'H3K4me1 ChIP-seq') +
  scale_y_continuous(expand = expansion(c(0, .05)), breaks = c(0, 1, 2)) +
  scale_x_continuous(breaks = c(20.9e6, 20.94e6)) +
  labs(x = as.character(seqnames(reg)), y = 'Normalized coverage') +
  thm +
  theme(panel.grid.major = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.ticks.y = element_line(color = 'black'),
        axis.line.y = element_line(color = 'black'))


load('../data/peaks/ATAC/fpkm.rda')

p2 <- atac.my.pks %>%
  ungroup() %>% 
  mutate(cell = factor(rnm(as.character(cell)), samps)) %>%
  na.omit() %>%
  ggplot(aes(x = mean, fill = cell, color = cell)) +
  geom_density(alpha = .5) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(cell ~ 'Openness') +
  thm +
  scale_y_continuous(expand = expansion(c(0,.05))) +
  labs(x = 'log2(FPKM + 1)', y = 'Density') +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())


load('../data/mcore.rda')

d <- d %>%
  mutate(logshift = log10(shift + 1),
         type = factor(rnm(type), samps)) %>%
  na.omit() %>%
  filter(mark == 'K4me1')
spls <- split(d, d$type) %>%
  lapply(function(x) {
    predict(smooth.spline(x$logshift, x$ctrb),
            seq(min(x$logshift), max(x$logshift),
                length.out = 100)) %>%
      bind_cols() %>%
      `names<-`(c('logshift', 'ctrb'))
  }) %>% 
  bind_rows(.id = 'type') %>%
  mutate(type = factor(type, samps),
         shift = 10^logshift - 1)

p4 <- ggplot(d, aes(x = logshift, y = ctrb, color = type, fill = type)) +
  geom_line(data = spls) +
  geom_area(data = spls, alpha = .5) +
  geom_point(aes(color = type), size = .2) +
  facet_grid(type ~ 'Domain width') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  coord_cartesian(xlim = c(1.2, 7.5), ylim = c(0, 1.2), expand = F) + 
  scale_y_continuous(breaks = c(0,.5,1)) +
  scale_x_continuous(breaks = c(2,4,6),
                     labels = c('100b', '10kb', '1mb')) +
  labs(x = 'Width', y = 'Contribution to global H3K4me1 enrichment') +
  thm

g <- ggplot_gtable(ggplot_build(p4))
strip <- which(grepl('strip-r', g$layout$name))
fills <- clrs
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}


{wrap_plots(p1, p2, p3, g, nrow = 1, widths = c(2.5,1,2.5,1.8)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('f2_c.pdf',., height = 5,width=10)
