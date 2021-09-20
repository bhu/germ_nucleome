library(tidyverse)
library(ggh4x)
library(rtracklayer)
library(pals)
library(patchwork)

samps <- c('ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1, 9, 2)], samps)

sess <- browserSession("UCSC")
genome(sess) <- "mm10"

reg <- GRanges('chr6', IRanges(85620e3, 85950e3))
rmsk.tbl <- getTable(ucscTableQuery(sess, track = "rmsk", range = reg, table = "rmsk"))

xlims <- c(85620e3, 85950e3)
xbrks <- c(85.7e6, 85.9e6)
xlabs <- c('85.7mb',  '85.9mb')

bws <- list.files('../data/tracks/K9me3', full.names = T) %>%
  setNames(., sub('_.*', '', basename(.))) %>%
  .[samps] %>%
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
}) %>% 
  bind_rows(.id = 'samp') %>%
  mutate(samp = factor(samp, samps),
         nm = sub('ESC', 'mESC', sub('PGC', ' mPGC', samp))) %>%
  na.omit()

p1 <- ggplot(scs, aes(x = idx, y = score)) +
  geom_line(aes(color = samp)) +
  geom_area(aes(fill = samp), alpha = .5) +
  geom_text(aes(label = nm, color = samp), data = distinct(scs, samp, nm), inherit.aes = F,
            hjust = 0, vjust = 1, fontface = 'bold', size = 4.4, #alpha = .5,
            x = start(reg) + width(reg) * .01, y = Inf) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(samp ~ 'H3K9me3 ChIP-seq') +
  scale_y_continuous(expand = expansion(c(0, .05)), breaks = c(0, .7)) +
  scale_x_continuous(as.character(seqnames(reg)), breaks = xbrks, labels = xlabs,
                   #  sec.axis = dup_axis(),
                     expand = expansion(0)) +
  coord_cartesian(xlim = xlims, clip = 'off') +
  labs(x = as.character(seqnames(reg)), y = 'Normalized coverage') +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x.bottom = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'black'),
        axis.line.y = element_line(color = 'black'))

rmsk.gr <- split(rmsk.tbl, rmsk.tbl$repClass) %>%
  lapply(function(x) {
    x %>%
      mutate(type = 'gene') %>%
      dplyr::select(chr = genoName, start = genoStart, end = genoEnd, strand, name = repName, type) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)
  })

pd <- rmsk.gr[c('LINE','SINE','LTR')] %>%
  lapply(as_tibble) %>%
  bind_rows(.id = 'family') %>% 
  mutate(idx = as.numeric(as.factor(family)))
p2 <- ggplot(pd, aes(xmin = start, xmax = end)) + 
  geom_rect(aes(fill = strand, ymin = idx - .45, ymax = idx + .45)) +
  scale_y_continuous(breaks = distinct(pd, idx, family)$idx,
                     labels = distinct(pd, idx, family)$family) +
  scale_fill_manual(values = c('#ff4500', '#00ced1')) +
  scale_x_continuous('Width',
                     breaks = xbrks, labels = xlabs) +
  coord_cartesian(xlim = xlims, expand = F) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        axis.title.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(size = 11),
        #legend.position = 'left',
        #legend.margin=margin(0, -30, 0, 0),
        axis.title.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks.x.top = element_blank())

load('../data/mcore.rda')

d <- d %>%
  mutate(logshift = log10(shift + 1),
         type = factor(type, samps)) %>%
  na.omit() %>%
  filter(mark == 'K9me3')
spls <- split(d, d$type) %>%
  lapply(function(x) {
    predict(smooth.spline(x$logshift, x$ctrb, spar = .5),
            seq(min(x$logshift), max(x$logshift),
                length.out = 100)) %>%
      bind_cols() %>%
      `names<-`(c('logshift', 'ctrb'))
  }) %>% 
  bind_rows(.id = 'type') %>%
  mutate(type = factor(type, samps),
         shift = 10^logshift - 1)

p3 <- ggplot(d, aes(x = logshift, y = ctrb, color = type, fill = type)) +
  geom_line(data = spls) +
  geom_area(data = spls, alpha = .5) +
  geom_point(aes(color = type), size = .2) +
  facet_grid(type ~ 'Domain width') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  coord_cartesian(xlim = c(1.2, 7.5), ylim = c(0, 1.2), expand = F, clip = F) + 
  scale_y_continuous(breaks = c(0,.5,1), labels = c('0', '0.5', '1')) +
  scale_x_continuous(breaks = c(2,4,6),
                     labels = c('100b', '10kb', '1mb')) +
  
  labs(x = 'Width', y = 'Contribution to global H3K9me3 enrichment') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        legend.position = 'none') 

wrap_plots(p1, p3,nrow=1, widths = c(1.5,1)) &
  theme(plot.background = element_blank()) &
  ggsave('f5_a_top.pdf', height = 3.72, width = 4.5)

ggsave('f5_a_bot.pdf', p2, height = 1, width = 3.2)

