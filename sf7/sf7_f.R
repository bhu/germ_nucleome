library(tidyverse)
library(rtracklayer)
library(gggenes)
library(patchwork)
library(pals)

samps <- c('d4c7PGCLC',"GSC",  'GSCLC') 
rnm <- function(x) sub('PGC', '\nmPGC', x)
samps <- rnm(samps)
clrs <- setNames(c("#D62728", '#9467bd', '#499894'), samps)
sclrs <- c('+' = 'forestgreen', '-' =  'firebrick3')

g <- import.gff3('../data/resources/gencode.vM25.basic.protein_coding.annotation.gff3.gz')

reg <- GRanges('chr12', IRanges(102445000,102565000))

sg <- g[g$transcript_name %in% c('Golga5-202', 'Chga-201')] %>%
  as_tibble() %>%
  filter(type %in% c('exon', 'five_prime_UTR', 'three_prime_UTR', 'transcript')) %>%
  mutate(utr = ifelse(grepl('UTR', type), 'UTR', 'not'),
         forward = strand == '+',
         position = case_when(strand == '-' ~ end, T ~ start))
xbrks <- seq(102.45e6,102.55e6,5e4)
xlabs <- sprintf('%.2fmb', xbrks / 1e6)
ylims <- c(0, .25)
ybrks <- c(0, .2)

bws <- list.files('../data/tracks/K9me2', pattern = 'bw$', full.names = T) %>%
  setNames(., sub('_.*', '', basename(.))) %>%
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
  mutate(samp = rnm(samp))


p1 <- ggplot(scs, aes(x = idx, y = score, color = samp, fill = samp)) +
  #geom_line() +
  geom_area(alpha = .5) +
  #geom_hline(yintercept = ylims[2], color = 'black', size = 1) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(samp ~ 'H3K9me2 up-regulation') +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0), breaks = ybrks) +
  coord_cartesian(ylim = ylims) +
  ylab('Normalized coverage') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        strip.clip = 'off',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size  = 11))

g1 <- ggplot_gtable(ggplot_build(p1))
strip <- which(grepl('strip-r', g1$layout$name))
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- clrs[k]
  j <- which(grepl('title', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- clrs[k]
  k <- k+1
}


enhs <- list(dELS = '../data/resources/ccre/ccre/regions/dELS.bed',
     Enhancer = '../data/resources/ensembl/ensembl/regions/Enhancer.bed') %>%
  lapply(function(x) {
    import.bed(x) %>%
      subsetByOverlaps(reg) %>%
      as_tibble() %>%
      dplyr::select(1:3)
  }) %>%
  bind_rows(.id = 'ann') %>%
  mutate(y = as.numeric(factor(ann)))

ynms <- distinct(enhs, y, ann) 

p2 <- ggplot(enhs, aes(xmin = start, xmax = end, ymin = y - .2, ymax = y + .2)) +
  geom_rect(fill = 'black', color = 'black') +
  scale_y_continuous(breaks = ynms$y, labels = ynms$ann) +
  coord_cartesian(expand = F, clip = 'off', xlim = range(xs), ylim = c(.6,2.3)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        #axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color= 'black', size = 11),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text =element_text(color = 'black', size = 11),
        plot.margin = margin(t = 0))
  
p3 <- ggplot(sg, aes(xmin = start, xmax = end, y = 1, color = strand, fill = strand, forward = forward)) +
  # geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = 'gold2',
  #         color = NA, data = distinct(rcts, xmin, xmax), inherit.aes = F, alpha = .2) +
  geom_text(aes(x = position + diff(range(xs)) * .01, y = 1, label = gene_name), hjust = 1.2,
            data = ~subset(., type == 'transcript'), color = 'black') +
  geom_linerange(data = ~subset(., type == 'transcript')) +
  geom_gene_arrow(aes(alpha = utr), color = NA, data = ~subset(., type != 'transcript')) +
  scale_alpha_manual(values = c('UTR' = .5, 'not' = 1)) +
  geom_feature(aes(x = position, y = 1, forward = forward, color = strand),
               data = ~subset(., type == 'transcript')) +
  scale_color_manual(values = c('black', 'black')) +
  scale_fill_manual(values = c('black', 'black')) +
  scale_x_continuous(labels = xlabs, breaks = xbrks,
                     expand = expansion(0)) +
  xlab(unique(sg$seqnames)) +
   scale_y_continuous(breaks = 1, labels = 'Gene') +
  coord_cartesian( xlim = range(xs), clip = 'off') +
  #facet_grid('Gene'~.) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color= 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        axis.text =element_text(color = 'black', size = 11),
        plot.margin = margin(t = 0))

{ wrap_plots(wrap_ggplot_grob(g1), p2, p3, ncol = 1, heights = c(9,2,1)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf7_f.pdf', ., height = 4.9, width = 9)
