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
sg <- g[g$transcript_name %in% c('Dmrt1-201', 'Dmrt3-201')] %>%
  as_tibble() %>%
  filter(type %in% c('exon', 'five_prime_UTR', 'three_prime_UTR', 'transcript')) %>%
  mutate(utr = ifelse(grepl('UTR', type), 'UTR', 'not'),
         forward = strand == '+',
         position = case_when(strand == '-' ~ end, T ~ start),
         grp = case_when(grepl('Dmrt', gene_name) ~ 'Dmrt cluster',
                         T ~ gene_name)) 

reg <- GRanges('chr19', IRanges(25475e3, 25665e3))
xbrks <- c(25.5e6,25.6e6)
xlabs <- c('25.5', '25.6mb')
ylims <- c(0, .6)
ybrks <- c(0, .2, .4)

bws <- list.files('../data/tracks/Dmrt', pattern = 'bw$', full.names = T) %>%
  setNames(., sub('_.*', '', basename(.))) %>%
  lapply(function(x) {
    import.bw(x, selection = BigWigSelection(reg))
  })

bins <- tile(reg, 300)[[1]]
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

rcts <- sg %>%
  filter(type == 'transcript') %>% 
  mutate(xmin = case_when(strand == '+' ~ start - 2000,
                          T ~ end -2000),
         xmax = case_when(strand == '+' ~ start + 2000,
                          T ~ end + 2000)) %>%
  {lapply(unique(scs$samp), function(y) mutate(., samp = y))} %>%
  bind_rows()


p1 <- ggplot(scs, aes(x = idx, y = score, color = samp, fill = samp)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = 'gold2',
            color = NA, data = rcts, inherit.aes = F, alpha = .4) +
  geom_area(alpha = .5) +
  geom_hline(yintercept = ylims[2], color = 'black', size = 1) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  facet_grid(samp ~ 'Residual\nH3K27me3') +
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
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.text.x = element_text(color = 'black'),
        strip.background.x = element_rect(fill =NA),
        axis.text = element_text(color = 'black'))

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


p2 <- ggplot(sg, aes(xmin = start, xmax = end, y = 1, color = strand, fill = strand, forward = forward)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = 'gold2',
            color = NA, data = distinct(rcts, xmin, xmax), inherit.aes = F, alpha = .4) +
  geom_linerange(data = ~subset(., type == 'transcript')) +
  geom_gene_arrow(aes(alpha = utr), color = NA, data = ~subset(., type != 'transcript')) +
  scale_alpha_manual(values = c('UTR' = .5, 'not' = 1)) +
  geom_feature(aes(x = position, y = 1, forward = forward, color = strand),
               data = ~subset(., type == 'transcript')) +
  scale_color_manual(values = c('black', 'black')) +
  scale_fill_manual(values = c('black', 'black')) +
  scale_x_continuous(limits = range(xs), labels = xlabs, breaks = xbrks,
                     expand = expansion(0)) +
  xlab(unique(sg$seqnames)) +
  coord_cartesian(clip = 'off') +
  #facet_grid('Gene'~.) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color= 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        axis.text =element_text(color = 'black', size = 11),
        plot.margin = margin(t = 0)) +
  geom_text(aes(x = start - diff(range(xs)) * .01, y = 1, label = gene_name), hjust = 1,
            data = ~subset(., type == 'transcript' & gene_name == 'Dmrt1'), color = 'black') +
  geom_text(aes(x = end + diff(range(xs)) * .01, y = 1, label = gene_name), hjust = 0,
            data = ~subset(., type == 'transcript' & gene_name == 'Dmrt3'), color = 'black')


wrap_plots(wrap_ggplot_grob(g1), p2, ncol = 1, heights = c(5,1)) &
  theme(plot.background = element_blank()) &
  ggsave('f7_i.pdf', height = 5.25, width = 3.5)

