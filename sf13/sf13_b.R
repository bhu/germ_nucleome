library(tidyverse)
library(pals)
library(rtracklayer)
library(ggh4x)
library(patchwork)

samps <- c('EpiLC','GSC')
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

bins <- tile(reg, 500)[[1]]
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

p1 <- scs %>%
  dplyr::rename(y = score) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           sub('b1', ' B1', .),
         pos = ifelse(y > 0, y, 0),
         neg = ifelse(y < 0, y, 0)) %>%
  ggplot(aes(x = idx, y = y)) +
  #geom_line(aes(alpha = y > 0, group = 1)) +
  scale_alpha_manual(values = c(.1,.5)) +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'forestgreen', alpha = .5,
              data = ~ subset(., mark == 'Lamin B1')) +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'darkorange3', alpha = .5,
              data = ~ subset(., mark == 'H3K9me3')) +
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'grey70', alpha = .1,
              data = ~ subset(., mark == 'Lamin B1')) +
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'grey70', alpha = .1,
              data = ~ subset(., mark == 'H3K9me3')) +
  geom_vline(xintercept = Inf, color = 'black', size = 1) +
  scale_color_manual(values = c('Lamin B1' = 'forestgreen', H3K9me3 = 'darkorange3', 
                                 bg = 'grey70')) +
  facet_nested(mark + samp ~ ., scales = 'free_y', nest_line = T) +
  scale_x_continuous('chrY', expand = expansion(0),
                     breaks = c(2.5e7, 5e7, 7.5e7),
                     labels = c('25mb', '50mb', '75mb')) +
  coord_cartesian(xlim = c(5e6, 9e7)) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 2), breaks = c(0,1), expand = expansion(0)),
    scale_y_continuous(limits = c(0, 2), breaks = c(0,1), expand = expansion(0)),
    scale_y_continuous(limits = c(0, .8), breaks = c( 0, .5), expand = expansion(0)),
    scale_y_continuous(limits = c(0, .8), breaks = c( 0, .5), expand = expansion(0))
  )) +
  ylab('log2 enrichment over input') +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.clip = 'off',
        axis.title.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        plot.margin = margin(b = 0))
g1 <- ggplot_gtable(ggplot_build(p1))
strip <- which(grepl('strip-r', g1$layout$name))[c(1,3,4,6)]
fills <- clrs[c(samps,samps)]
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}

load('../data/tracks/mCG.10kb.rda')

w <- median(diff(sort(unique(scs$idx))))
p2 <- o %>% 
  filter(chr == 'chrY') %>%
  mutate(x = (start + end) / 2) %>%
  select(-chr, -start, -end) %>%
  pivot_longer(-x, names_to = 'samp', values_to = 'y') %>%
  filter(samp %in% samps) %>%
  split(.,.$samp) %>%
  lapply(function(x) {
    x %>%
      na.omit() %>%
      {smooth.spline(.$x, .$y, nknots = 1000)} %>% 
      predict(x = sort(unique(scs$idx))) %>% 
      {tibble(x = .$x, y = .$y)}
  }) %>%
  bind_rows(.id = 'samp') %>%
  mutate(mark = 'DNAme') %>%
  ggplot(aes(x = x, y = y)) +
  geom_rect(aes(xmin = x - w, xmax = x + w, ymin = -Inf, ymax = Inf, fill = y)) +
  facet_nested(mark + samp~., nest_line = T) +
  scale_fill_gradientn('mCG/CG', limits = 0:1, colors = rev(brewer.rdbu(10)),
                       breaks = c(0,.5,1)) +
  scale_x_continuous('chrY', expand = expansion(0),
                     breaks = c(4e7, 8e7),
                     labels = c('40mb', '80mb'),
                     sec.axis = dup_axis()) +
  coord_cartesian(xlim = c(5e6, 9e7)) +
  guides(fill = guide_colorbar(title.vjust = 1, barheight = 3.2, barwidth = .5, title.position = 'left')) +
  theme(plot.background = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(1,-.1),
        plot.margin = margin(l=20, t=1),
        strip.clip = 'off',
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        axis.text.x.top = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.title.x.top = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(angle = 90),
        #axis.line = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'))
p2

g2 <- ggplot_gtable(ggplot_build(p2))
strip <- which(grepl('strip-r', g2$layout$name))[c(1,3)]
fills <- clrs[samps]
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}

wrap_plots(wrap_ggplot_grob(g1), wrap_ggplot_grob(g2), nrow = 2, heights = c(2.2,1)) &
  theme(plot.background = element_blank()) &
  ggsave('sf13_b.pdf', height = 3.8, width = 3.2)

ggsave('sf13_b_top.pdf', g1, height = 2, width = 3.6)
ggsave('sf13_b_bot.pdf', g2, height = 1.75, width = 3.98)

