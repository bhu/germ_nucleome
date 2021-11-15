library(tidyverse)
library(pals)
library(ggh4x)
library(eulerr)
library(rtracklayer)

load('../data/tracks/inpnorm.50kb.rda')
load('../data/tracks/mCG.10kb.rda')

chrom <- 'chr3'
i <- 32e6
j <- 47e6

cd <- d50 %>%
  filter(samp %in% c('d4c7PGCLC', 'GSC')) %>%
  filter(chr == chrom & start > i & end < j) %>%
  select(samp, K9me3, K9me2, Laminb1, start) %>%
  pivot_longer(-c(samp, start), names_to = 'mark', values_to = 'y')

p <- o %>% 
  filter(chr == chrom & start > i & end < j) %>% 
  select(start, GSC) %>%
  mutate(start = floor(start / 5e4) * 5e4) %>%
  group_by(start) %>%
  summarise(GSC = mean(GSC)) %>%
  pivot_longer(-start, names_to = 'samp', values_to = 'y') %>%
  mutate(mark = 'mCG') %>%
  select(samp, start, mark, y) %>%
  rbind(cd) %>%
  mutate(start = start + 25e3) %>%
  filter(!((samp == 'd4c7PGCLC' & mark == 'Laminb1') | (samp == 'GSC' & mark == 'K9me2'))) %>%
  na.omit() %>%
  split(., .$samp) %>%
  lapply(function(x) {
    split(x, x$mark) %>%
      lapply(function(y) {
        smooth.spline(y$start, y$y, nknots = 295) %>% 
          predict(x = seq(i, j, length.out = 1e3)) %>% 
          {tibble(x = .$x, y = .$y)}
      }) %>% bind_rows(.id = 'mark')
  }) %>% bind_rows(.id = 'samp') %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           #mark %>%
           sub('b1', ' B1', .),
         pos = ifelse(y > 0, y, 0),
         neg = ifelse(y < 0, y, 0),
         yint = 0) %>%
  filter(samp == 'GSC') %>%
  ggplot(aes(x = x, y = y)) +
  geom_hline(aes(yintercept = yint), color = 'grey70', data = ~distinct(subset(., mark != 'mCG'))) +
  geom_line(aes(alpha = y > 0, group = 1), size = .5) +
  scale_alpha_manual(values = c(.1,.5)) +
  scale_x_continuous(chrom, breaks = c(36e6,40e6,44e6),
                     labels = c('36mb', '40mb', '44mb'),
                     expand = expansion(0)) +
  facet_nested(samp + mark ~ ., scales = 'free_y') +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'forestgreen', alpha = .5,
              data = ~ subset(., mark == 'Lamin B1')) +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'darkorange3', alpha = .5,
              data = ~ subset(., mark == 'H3K9me3')) +
  
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'gold2', alpha = .5,
              data = ~ subset(., mark == 'mCG')) +
  
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'grey70', alpha = .1,
              data = ~ subset(., mark == 'Lamin B1')) +
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'grey70', alpha = .1,
              data = ~ subset(., mark == 'H3K9me3')) +
  
  geom_rug(aes(color = ifelse(y > 0, "Lamin B1", "bg")), sides = "b",
           data = ~ subset(., mark == 'Lamin B1')) +
  geom_rug(aes(color = ifelse(y > 0, "H3K9me3", "bg")), sides = "b",
           data = ~ subset(., mark == 'H3K9me3')) +
  geom_rug(aes(color = ifelse(y > .75, "mCG", "bg")), sides = "b",
           data = ~ subset(., mark == 'mCG')) +
  
  scale_color_manual(values = c('Lamin B1' = 'forestgreen', H3K9me3 = 'darkorange3', 
                                H3K9me2 = 'darkorchid3', bg = 'grey70',
                                mCG = 'gold2')) +
  
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(-1, .5), breaks = c(-1,0,.5), expand = expansion(.1)),
    scale_y_continuous(limits = c(-.1, .2), breaks = c(-.1, 0, .2), expand = expansion(.1)),
    scale_y_continuous(limits = c(.5, 1), breaks = c(.5, 1), expand = expansion(.1))
  )) +
  ylab('log2 enrichment over input') +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) 

g <- ggplot_gtable(ggplot_build(p))
strip <- which(grepl('strip-r', g$layout$name))[2]
fills <- tableau20()[9]
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}  
ggsave('sf12_d_bot.pdf', g, height = 3.0, width = 3.5)  

pmds <- c(GSC = 'GSC', SSC = 'SSC') %>%
  lapply(function(x) {
    import.bed(sprintf('../data/PMD/%s.bed', x)) %>%
      {.[seqnames(.) != 'chrY']}
  })

u <- distinct(o, chr, start, end) %>% makeGRangesFromDataFrame()

pdf('sf12_d_top.pdf', bg = 'transparent', width = 2.8,height = 1.3)
lapply(pmds, function(x) overlapsAny(u, x)) %>%
  bind_cols() %>%
  euler() %>%
  {plot(., quantities = sprintf('%.1fMb', .$original.values / 100))}
dev.off()


  
