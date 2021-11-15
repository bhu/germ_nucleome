library(tidyverse)
library(rtracklayer)
library(pals)
library(patchwork)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', '\nmPGC', x))
samps <- rnm(samps)
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

reg <- GRanges('chr6', IRanges(124774268, 124852544))

bws <- list.files('../data/tracks/sf3', pattern = 'bw$', full.names = T) %>%
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
  mutate(samp = factor(rnm(samp), samps)) %>%
  na.omit()

kept <- lapply(bws, function(x) {
  tibble(start = c(124846815,124840705,124808980),
         end = c(124847990, 124841749, 124810233)) 
}) %>% bind_rows(.id = 'samp') %>%
  mutate(samp = factor(rnm(samp), samps)) %>% 
  na.omit()

lost <- lapply(bws, function(x) {
  tibble(start = c(124850660,124827974,124779051),
         end = c(124851964, 124834714, 124807027)) 
}) %>% bind_rows(.id = 'samp') %>%
  mutate(samp = factor(rnm(samp), samps)) %>% 
  na.omit()


p <- ggplot(scs, aes(x = idx, y = score)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'green',
            alpha = .3, data = kept, inherit.aes = F) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'red',
            alpha = .15, data = lost, inherit.aes = F) +
  geom_line(aes(color = samp)) +
  geom_area(aes(fill = samp), alpha = .5) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  coord_cartesian(ylim = c(0,6)) +
  facet_grid(samp ~ 'Effect on weaker CTCF sites') +
  scale_y_continuous(expand = expansion(c(0, .05)), breaks = c(0, 2.5, 5)) +
  scale_x_continuous(breaks = c(124780000, 124840000),
                     labels = c('124.78mb', '124.84mb')) +
  labs(x = as.character(seqnames(reg)), y = 'CPM') +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major = element_blank(),
        panel.spacing = unit(.9, "lines"),
        strip.clip = 'off',
        axis.ticks.y = element_line(color = 'black'),
        axis.line.y = element_line(color = 'black')) 

g <- ggplot_gtable(ggplot_build(p))
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


ggsave('sf8_c.pdf', g, height = 6.2, width = 3.6)

