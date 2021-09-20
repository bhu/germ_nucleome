library(tidyverse)
library(pals)
library(GenomicRanges)
library(broom)
library(ggsignif)
library(patchwork)
library(boot)

load('../data/tracks/inpnorm.50kb.rda')
load('../data/calder.rda')

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)
mrks <- c('PC1','ATAC','CTCF','Rad21',
          'Ring1b','H2Aub','K27me3', 'K27ac', 'K4me1','K4me3',
          'K36me2','K36me3','K9me2','K9me3', 'Laminb1')

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}


sc <- lapply(o$`50000`, makeGRangesFromDataFrame, keep.extra.columns = T)

bins <- distinct(d50, chr, start, end) %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame()

lbs <- o$`50000` %>%
  lapply(function(x) {
    makeGRangesFromDataFrame(x) %>%
      findOverlaps(bins, .) %>%
      as("List") %>%
      extractList(x$lab, .) %>%
      sapply(`[`, 1)
  }) %>%
  bind_cols()


p <- d50 %>%
  split(., .$samp) %>%
  lapply(function(x) {
    x$sc <- x %>% 
      mutate(start = start + 1) %>%
      makeGRangesFromDataFrame() %>%
      findOverlaps(sc[[x$samp[1]]]) %>%
      as("List") %>%
      extractList(sc[[x$samp[1]]]$lab, .) %>%
      sapply(`[`, 1)
    x
  }) %>%
  bind_rows() %>%
  filter(samp %in% samps) %>%
  na.omit() %>%
  select(-chr, -start, -end) %>%
  mutate_if(is.numeric, function(x) (x - mean(x)) / sd(x)) %>%
  group_by(samp, sc) %>%
  summarise(across(everything(), mean), .groups = 'drop') %>%
  pivot_longer(-c(samp, sc), names_to = 'mark', values_to = 'mu') %>%
  mutate(mark = factor(mark, mrks)) %>%
  na.omit() %>%
  arrange(mark) %>%
  mutate(mark = sub('b1', ' B1', mark) %>%
           sub('^K', 'H3K', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x = sc, y = mark, fill = mu)) +
  geom_tile() +
  scale_fill_gradientn('Z-score', colors = turbo(25),
                       c(-1,0,1)) +
  facet_wrap(. ~ samp) +
  coord_cartesian(expand = F) +
  xlab('Subcompartment') +
  guides(fill = guide_colorbar(title.position = 'top',
                               barheight = .5,
                               barwidth = 2.5)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(1.4,1),
        legend.direction = 'horizontal',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
  
g <- ggplot_gtable(ggplot_build(p))
strip <- which(grepl('strip-t', g$layout$name))
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- clrs[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- clrs[k]
  k <- k+1
}
ggsave('sf6_c.pdf', g, height = 3.6, width = 3.9)

