library(tidyverse)
library(pals)
library(LOLA)
library(GenomicRanges)
library(gratia)
library(mgcv)
library(pals)
load('../data/tracks/inpnorm.50kb.rda')

regionDB <- loadRegionDB('../data/resources/reps')
regs <- regionDB$regionGRL %>% 
  setNames(sub('.bed', '', regionDB$regionAnno$filename)) %>%
  .[c('L1', 'ERVK', 'ERV1')] %>%
  unlist()

bsz <- 1e6
r <- d50 %>% 
  distinct(chr, start, end) %>%
  filter((start %% bsz) == 0) %>%
  mutate(end = start + bsz,
         start = start + 1) %>%
  makeGRangesFromDataFrame() %>%
  {mutate(as_tibble(.), rep = countOverlaps(., regs))} %>% 
  mutate(start = start - 1, chr = as.character(seqnames)) %>%
  dplyr::select(chr, start, rep)

  
pd <- d50 %>%
  filter(samp == 'GSC') %>% 
  mutate(start = floor(start / bsz) * bsz) %>%
  group_by(chr, start) %>%
  summarise(Laminb1 = mean(Laminb1),
            K9me3 = mean(K9me3)) %>%
  ungroup() %>%
  merge(r) 

r <- sprintf('Pearson\'s r = %s', round(cor(pd$Laminb1, pd$rep, method = 'pearson'), 2))

ggplot(pd, aes(x = Laminb1, y = rep, color = K9me3)) +
  geom_point(alpha = .5, shape = 16) +
  geom_smooth(method = 'lm', color = 'firebrick', size = .5) +
  annotate('label', x = -Inf, y = Inf, hjust = 0, vjust = 1, label = r) +
  scale_color_gradientn('H3K9me3', colors = c(brewer.brbg(21)[21:11],brewer.brbg(11)[5:1]),
                        breaks = c(-1,0)) +
  labs(x = 'Lamin b1 enrichment in GSCs', y = '# of L1/ERV1/ERVK per Mb') +
  facet_wrap(~'Determinants of lamin B1 binding in GSCs') +
  coord_cartesian(clip = F, expand = F) +
  guides(color = guide_colorbar(title.vjust = 1, barwidth = 2, barheight = .5)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 11),
        #legend.position = c(1,0),
        legend.position = 'none',
        legend.justification = c(1,0),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = '#ffffff66'),
        strip.text = element_text(color = 'black', size = 13, face= 'bold'),
        axis.text = element_text(size = 11, family = 'Arial', color = 'black')) +
  ggsave('f5_g.pdf', height = 3.4, width=4, device = cairo_pdf, bg = 'transparent')
