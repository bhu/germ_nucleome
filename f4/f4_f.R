library(tidyverse)
library(pals)
library(rtracklayer)
library(UpSetRm)
library(patchwork)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

uSet <- read.table('../data/resources/mm10.chrom.sizes', header = F) %>% 
  {Seqinfo(.$V1, .$V2, genome = 'mm10')} %>%
  tileGenome(tilewidth = 1000) %>%
  subsetByOverlaps(import.bed('../data/resources/blacklist.bed'), invert = T) %>%
  GRangesList() %>%
  unlist() %>%
  {.[seqnames(.) != 'chrY']}

qSet <- list.files('../data/peaks/LADs', full.names = T)  %>% 
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(function(x) {
    import.bed(x) %>%
      subsetByOverlaps(uSet, .)
  })

ol <- tibble(samp = names(qSet)) %>%
  mutate(study = case_when(grepl('Poleshko', samp) ~ 'Poleshko 2017',
                         grepl('BDF121|AAG', samp) ~ 'This study',
                         grepl('OPC', samp) ~ 'Yattah 2020',
                         grepl('diffe', samp) ~ 'Robson 2016',
                         T ~ 'Peric-Hupkes 2010') %>%
         factor(rev(c('This study', 'Yattah 2020', 'Poleshko 2017', 'Robson 2016', 'Peric-Hupkes 2010')))) %>%
  split(., .$study) %>%
  lapply(function(x) {
    qSet[x$samp] %>%
      GRangesList() %>%
      unlist() %>%
      unique() %>%
      overlapsAny(uSet, .)
  }) 


ps <- lapply(ol, which) %>%
  fromList() %>%
  upset(order.by = 'freq', show.numbers = 'no', ncut = 10)

p1 <- ps$Main_bar + theme_gray() +
  scale_y_continuous('Set coverage', breaks = seq(0,8e5,2e5),
                     labels = c('0', '200Mb', '400Mb', '600Mb', '800Mb')) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        plot.margin = margin(0),
        panel.grid.major.y = element_line(color=  'grey70', linetype = 'dashed')) 

p3 <- ps$Matrix + theme_gray() +
  xlab('Intersection') +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(hjust = .5),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0))

p2 <- ps$Sizes+ theme_gray() +
  scale_y_reverse('Total coverage', breaks = c(0,1e6),
                     labels = c(0,'1Gb')) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0))

wrap_plots(plot_spacer() + theme_void(),
           p1 + theme(axis.title.y = element_text(margin = margin(r = -80, unit = 'pt'))),
           p2, p3, nrow = 2, heights = c(2,1), widths = c(1,2)) &
  theme(plot.background = element_blank()) &
  ggsave('f4_f.pdf', height = 2.85, width = 5.7)

