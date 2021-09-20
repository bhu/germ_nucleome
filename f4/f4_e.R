library(tidyverse)
library(pals)
library(rtracklayer)
library(eulerr)

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


fit <- list('All other studies' = ol[names(ol) != 'This study'] %>%
              bind_cols() %>% 
              apply(1, any),
            'ESC/EpiLC/d2/d4c7PGCLC' = qSet[grep('BDF121', names(qSet))] %>%
              GRangesList() %>%
              unlist() %>%
              unique() %>%
              overlapsAny(uSet, .),
            'GSC' = overlapsAny(uSet, qSet$GSC_AAG)) %>%
  bind_cols() %>%
  euler(shape = 'ellipse')

pdf('f4_e.pdf', height = 3.3, width = 5.5)
plot(fit, fills = c('grey70', 'white', clrs['GSC']),
     quantities = sprintf('%.0fMb', fit$original.values / 1000))
dev.off()
