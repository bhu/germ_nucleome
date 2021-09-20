library(tidyverse)
library(LOLA)
library(rtracklayer)
library(ggrepel)
load('../scripts/runLOLA2.rda')
regionDB <- loadRegionDB('../data/resources/ensembl')
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}
uSet <- read.table('../data/resources/mm10.chrom.sizes', header = F) %>% 
  {Seqinfo(.$V1, .$V2, genome = 'mm10')} %>%
  tileGenome(tilewidth = 1000) %>%
  subsetByOverlaps(import.bed('../data/resources/blacklist.bed'), invert = T) %>%
  GRangesList() %>%
  unlist() %>%
  {.[seqnames(.) != 'chrY']}


qSet <- readRDS('../data/csaw/K27me3.rda') %>%
  .[.$direction == 'up' & .$FDR < .05] %>%
  subsetByOverlaps(uSet, .)

res1 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "two.sided")
res2 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "greater")

res1[, setdiff(colnames(res1), c('pValue', 'pValueLog', 'qValue')), with = F] %>%
  merge(res2[,c('userSet', 'dbSet', 'pValue', 'qValue', 'pValueLog')]) %>%
  mutate(reg = sub('.bed', '', filename),
         sig = as.character(signif.num(qValue))) %>%
  filter(support > 10) %>%
  ggplot(aes(x = -log10(qValue), y = oddsRatio, ymin = cLo, ymax = cHi)) +
  geom_linerange() +
  geom_point() +
  geom_label_repel(aes(label = reg), data = ~subset(., oddsRatio > 2.5),
                   box.padding = 1, 
                   color= 'firebrick', alpha = .8) +
  facet_grid(.~'H3K27me3: GSCLC > GSC') +
  labs(x = expression(-'log'[10]~p['adj']), y = 'Overlap odds ratio') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        axis.ticks = element_blank()) +
  ggsave('sf6_d.pdf', height = 3.7, width = 2.713)
