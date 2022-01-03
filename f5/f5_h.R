library(data.table)
library(tidyverse)
library(rtracklayer)
library(LOLA)
library(ggrepel)

load('../scripts/runLOLA2.rda')
regionDB <- loadRegionDB('../data/resources/reps')
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}

uSet <- fread('../data/resources/mm10.chrom.sizes') %>% 
  {Seqinfo(.$V1, .$V2, genome = 'mm10')} %>%
  tileGenome(tilewidth = 1000) %>%
  subsetByOverlaps(import.bed('../data/resources/blacklist.bed'), invert = T) %>%
  GRangesList() %>%
  unlist() %>%
  {.[seqnames(.) != 'chrY']}

qSet <- list.files('../data/peaks/K9me3', pattern = 'bed$', full.names = T) %>% 
  setNames(.,sub('.bed', '', basename(.))) %>%
  .[c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC')] %>%
  lapply(function(x) {
    fread(x, col.names = c('chr','start','end')) %>%
      mutate(start = start + 1, end = end + 1) %>% 
      makeGRangesFromDataFrame()
  }) %>% 
  Reduce(intersect, .) %>%
  subsetByOverlaps(uSet, .)

res1 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "two.sided")
res2 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "greater")
pd <- res1[,setdiff(colnames(res1), c('pValue', 'pValueLog', 'qValue')),with=F] %>%
  merge(res2[,c('userSet', 'dbSet', 'pValue', 'qValue', 'pValueLog')]) %>%
  mutate(reg = sub('.bed', '', filename),
         sig = as.character(signif.num(qValue))) %>%
  filter(support > 1)

ggplot(pd, aes(x = qValue, y = oddsRatio)) +
  geom_hline(yintercept = 1, color = 'grey70') +
  geom_point() +
  geom_linerange(aes(ymin = cLo, ymax = cHi)) +
  geom_label_repel(aes(label = reg), data = ~subset(., qValue < .05),
                   nudge_x = .3) +
  facet_grid(.~'H3K9me3-enriched\nrepeat families') +
  scale_x_continuous(expression("Fisher's exact test"~p[adj]), breaks = c(0, .5, 1)) +
  #scale_y_continuous(position = 'left') +
  ylab('Overlap odds ratio') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.clip = 'off',
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black')) -> p

ggsave(file = 'f5_h.pdf', height = 3.1, width = 2.3)
