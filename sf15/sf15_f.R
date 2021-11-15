library(tidyverse)
library(LOLA)
library(rtracklayer)
library(ggrepel)
library(patchwork)
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

load('../data/csaw/K9me2.rda')
qSet <- out.ranges %>%
  .[.$direction == 'up' & .$FDR < .05] %>%
  subsetByOverlaps(uSet, .)

res1 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "two.sided")
res2 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "greater")

pd1 <- res1[, setdiff(colnames(res1), c('pValue', 'pValueLog', 'qValue')), with = F] %>%
  merge(res2[,c('userSet', 'dbSet', 'pValue', 'qValue', 'pValueLog')]) %>%
  mutate(reg = sub('.bed', '', filename),
         sig = as.character(signif.num(qValue))) %>%
  filter(support > 100) 

p1 <- ggplot(pd1, aes(x = -log10(qValue), y = oddsRatio, ymin = cLo, ymax = cHi)) +
  geom_linerange() +
  geom_point() +
  geom_label_repel(aes(label = reg), data = ~subset(., qValue < 1e-15),
                   box.padding = 1, nudge_x = 1, color = 'firebrick', alpha = .8) +
  facet_grid(.~'Ensembl regulatory build') +
  labs(x = expression(-'log'[10]~p['adj']), y = 'Overlap odds ratio') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_rect(fill = NA),
        axis.ticks = element_blank()) 

regionDB <- loadRegionDB('../data/resources/ccre')
res1 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "two.sided")
res2 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "greater")

pd2 <- res1[, setdiff(colnames(res1), c('pValue', 'pValueLog', 'qValue')), with = F] %>%
  merge(res2[,c('userSet', 'dbSet', 'pValue', 'qValue', 'pValueLog')]) %>%
  mutate(reg = sub('.bed', '', filename),
         sig = as.character(signif.num(qValue))) 
p2 <- ggplot(pd2, aes(x = -log10(qValue), y = oddsRatio, ymin = cLo, ymax = cHi)) +
  geom_linerange() +
  geom_point() +
  geom_label_repel(aes(label = reg), data = ~subset(., qValue < 1e-2),
                   box.padding = 1, nudge_x = 1, color = 'firebrick', alpha = .8) +
  facet_grid('H3K9me2: GSCLC > GSC regions'~'ENCODE cCREs') +
  labs(x = expression(-'log'[10]~p['adj']), y = 'Overlap odds ratio') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_rect(fill = NA),
        axis.ticks = element_blank()) 

{ wrap_plots(p1, p2, nrow = 1, widths = c(2,1)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf15_f.pdf', ., width = 6.3, height = 4)
