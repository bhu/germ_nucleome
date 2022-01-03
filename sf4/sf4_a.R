library(data.table)
library(tidyverse)
library(pals)
library(LOLA)
library(GenomicRanges)

odr <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
odrr <- c('repressed', 'CTCF', 'enhancer',  'bivalent', 'promoterCTCF', 'promoter')
clrs <- setNames(tableau20(12)[seq(1, 10, 2)], odr)
cclrs <- setNames(kelly(22)[2:7], odrr)

load('../scripts/runLOLA2.rda')
regionDB <- loadRegionDB('../data/resources/ensembl')
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}


load('../data/clust/open.rda')

qSet <- d %>%
  filter(cl != -1 & type %in% odr) %>% 
  mutate(V1 = cl,
         cl = case_when(V1 %in% c(1,10,11) ~ 'promoterCTCF',
                        V1 %in% c(7,8,9,12) ~ 'bivalent',
                        V1 %in% c(6,13) ~ 'promoter',
                        V1 %in% c(4,5,14,17) ~ 'enhancer',
                        V1 %in% c(15,16,18) ~ 'repressed',
                        V1 %in% c(0,2,3) ~ 'CTCF',
                        T ~ NA_character_) %>%
           factor(odrr)) %>% 
  na.omit() %>%
  distinct(cl, crd) %>%
  separate(crd, c('chr','start','end'), '[:-]') %>% split(., .$cl) %>% 
  lapply(makeGRangesFromDataFrame) %>%
  GRangesList()
  
uSet <- unlist(qSet) %>% unique()
  
res1 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "two.sided")
res2 <- runLOLA2(qSet, uSet, regionDB, cores = 4,
                 direction = "greater")

pd <- res1[,setdiff(colnames(res1), c('pValue', 'pValueLog', 'qValue')),with=F] %>%
  merge(res2[,c('userSet', 'dbSet', 'pValue', 'qValue', 'pValueLog')]) %>%
  {.[.$filename %in% .$filename[.$support > 300 & .$qValue < .25],]} %>%
  mutate(reg = sub('.bed', '', filename),
         sig = as.character(signif.num(qValue))) %>%
  dplyr::rename(oddsLower = cLo, oddsUpper = cHi) %>%
  group_by(reg) %>%
  mutate(msig = mean(-log10(qValue))) %>%
  ungroup() %>%
  arrange(msig) %>%
  mutate(userSet = factor(userSet),
         reg = fct_inorder(reg),
         y = as.numeric(reg),
         x = as.numeric(userSet),
         z = log2(oddsRatio),
         z = case_when(!is.finite(z) ~ -7.5, T ~ z)) 
xl <- distinct(pd, userSet, x)
yl <- distinct(pd, reg, y)
ggplot(pd, aes(x = y, y = x, fill = z)) +
  geom_tile() +
  geom_text(aes(label = sig)) +
  scale_fill_gradientn('log2(odds ratio)', limits = c(-8.2, 8.2), colors = rev(brewer.rdgy(25))) +
  coord_cartesian(expand = F) +
  scale_x_continuous('Annotation', breaks = yl$y, labels = yl$reg) +
  scale_y_continuous('Open site cluster', breaks = xl$x, labels = xl$userSet) +
  facet_grid(.~'Genomic localization of open site clusters') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.clip = 'off',
        legend.key = element_blank(),
        legend.position = c(1,0),
        legend.text = element_text(size = 11),
        legend.direction = 'horizontal',
        legend.justification = c(1,1.8)) +
  guides(fill = guide_colorbar(barheight = .5, title.position = 'top')) -> p
  ggsave('sf3_a.pdf', p, height = 3.7, width = 6.3)
