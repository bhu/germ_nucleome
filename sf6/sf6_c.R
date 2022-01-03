library(tidyverse)
library(pals)
library(rtracklayer)

load('../data/tracks/inpnorm.50kb.rda')
load('../data/tracks/mCG.10kb.rda')
lad <- d50 %>% mutate(start = start + 1) %>% 
  distinct(chr, start, end) %>%
  makeGRangesFromDataFrame() %>% 
  overlapsAny(import.bed('../data/PMD/GSC.bed'))

mcg <- o %>%
  mutate(start = round(start / 5e4) * 5e4) %>%
  group_by(chr, start) %>%
  summarise(mCG = mean(GSC), .groups = 'drop') %>%
  na.omit() 

mrks <- names(d50)[-c(17:21)]

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

res <- d50 %>%
  split(., .$samp) %>%
  lapply(function(x) {
    d <- x %>%
      mutate(lad = lad) %>%
      filter(lad) %>%
      merge(mcg)
    mrks %>%
      lapply(function(y) {
        cor(d$mCG, d[[y]], method = 'spearman') %>%
          {tibble(rho = ., mark = y)}
      }) %>%
      bind_rows()
  }) %>%
  bind_rows(.id = 'samp')



mrks <- c('ATAC','CTCF','Rad21',
  'Ring1b','H2Aub','K27me3', 'K27ac', 'K4me1','K4me3',
  'K36me2','K36me3','K9me2','K9me3', 'Laminb1')


res %>%
  mutate(samp = factor(samp, samps),
         mark = factor(mark, mrks)) %>%
  na.omit() %>%
  arrange(samp, mark) %>%
  mutate(mark = sub('ATAC', 'ATAC-seq', mark) %>%
           sub('b1', ' B1', .) %>%
           sub('^K', 'H3K', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x = samp, y = mark)) +
  geom_tile(aes(fill = rho)) +
  scale_fill_gradientn(expression('Spearman\'s'~rho~'vs GSC DNAme in GSC PMDs'),
                       colors = coolwarm(25), limits = c(-.6,.6), breaks = c(-.5, 0, .5)) +
  #guides(fill = guide_colorbar(barheight = .5, barwidth = 4, title.vjust = 1)) +
  guides(fill = guide_colorbar(barwidth = .5, barheight = 2, title.position = 'bottom', title.vjust = 0)) +
  coord_cartesian(expand = F) +
  scale_x_discrete(labels = rnm) +
  #facet_grid(.~'Corr. vs GSC DNAme\n in GSC PMDs') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.background = element_rect(fill = NA),
        axis.title = element_blank(),
        #legend.direction = 'horizontal',
        #legend.justification = c(1,0),
        #legend.position =  c(1,1),
        legend.position = 'right',
        legend.justification = 'top',
        legend.title = element_text(angle = 270),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.clip = 'off',
        axis.text.x = element_text(angle = 30, hjust = 1)) -> p

ggsave(file = 'sf6_c.pdf', p, height = 4.5, width = 3.2, device = cairo_pdf, bg = 'transparent')
  
