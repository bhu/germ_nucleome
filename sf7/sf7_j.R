library(tidyverse)
library(colorspace)
library(rtracklayer)

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)

source('../scripts/diverging_map.R')
cmap <- diverging.colormap(seq(0, 1, .01),
                           rgb2 = hex2RGB(clrs['GSC']),
                           rgb1 = hex2RGB(clrs['GSCLC']),
                           outColorspace = "sRGB") %>%
  {.[. > 1] <- 1; .} %>%
  {rgb(.[,1], .[,2], .[,3])}

load('../data/tracks/inpnorm.50kb.rda')
load('../data/diffbind/CTCF.GSCvGSCLC.rda')

db <- o[o$FDR < 1e-3 & o$Fold < -3] 
dat <- cbind(d50[d50$samp == 'GSC',c('chr','start','end')],
            d50[d50$samp == 'GSC',-c(17:21)] - 
              d50[d50$samp == 'GSCLC',-c(17:21)]) %>%
  mutate(start = start + 1) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

hits <- findOverlaps(db, dat) %>% as("List")

mcols(dat) %>% 
  colnames() %>% 
  setNames(.,.) %>% 
  lapply(function(x) {extractList(mcols(dat)[[x]], hits) %>% mean()}) %>%
  bind_cols() %>%
  mutate(idx = 1:n()) %>% 
  pivot_longer(-idx, names_to = 'mark', values_to = 'v') %>%
  group_by(mark) %>% 
  mutate(m = median(v, na.rm = T)) %>%  
  ungroup() %>% 
  arrange(m) %>% 
  filter(!grepl('Stag', mark)) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>% fct_inorder()) %>% 
  ggplot(aes(x = mark, y = -v, color = m, fill = m)) + 
  geom_hline(yintercept = 0, color = 'black') +
  geom_boxplot(outlier.color = NA) +
  stat_summary(geom = "crossbar", width = 0.4, fatten = 0, color = "white",
              fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  coord_cartesian(ylim = c(-.4,.4)) +
  scale_fill_gradientn(colors = cmap, limits = c(-.1,.1)) +
  scale_color_gradientn(colors = cmap, limits = c(-.1,.1)) +
  ylab(expression(log[2]~frac(GSCLC,GSC)~'at GSCLC > GSC CTCF peaks')) +
  facet_wrap(.~'Epigenetic determinats of retained CTCF binding') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = 'bold', color = 'black'),
        axis.text = element_text(size =11, color = 'black'),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none',
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color ='grey70', linetype = 'dashed'))  -> p

ggsave('sf7_j.pdf', p, height = 4, width = 7)