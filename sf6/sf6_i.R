library(tidyverse)
library(pals)
library(GenomicRanges)
library(ggrastr)
library(colorspace)
library(ggnewscale)

load('../data/diffbind/CTCF.GSCvGSCLC.rda')

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)
source('../scripts/diverging_map.R')
cmap <- diverging.colormap(seq(0, 1, .01),
                           rgb1 = hex2RGB(clrs['GSC']),
                           rgb2 = hex2RGB(clrs['GSCLC']),
                           outColorspace = "sRGB") %>%
  {.[. > 1] <- 1; .} %>%
  {rgb(.[,1], .[,2], .[,3])}


anns <- o %>% 
  as_tibble() %>% 
  filter(FDR < .01 & abs(Fold) > 2) %>% 
  mutate(p = Fold > 0) %>%
  count(p) %>%
  mutate(x = ifelse(p, 1, -1) * 5,
         y = 100,
         samp = ifelse(p, 'GSC', 'GSCLC'))

o %>%
  as_tibble() %>%
  ggplot(aes(x = Fold, y = -log10(FDR))) +
  rasterize(geom_point(aes(color = Fold, alpha = -log10(FDR)), shape = 16),
            dpi = 600) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = c(-2, 2)) +
  scale_color_gradientn(expression(log[2]~'(GSC / GSCLC)'),
                        colors = rev(cmap), limits = c(-8,8)) +
  scale_alpha(expression(-log[10]~p[adj])) +
  new_scale_color() +
  geom_label(aes(x = x, y = y, label = n, color = samp), hjust = .5, vjust = .5,
             data = anns, show.legend = F) +
  scale_color_manual(values = clrs) +
  labs(x = expression(log[2]~'(GSC / GSCLC)'), y = expression(-log[10]~p[adj])) +
  facet_wrap(~'Differential CTCF binding') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = 'bold', color = 'black'),
        axis.text = element_text(size =11, color = 'black'),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color ='grey70', linetype = 'dashed')) +
  ggsave('sf6_i.pdf', height = 4, width = 7)

