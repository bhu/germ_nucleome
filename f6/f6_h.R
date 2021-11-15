library(tidyverse)
library(pals)
library(colorspace)
library(patchwork)
source('../scripts/diverging_map.R')

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)
cmap <- diverging.colormap(seq(0, 1, .01),
                           rgb1 = hex2RGB(clrs['GSC']),
                           rgb2 = hex2RGB(clrs['GSCLC']),
                           outColorspace = "sRGB") %>%
  {.[. > 1] <- 1; .} %>%
  {rgb(.[,1], .[,2], .[,3])}


load('../data/tads/consensus.tad.pileup.rda')

pd <- lapply(d[c('GSCLC','GSC')], function(z) {
  z %>%
    rownames_to_column('y') %>%
    pivot_longer(-y, names_to = 'x', values_to = 'z') %>%
    mutate(y = as.numeric(y), x = as.numeric(sub('^X', '', x)))
}) %>%
  bind_rows(.id = 'samp') %>%
  pivot_wider(names_from = 'samp', values_from = 'z') %>%
  mutate(z = log2(GSCLC) - log2(GSC)) 
 
m <- max(abs(quantile(pd$z, .001, .999)))

p <- pd %>%
  mutate(z = case_when(abs(x - y) < 3 ~ NA_real_, T ~ z)) %>%
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_y_reverse(breaks= c(34,67)) +
  scale_fill_gradientn('log2(GSCLC / GSC)', colors = cmap, limits = c(-m2, m2),
                       breaks = c(-.1,0,.1), oob = scales::squish) +
  facet_grid('Difference'~ .) +
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = c(33, 66)) +
  theme(legend.position = 'none',
        legend.justification = 'bottom',
        #legend.position = c(1,0),
        #legend.justification = c(1,1.2),
        #legend.direction = 'horizontal',
        axis.title = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        # plot.margin = margin(l = 5),
        plot.margin = margin(t = 0),
        axis.text = element_blank(),
        legend.text = element_text(size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color= 'black', size = 13, face = 'bold'),
        legend.title = element_text(angle = 90)) +
  guides(fill = guide_colorbar(barwidth = .5, barheight = 6.7, title.position = 'left')) -> p
ggsave('f6_h.pdf', p, height = 2, width = 2)

pd %>%
  ggplot(aes(x, y, fill = z)) +
  geom_raster() +
  scale_y_reverse(breaks= c(34,67)) +
  scale_fill_gradientn(expression(log[2]~frac('GSCLC','GSC')), 
                       colors = cmap, limits = c(-m2, m2),
                       breaks = c(-.1,0,.1), oob = scales::squish) +
  guides(fill = guide_colorbar(barheight  = .5, barwidth = 4.4, title.position = 'top',
                               title.hjust = .5)) +
  theme(legend.background = element_blank(),
        plot.background = element_blank(),
        legend.position = 'bottom') -> l
ggsave('f6_h.leg.pdf', l, width = 2)

