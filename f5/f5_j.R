library(tidyverse)
library(pals)
library(ggrastr)

load('../data/tracks/inpnorm.50kb.rda')
d50 %>%
  filter(samp == 'd4c7PGCLC') %>%
  dplyr::select(K9me3, K27me3, K9me2, K36me2) %>%
  mutate(EpiLC = d50$K36me2[d50$samp == 'EpiLC'],
         dif = EpiLC - K36me2) %>%
  ggplot(aes(x = K9me3, y = K27me3, color = EpiLC - K36me2)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(alpha = .5, size = .25,  shape = 16) +
  scale_color_gradientn(expression(atop(Delta*'H3K36me2', log[2]~frac('EpiLC', 'd4c7 mPGCLC'))) ,
                        colors = rev(brewer.piyg(10)), 
                        limits = c(-1.3,1.3), breaks= c(-1,0,1)) +
  labs(x = 'H3K9me3 in d4c7 mPGCLCs', y = 'H3K27me3 in d4c7 mPGCLCs') +
  guides(color = guide_colorbar(title.position = 'top',barwidth = 5.5,
                                barheight = .5, title.hjust = .5)) +
  facet_grid(. ~'Determinants of H3K36me2 loss') +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.direction = 'horizontal',
        legend.key = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.clip = 'off',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.background = element_rect(fill = '#ffffff66', color = NA),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        panel.background = element_blank()) -> p
ggsave('f5_j.pdf', p, height = 4.1, width = 3.3, device = cairo_pdf, bg = 'transparent')
