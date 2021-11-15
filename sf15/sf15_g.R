library(tidyverse)
library(rtracklayer)
library(rGREAT)
library(patchwork)
load('../data/csaw/K9me2.rda')
db <- out.ranges %>%
  .[.$direction == 'up' & .$FDR < .05]
dels <- import.bed('../data/resources/ccre/ccre/regions/dELS.bed')
ol <- dels %>%
  subsetByOverlaps(db)

job <- submitGreatJob(ol, dels, 'mm10')
et <- getEnrichmentTables(job)

pd <- et$`GO Biological Process` %>%
  slice_min(Hyper_Adjp_BH, n = 10) %>%
  arrange(Hyper_Fold_Enrichment) 

p <- ggplot(pd, aes(x = fct_inorder(name), y = Hyper_Fold_Enrichment, 
             color = -log10(Hyper_Adjp_BH), ymin = 0, ymax = Hyper_Fold_Enrichment)) +
  geom_point(aes(size = Hyper_Foreground_Region_Hits)) +
  geom_linerange() +
  coord_flip() +
  scale_size('|Hits|', breaks= c(100, 300, 500)) +
  scale_color_viridis_c(expression('-log'[10]~'p'['adj']),
                        breaks = c(110,130,150)) +
  facet_grid(.~'Enriched GO:BP terms') +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 95)) +
  guides(color = guide_none()) +
  ylab('Fold enrichment') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.direction = 'horizontal',
        legend.position = c(0,0),
        legend.text = element_text(size = 11),
        legend.justification = c(1.2,1),
        strip.clip = 'off',
        axis.line.x = element_line(color = 'black'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'))

l <- ggplot(pd, aes(x = name, y = Hyper_Fold_Enrichment, color = -log10(Hyper_Adjp_BH))) +
  geom_point() +
  scale_color_viridis_c(expression('-log'[10]~'p'['adj']),
                        breaks = c(110,130,150)) +
  guides(color = guide_colorbar(barheight = 3.1, barwidth = .6)) +
  theme(legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0,1),
        legend.text = element_text(size =11),
        legend.justification = c(0,1))
  

layout <- c(
  area(t = 1, l = 1, b = 20, r = 20),
  area(t = 1, l = 1, b = 1, r = 1)
)
{ {p + cowplot::get_legend(l) + 
  plot_layout(design = layout) } &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf15_g.pdf', ., height = 3.9, width = 7.5)

