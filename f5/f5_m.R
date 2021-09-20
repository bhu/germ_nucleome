library(data.table)
library(tidyverse)
library(rtracklayer)
library(ggdist)
library(ggpubr)
library(gghalves)
library(pals)
library(patchwork)
library(ggrastr)

g <- import.gff3('../data/resources/gencode.vM25.basic.protein_coding.annotation.gff3.gz')
gg <- readLines('../data/resources/Kurimoto_germline_genes.txt')

gs <- list.files('../data/resources/', pattern = 'QuickGO', full.names = T) %>%
  lapply(function(x) fread(x)$SYMBOL) %>%
  unlist() %>%
  unique()
ko <- readLines("../data/resources/setdb1_KO_up.txt") %>%
  unique() %>%
  setdiff('') %>%
  grep('^NA$', ., invert = T, value = T)

sg <- intersect(ko, gs) 

load('../data/tss.1kb.MSnorm.rda')

s <- tss.1kb.norm.mass[c('d2PGCLC', 'd4c7PGCLC')] %>% 
  bind_rows(.id = 'samp') %>%
  filter(name %in% unique(g$gene_name)) %>%
  select(samp, name, K9me3) %>%
  group_by(samp, name) %>%
  summarise(K9me3 = mean(K9me3)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'samp', values_from = 'K9me3') %>%
  mutate(dif = d2PGCLC - d4c7PGCLC) %>%
  select(name, dif) %>%
  na.omit() %>%
  {rbind(mutate(.[.$name %in% gg,], cell = 'Germline genes\n(Kurimoto)'),
         mutate(.[.$name %in% sg,], cell = 'Setdb1-repressed\ngermline genes'),
         mutate(.[!(.$name %in% c(gg, sg)),], cell = 'Other'))} %>%
  mutate(cell = fct_inorder(cell)) 

stat.dif <- compare_means(dif ~ cell, data = s, method = 'wilcox.test') %>%
  mutate(y.position = 4 + (1:n()) * 1.4)

ggplot(s, aes(x = cell, y = -dif)) + 
  geom_hline(yintercept = 0, color= 'grey70') +
  rasterize(geom_half_point(aes(color = cell), range_scale = 0.4, size = .1, side = 'l',
                  position = position_nudge(x = -.05)), dpi = 600) +
  stat_slab(aes(fill = cell), alpha = .5, scale = .6) +
  stat_pointinterval(aes(color = cell)) +
  stat_pvalue_manual(stat.dif, coord.flip = T, label = "p.signif",
                     tip.length = .01) +
  scale_y_continuous(expand = expansion(c(0, .05)),
                     breaks= c(-5,0,5)) +
  ylab(expression("Promoter H3K9me3 log"[2]~"(d4c7 / d2 mPGCLC)")) +
  facet_grid(.~'\u0394H3K9me3 at germline genes') +
  coord_flip(clip = 'off') +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank()) +
  ggsave('f5_m.pdf', device = cairo_pdf, bg = 'transparent', height = 3.68, width = 5)
