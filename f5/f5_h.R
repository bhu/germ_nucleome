library(tidyverse)
library(tximport)
library(rtracklayer)
library(pals)
library(ggh4x)
library(ggtext)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC') 
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
samps <- rnm(samps)
clrs <- setNames(pals::tableau20(9)[seq(1, 9, 2)], samps)

g <- import.gff3('../data/resources/gencode.vM25.annotation.gff3.gz')
e2g <- g[g$type == 'gene'] %>%
  {setNames(.$gene_name, .$gene_id)}
tx2gene <- g[g$type == 'transcript'] %>%
  as_tibble() %>%
  dplyr::select(transcript_id, gene_id) %>%
  distinct() 

fs <- list.files('../data/rna', full.names = T) %>%
  {setNames(file.path(., 'quant.sf'), basename(.))}
mdat <- data.frame(s = names(fs)) %>%
  separate(s, c('type', 'rep'), '_', remove = F) %>%
  column_to_rownames('s')
salmon <- tximport(fs, type = "salmon", tx2gene = tx2gene)

scls <- lapply(1:5, function(x) {
  scale_x_continuous(breaks = 1:5, labels = case_when(x == 1 ~samps, T ~ rep('', 5)))
})

salmon$abundance %>% 
  as.data.frame() %>%
  rownames_to_column('name') %>%
  mutate(name = e2g[name]) %>%
  filter(name %in% c('Suv39h1', 'Suv39h2', 'Setdb1', 'Ehmt1', 'Ehmt2')) %>%
  pivot_longer(-name, names_to = 'samp', values_to = 'tpm') %>%
  mutate(type = sub('_.*', '', samp) %>%
           rnm() %>%
           factor(samps),
         name = sub('Ehmt1', 'GLP', name) %>%
           sub('Ehmt2', 'G9a', .) %>%
           factor(., c('Suv39h1', 'Suv39h2', 'Setdb1','GLP', 'G9a')),
         x = as.numeric(type)) %>%
  na.omit() %>%
  ggplot(aes(x = x, y = log10(tpm), color = type)) +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), color = 'grey70', alpha = .5) +
  geom_point() +
  facet_grid(.~name, scales = 'free_x') +
  ylab('RPM') +
  scale_color_manual(values = clrs) +
  facetted_pos_scales(x = scls) +
  scale_y_continuous(breaks = log10(c(10,30,100)),
                     labels = c(10,30,100)) +
  scale_x_continuous(labels = samps) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black', size = 11),
        #axis.text.x = element_markdown(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold')) + 
  ggsave('f5_h.pdf', height = 2.55, width = 5.5)
 