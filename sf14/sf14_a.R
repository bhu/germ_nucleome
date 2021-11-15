library(tidyverse)
library(rtracklayer)
library(pals)
library(colorspace)
library(tximport)

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)

g <- import.gff3('../data/resources/gencode.vM25.annotation.gff3.gz')
gg <- g[g$type == 'gene']
e2g <- gg %>%
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


source('../scripts/diverging_map.R')
cmap <- diverging.colormap(seq(0, 1, .01),
                           rgb1 = hex2RGB(clrs['GSC']),
                           rgb2 = hex2RGB(clrs['GSCLC']),
                           outColorspace = "sRGB") %>%
  {.[. > 1] <- 1; .} %>%
  {rgb(.[,1], .[,2], .[,3])}

pd <- salmon$abundance %>% 
  as.data.frame() %>%
  rownames_to_column('ID') %>% 
  select(contains('GSC'), ID) %>% 
  filter(ID %in% gg$gene_id[gg$gene_type == 'protein_coding']) %>%
  pivot_longer(-ID, names_to = 'samp', values_to = 'TPM') %>% 
  mutate(samp = sub('_.*', '', samp)) %>%
  group_by(samp, ID) %>% 
  summarise(TPM = mean(TPM)) %>% 
  pivot_wider(names_from = 'samp', values_from = 'TPM') %>%
  mutate(d = log2(GSCLC / GSC)) 

rho <- sprintf('Spearman\'s \u03c1: %.3f', cor(pd$GSC, pd$GSCLC, method = 'spearman')) 

mx <- max(abs(pd$d[is.finite(pd$d)]), na.rm = T)

ggplot(pd, aes(x = GSC, y = GSCLC, color = d)) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = .5, size = .5, shape = 16) + 
  scale_color_gradientn(expression(log[2]~'GSCLC/GSC'), colors = cmap,
                        limits = c(-mx, mx), breaks = c(-5, 0, 5)) +
  scale_x_continuous(trans = 'log1p', breaks = 10^(1:4)) +
  scale_y_continuous(trans = 'log1p', breaks = 10^(1:4)) +
  coord_cartesian(expand = F, clip = F) +
  annotate('label', x = 0, y = max(pd$GSCLC), hjust = -0.01, vjust = 1.1,
           label = rho, alpha = .7) +
  guides(color = guide_colorbar(title.position = 'top', barheight = .5)) +
  facet_grid(.~'Gene expression (RPM)') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.direction = 'horizontal',
        legend.key = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.background = element_rect(fill = '#ffffff66', color = NA),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(color = 'black', size = 11)) -> p
ggsave('sf14_a.pdf', height = 3.6, width = 3.8, device = cairo_pdf, bg = 'transparent')

