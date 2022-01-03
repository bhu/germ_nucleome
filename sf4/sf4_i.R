library(data.table)
library(tidyverse)
library(InteractionSet)
library(rtracklayer)
library(DESeq2)
library(tximport)
library(patchwork)
library(RRHO2)

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

dds <- DESeqDataSetFromTximport(txi = salmon,
                                colData = mdat,
                                design = ~type)
dds <- DESeq(dds)

de <- results(dds, contrast = c('type', 'd4c7PGCLC', 'GSC')) %>%
  as.data.frame() %>%
  rownames_to_column('ID') %>%
  mutate(gene_name = e2g[ID]) %>%
  select(gene_name, DE.pvalue = pvalue, DE.padj = padj, DE.stat = stat,
         DE.baseMean = baseMean, DE.log2FoldChange = log2FoldChange) %>%
  na.omit()

load('../data/abc/diffint.d4c7vGSC.rda')
ctcf.dn <- import.bed('../data/peaks/CTCF/d4c7vGSC.d4c7_up.bed')

round_any = function(x, accuracy, f = round) f(x / accuracy) * accuracy
ep <- c('GSC','d4c7PGCLC') %>%
  sprintf('../data/abc/%s.csv.gz', .) %>%
  lapply(fread) %>%
  rbindlist() %>% 
  mutate(i = round_any((start + end) / 2, 1e4, floor), 
         j = round_any(TargetGeneTSS, 1e4, floor), 
         startI = case_when(i > j ~ j, T ~ i), 
         startJ = case_when(i > j ~ i, T ~ j)) %>%
  select(chr, startI, startJ, TargetGene) %>% 
  distinct() %>% 
  arrange(chr, startI, startJ) %>% 
  filter(startI != startJ)

pd <- c('E-P pairs straddling \nGSC-low CTCF sites' = F, 'Other E-P pairs' = T) %>%
  lapply(function(x) {
    vs <- di %>%
      mutate(start = as.numeric(startI) + 2e4,
             end = as.numeric(startJ) - 1e4) %>%
      filter(end > start) %>%
      makeGRangesFromDataFrame() %>%
      subsetByOverlaps(ctcf.dn, invert = x) %>%
      as_tibble() %>%
      mutate(idx = sprintf('%d:%d', start - 2e4, end + 1e4),
             chr = as.character(seqnames)) %>%
      select(chr, idx) %>%
      merge(di) %>%
      merge(ep) %>%
      merge(de, by.x = 'TargetGene', by.y = 'gene_name') %>%
      na.omit() %>%
      group_by(TargetGene) %>%
      summarise(DI.stat = -mean(stat),
                DE.stat = -mean(DE.stat))
    v1 <- as.data.frame(vs[,1:2])
    v2 <- as.data.frame(vs[,c(1,3)])
    
    RRHO_obj <-  RRHO2_initialize(v1, v2, labels = c("di", "de"), log10.ind=TRUE,
                                  boundary = .02, multipleTesting = 'BH')
    
    RRHO_obj$hypermat %>%
      as.data.frame() %>%
      rownames_to_column('x') %>%
      pivot_longer(-x, names_to = 'y', values_to = 'z') %>%
      mutate(x = as.numeric(x),
             y = as.numeric(sub('^V', '', y)))
  }) %>%
  bind_rows(.id = 'kind')


pd %>%
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_tile() +
  #annotate('label', x = c(-Inf, Inf), y = c(-Inf, Inf), hjust = c(-.01, 1.01), vjust = c(-.1, 1.1),
  #         label = c('d4c7PGCLC > GSC', 'd4c7PGCLC < GSC'), alpha = .8) +
  scale_fill_viridis_c(name = expression(-log[10]~p[adj]), 
                       na.value = '#00000000', breaks = c(0,10,20)) +
  labs(x = 'Differential E-P interaction', y = 'Differential expression') +
  facet_wrap(~kind, scales = 'free', strip.position = 'right', ncol=1)+
  coord_cartesian(expand = F) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(angle = 270),
        #legend.position = 'right',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.clip = 'off',
        #legend.direction = 'horizontal',
        strip.text = element_text(size = 13, color = 'black', face = 'bold'),
        legend.key = element_blank()) +
  guides(fill = guide_colorbar(barwidth = .5, barheight = 4, title.position = 'top', title.vjust = 0)) -> p
  ggsave('sf4_i.pdf', p, height = 3.8, width = 3.1) 
