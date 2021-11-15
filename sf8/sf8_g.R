library(data.table)
library(tidyverse)
library(InteractionSet)
library(rtracklayer)
library(DESeq2)
library(tximport)
library(patchwork)
library(RRHO2)
library(pals)

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

pds <- list('ctcf.dn' = F, 'other' = T) %>%
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
    
    pd <- RRHO_obj$hypermat %>%
      as.data.frame() %>%
      rownames_to_column('x') %>%
      pivot_longer(-x, names_to = 'y', values_to = 'z') %>%
      mutate(x = as.numeric(x),
             y = as.numeric(sub('^V', '', y)))
    nm <- max(pd$x)
    vss <- lapply(list(DI = v1, DE = v2), function(x) {
      x %>% 
        `colnames<-`(c('gene', 'stat')) %>%
        mutate(num = ntile(-stat, n = nm)) %>%
        group_by(num) %>%
        summarise(stat = mean(stat, na.rm = T))
    }) 
    list(pd = pd, vss = vss)
  })



pl <- ggplot(pds$other$vss$DE, aes(x = num, fill = stat, y = 1)) + 
  geom_tile() + 
  scale_fill_gradient2(low = clrs[2], 
                       high = clrs[1], 
                       limits = c(-10,10), 
                       oob = scales::squish) +
  coord_flip(expand = F) +
  xlab('Genes ranked by differential expression') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none',
        plot.margin = margin(0,0,0,0),
        axis.ticks = element_blank())

pb <- ggplot(pds$other$vss$DI, aes(x = num, fill = stat, y = 1)) + 
  geom_tile() + 
  scale_fill_gradient2(low = clrs[2], 
                       high = clrs[1], 
                       limits = c(-1,1), 
                       oob = scales::squish) +
  coord_cartesian(expand = F) +
  xlab( 'Genes ranked by differential E-P interaction') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none',
        plot.margin = margin(0,0,0,0),
        axis.ticks = element_blank())

pm <- ggplot(pds$other$pd, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  #annotate('label', x = c(-Inf, Inf), y = c(-Inf, Inf), hjust = c(-.01, 1.01), vjust = c(-.1, 1.1),
  #         label = c('d4c7PGCLC > GSC', 'd4c7PGCLC < GSC'), alpha = .8) +
  scale_fill_viridis_c(name = expression(-log[10]~p[adj]~'of overlap'), 
                       na.value = '#00000000', breaks = c(0,10,20),
                       limits = range(pds$ctcf.dn$pd$z, na.rm = T)) +
  coord_cartesian(expand = F) +
  facet_wrap(~'E-P not straddling\nGSC-depleted CTCF binding sites') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
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
  guides(fill = guide_colorbar(barwidth = .5, barheight = 4, title.position = 'top', title.vjust = 0))  


{ wrap_plots(pl, pm, plot_spacer(), pb, heights = c(20, 1),
           widths = c(1,20), nrow = 2) &
  theme(plot.background = element_blank()) } %>%
ggsave('sf8_g.pdf', ., width = 4.85, height = 4.5) 

