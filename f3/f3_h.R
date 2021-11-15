library(data.table)
library(tidyverse)
library(InteractionSet)
library(rtracklayer)
library(DESeq2)
library(tximport)
library(DiffBind)
library(fgsea)
library(RobustRankAggreg)
library(msigdbr)
library(patchwork)
library(ggdist)
library(ggnewscale)
library(stringr)
library(shadowtext)

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

ps <- msigdbr(species = "Mus musculus") %>% 
  mutate(cat = paste0(gs_cat, '.', gs_subcat) %>% 
             sub('\\.$', '', .)) %>% 
  split(., .$cat) %>% 
  lapply(function(x) split(x$gene_symbol, x$gs_name))

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

vs <- di %>%
  mutate(start = as.numeric(startI) + 2e4,
         end = as.numeric(startJ) - 1e4) %>%
  filter(end > start) %>%
  makeGRangesFromDataFrame() %>%
  subsetByOverlaps(ctcf.dn) %>%
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

vs2 <- di %>%
  mutate(start = as.numeric(startI) + 2e4,
         end = as.numeric(startJ) - 1e4) %>%
  filter(end > start) %>%
  makeGRangesFromDataFrame() %>%
  subsetByOverlaps(ctcf.dn) %>%
  as_tibble() %>%
  mutate(idx = sprintf('%d:%d', start - 2e4, end + 1e4),
         chr = as.character(seqnames)) %>%
  select(chr, idx) %>%
  merge(di) %>%
  merge(ep) %>%
  merge(de, by.x = 'TargetGene', by.y = 'gene_name') %>%
  na.omit() %>%
  group_by(TargetGene) %>%
  summarise(DI.stat = -mean(log2FoldChange),
            DE.stat = -mean(DE.log2FoldChange))

s <- vs %>%
  {list(DI = .[,c('TargetGene', 'DI.stat')],
        DE = .[,c('TargetGene', 'DE.stat')])} %>%
  lapply(function(x) {
    deframe(x) %>% sort() %>% names()
  }) %>%
  {inner_join(aggregateRanks(.),
              aggregateRanks(lapply(., rev)),
              by = 'Name')} %>%
  mutate(sc = log10(Score.x) - log10(Score.y)) %>%
  select(Name, sc) %>%
  deframe()

rr <- names(ps) %>%
  setNames(.,.) %>%
  lapply(function(x) {
  print(x)
  fgsea(ps[[x]], s, eps = 0.0)
})

rrr <- names(ps) %>%
  setNames(.,.) %>%
  lapply(function(x) {
    print(x)
    fgseaRes <- rr[[x]]
    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], ps[[x]], s)
    fgseaRes[pathway %in% collapsedPathways$mainPathways]
  })

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}

ss <- rrr$`C5.GO:BP` %>% 
  mutate(dir = ifelse(ES > 0, 'Up-regulated', 'Down-regulated')) %>% 
  ungroup() %>%
  arrange(NES) %>%
  mutate(pathway = sub('GOBP_', '', pathway) %>%
           gsub('_', ' ', .) %>%
           tolower() %>%
           sub('rna', 'RNA', .) %>%
           sub('nonsense mediated decay', 'NMD', .) %>%
           fct_inorder(),
         symb = as.character(signif.num(padj)))

rcts <- distinct(ss, pathway) %>%
  mutate(alt = as.character((1:n()) %% 2),
         x = as.numeric(pathway))

p1 <- ggplot(ss, aes(x = pathway, y = NES)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 90)) +
  geom_rect(aes(xmin = x - .5, xmax = x + .5, ymin = -Inf, ymax = Inf, fill = alt), 
            data = rcts, show.legend = F, inherit.aes = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  new_scale_fill() +
  geom_col(aes(fill = -log10(padj))) +
  geom_shadowtext(aes(label = symb, y = NES / 2), color = 'black', 
                  bg.colour = 'white', vjust = .8) +
  scale_fill_viridis_c(expression(-log[10]~p[adj]),
                       limits = c(0,max(-log10(ss$padj)))) +
  coord_flip(xlim = c(1, nrow(ss))) +
  scale_y_continuous('Normalized\nenrichment score', breaks = c(-2,0,2)) +
  facet_grid(.~'Enriched \nGO:BP terms') +
  guides(fill = guide_colorbar(barheight = .5, title.vjust = 1)) +
  geom_hline(yintercept = 0, color = 'black') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(1.2,1.2),
        legend.text = element_text(size = 11),
        legend.direction = 'horizontal',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        plot.margin = margin(2,2,10,2),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'))

vss <- split(ss, ss$pathway) %>%
  lapply(function(x) {
    tibble(TargetGene = x$leadingEdge[[1]]) %>%
      merge(vs)
  }) %>%
  bind_rows(.id = 'pathway') %>%
  mutate(pathway = fct_inorder(pathway))

vss <- split(ss, ss$pathway) %>%
  lapply(function(x) {
    tibble(TargetGene = x$leadingEdge[[1]]) %>%
      merge(vs2)
  }) %>%
  bind_rows(.id = 'pathway') %>%
  mutate(pathway = fct_inorder(pathway))

scl <- 10
shf <- .2
dark <- c('#77B3CA', '#C9B97D')
light <- c('#56BBCC', '#BDBC45' )

p2 <- ggplot(vss, aes(x = pathway)) +
  scale_x_discrete() +
  geom_rect(aes(xmin = x - .5, xmax = x + .5, ymin = -Inf, ymax = Inf, fill = alt), 
            data = rcts, show.legend = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  new_scale_fill() +
  geom_hline(yintercept = 0) +
  stat_pointinterval(aes(y = DI.stat * scl), side = 'right', scale = .5,
                     color = light[1], fill = dark[1], position = position_nudge(x = shf),
                     point_size = 2) +
  stat_pointinterval(aes(y = DE.stat), side = 'left', scale = .5,
                     color = light[2], position = position_nudge(x = -shf),
                     point_size = 2) +
  scale_y_continuous(expression(atop(log[2]*'FC (expression)', 'GSC - d4c7PGCLC')),
                     sec.axis = sec_axis(~ (. / 10), expression(log[2]*'FC (interaction)'),
                                         breaks = -1:1), breaks = c(-10,0,10)) +
  coord_flip(xlim = c(1,nrow(ss)), y = c(-10,10), clip = 'on') +
  facet_grid('Leading edge genes' ~ .) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x.top = element_text(color = light[1]),
        axis.title.x.top = element_text(color = light[1]),
        axis.title.x.bottom = element_text(color = light[2]),
        axis.text.x.bottom = element_text(color = light[2]),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'))
{wrap_plots(p1, p2, nrow = 1, widths = c(1.2,1)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('f3_h.pdf', ., height = 5.6, width = 8.5)
