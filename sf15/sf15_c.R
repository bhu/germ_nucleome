library(data.table)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(pals)
library(patchwork)
library(ggnewscale)
library(shadowtext)

load('../data/clust/promoter.rda')
s <- d %>%
  filter(type %in% c('GSCLC', 'GSC')) %>%
  filter(gene_type == 'protein_coding') %>%
  select(type, gene_name, K27me3) %>%
  group_by(type, gene_name) %>%
  summarise(K27me3 = mean(K27me3)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'type', values_from = 'K27me3') %>%
  mutate(dif = GSCLC - GSC) %>%
  select(gene_name, dif) %>%
  na.omit() %>%
  deframe()

ps <- msigdbr(species = "Mus musculus") %>% 
  mutate(cat = paste0(gs_cat, '.', gs_subcat) %>% 
           sub('\\.$', '', .)) %>% 
  split(., .$cat) %>% 
  lapply(function(x) split(x$gene_symbol, x$gs_name))

rr <- names(ps) %>%
  setNames(.,.) %>%
  lapply(function(x) {
    print(x)
    fgsea(ps[[x]], s, eps = 0.0)
  })

r <- fgsea(ps$`C5.GO:BP`, s, eps = 0.0)
clps <- collapsePathways(r[padj < .05], ps$`C5.GO:BP`, s)
r2 <- r[pathway %in% clps$mainPathways][order(-NES),]


signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}

ss <- r2 %>% 
  filter(NES > 0) %>% 
  arrange(padj) %>%
  head(15) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(pway = pathway,
         pathway = sub('GOBP_', '', pathway) %>%
           gsub('_', ' ', .) %>%
           tolower() %>%
           sub('dna', 'DNA', .) %>%
           sub('rna', 'RNA', .) %>%
           sub(' i ', ' I ', .) %>%
           fct_inorder(),
         symb = as.character(signif.num(padj)))

rcts <- distinct(ss, pathway) %>%
  mutate(alt = as.character((1:n()) %% 2),
         x = as.numeric(pathway))

ggplot(ss, aes(x = pathway, y = NES)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 90)) +
  geom_rect(aes(xmin = x - .5, xmax = x + .5, ymin = -Inf, ymax = Inf, fill = alt), 
            data = rcts, show.legend = F, inherit.aes = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  new_scale_fill() +
  geom_col(aes(fill = -log10(padj)), width = .7) +
  scale_fill_viridis_c(expression('-log'[10]~'p'['adj']),
                       limits = c(0,max(-log10(ss$padj))), breaks = c(0, 6)) +
  coord_flip(xlim = c(1, 15), clip = F) +
  geom_shadowtext(aes(label = symb, y = NES / 2), color = 'black', 
                  bg.colour = 'white', vjust = .8)  +
  scale_y_continuous('Normalized enrichment score', expand = expansion(c(0, .05)), breaks = 0:2) +
  facet_grid(.~'Enriched\nGO:BP pathways') +
  guides(fill = guide_colorbar(barheight = .5, title.vjust = 1, barwidth = 3)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.clip = 'off',
        axis.text = element_text(color = 'black', size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        legend.justification = c(1,0.2),
        legend.text = element_text(size = 11),
        strip.background = element_rect(fill = NA),
        legend.direction = 'horizontal',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        plot.margin = margin(5,30,5,5),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'))  -> p
ggsave('sf15_c.pdf', p, height = 3.8, width = 5)
