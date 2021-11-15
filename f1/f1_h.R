library(data.table)
library(tidyverse)
library(tidygraph)
library(igraph)
library(GenomicRanges)
library(boot)
library(pals)
library(ggraph)
library(patchwork)

load('../data/tads/insulscore.rda')
load('../data/compscore/25kb.tads.rda')
load('../data/cliques/cliq.int.rda')

samps <- c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
           'ESC (serum)', 'NPC', 'CN',
           'Day 6','Day 4','Day 2','B\u03B1',
           'CM', 'CPC', 'PCM', 'VCM')
clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(12)[-4]), samps)

rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}

res <- list(germ = c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC'),
            neural = c('ESC', 'NPC', 'CN'),
            stadhouders = c('D6', 'D4', 'D2', 'Ba'),
            cardio = c('CM', 'CPC', 'PCM', 'VCM')) %>%
  lapply(function(x) setNames(x, x)) %>%
  Map(function(ss, p) {
    tad <- tads[[p]] %>%
      `colnames<-`(c('chr','start','end')) %>%
      mutate(start = start + 1) %>%
      makeGRangesFromDataFrame()
    nds <- as_tibble(tad)[,1:3] %>%
      mutate(idx = sprintf('%s:%d-%d', seqnames, start-1, end)) %>%
      pull(idx)
    gs <- ss %>%
      lapply(function(s) {
        print(s)
        gr <- na.omit(eig[[p]][[s]]) %>%
          mutate(comp = case_when(E1 > 0 ~ 'A', T ~ 'B'),
                 start = start + 1) %>%
          dplyr::select(chrom, start, end, E1, comp) %>%
          makeGRangesFromDataFrame(keep.extra.columns = T)
        ab <- findOverlaps(tad, gr) %>%
          as("List") %>%
          extractList(gr$comp, .) %>%
          sapply(function(y) {
            ifelse(sum(y == 'A') > sum(y == 'B'), 'A', 'B')
          }) %>%
          setNames(nds)
        ints <-  o[[p]][[s]] %>%
          mutate(i = sprintf('%s:%d-%d', V1, V2, V3),
                 j = sprintf('%s:%d-%d', V4, V5, V6)) %>%
          select(i, j)
        g <- ints %>%
          as.matrix() %>%
          graph_from_edgelist() %>%
          max_cliques() %>%
          lapply(function(y) {
            tibble(node = names(y)) %>%
              add_count()
          }) %>%
          bind_rows(.id = 'idx') %>%
          filter(n > 2) %>% 
          separate(node, c('V1', 'V2', 'V3'), '[:-]', F) %>%
          arrange(idx, V1, V2) %>%
          select(idx, node) %>%
          setDT() %>%
          .[, {tmp <- combn(node, 2); .(i = tmp[1,], j = tmp[2,])}, idx] %>%
          select(-idx) %>%
          merge(ints) %>%
          distinct() %>%
          as_tbl_graph() 
        cs <- igraph::components(g) %>% {.$csize[.$membership]} 
        s <- rnm(s)
        g %>%
          activate(edges) %>%
          mutate(cs = cs[from],
                 smp = s) %>%
          activate(nodes) %>%
          mutate(cs = cs,
                 smp = s,
                 comp = ab[name]) 
      })
    gs
  }, ., names(.))


c.b <- 'firebrick3'
c.a <- 'darkolivegreen4'

ps <- res[c('germ','cardio')] %>%
  unname() %>%
  unlist(recursive = F) %>%
  Map(function(x, s) {
    p <- x %>%
      activate(edges) %>%
      mutate(proj = ifelse(s == 'ESC', 'Germline', 'Cardiac')) %>%
      tidygraph::filter(cs > 4) %>%
      activate(nodes) %>%
      mutate(proj = ifelse(s == 'ESC', 'Germline', 'Cardiac')) %>%
      tidygraph::filter(cs > 4) %>%
      ggraph(layout = 'kk') +
      scale_color_manual(name = 'Compartment',
                         values = c(B = c.b,
                                    A = c.a),
                         na.value = 'black') +
      geom_edge_link(color = 'black', width = .1) +
      geom_node_point(aes(color = factor(comp)), size = .5) +
      facet_graph(~smp) +
      scale_edge_color_gradientn(colors = viridis::cividis(10, end = .9),
                                 breaks = c(5,10,20,40),
                                 trans = 'log',
                                 name = 'Component size') +
      theme(legend.position = 'none',
            plot.background = element_blank(),
            panel.background = element_rect(fill = NA, color = clrs[s], size = 1),
            strip.background = element_rect(fill = NA),
            strip.clip = 'off',
            strip.text = element_text(color = 'black', size = 13, face = 'bold'),
            strip.text.x = element_text(color = clrs[s]))
    
    if (s %in% c('ESC','CM')) {
      p <- p +
        facet_graph(proj ~ smp, switch = 'y')
    } 
    p
  }, ., names(.))

l <- tibble(x = 1:2, y = c('A', 'B') %>% factor(c('A','B'))) %>%
  ggplot(aes(x, y, color = y)) + 
  geom_point()+
  scale_color_manual(name = 'Compartment',
                     values = c(A = c.a,
                                B = c.b)) +
  theme(legend.position = 'bottom',
        legend.background = element_blank(),
        legend.direction = 'vertical',
        legend.text = element_text(size = 11),
        legend.key = element_blank())

{ wrap_plots(c(ps, list(cowplot::get_legend(l))), nrow = 2) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f1_h.pdf', ., height = 4, width = 8)

