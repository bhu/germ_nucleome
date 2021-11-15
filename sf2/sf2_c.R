library(data.table)
library(tidyverse)
library(tidygraph)
library(igraph)
library(GenomicRanges)
library(boot)
library(pals)

load('../data/compscore/25kb.tads.rda')
load('../data/cliques/sig.ints.rda')

samps <- c('ESC (2i)','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
           'ESC (serum)', 'NPC', 'CN',
           'Day 6','Day 4','Day 2','B\u03B1',
           'CM', 'CPC', 'PCM', 'VCM')
rnm <- function(x) {
  #x %>%
  #  sub('C ', 'Cs ', .) %>%
  #  sub('^(.*[CMN])$', '\\1s', .) 
  x <- sub('ESC', 'mESC', x)
  x <- case_when(x == 'ESC' ~ 'mESC',
                 x == 'd2PGCLC' ~ 'd2 mPGCLC',
                 x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
                 T ~ x)
}
samps <- rnm(samps)

clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(12)[-4]), samps)

# res <- list(germ = c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC'),
#             neural = c('ESC', 'NPC', 'CN'),
#             stadhouders = c('D6', 'D4', 'D2', 'Ba'),
#             cardio = c('CM', 'CPC', 'PCM', 'VCM')) %>%
#   lapply(function(x) setNames(x, x)) %>%
#   Map(function(ss, p) {
#     tad <- tads[[p]] %>%
#       `colnames<-`(c('chr','start','end')) %>%
#       mutate(start = start + 1) %>%
#       makeGRangesFromDataFrame()
#     nds <- as_tibble(tad)[,1:3] %>%
#       mutate(idx = sprintf('%s:%d-%d', seqnames, start-1, end)) %>%
#       pull(idx)
#     ss %>%
#       lapply(function(s) {
#         print(s)
#         gr <- na.omit(eig[[p]][[s]]) %>%
#           mutate(comp = case_when(E1 > 0 ~ 'A', T ~ 'B'),
#                  start = start + 1) %>%
#           dplyr::select(chrom, start, end, E1, comp) %>%
#           makeGRangesFromDataFrame(keep.extra.columns = T)
#         ab <- findOverlaps(tad, gr) %>%
#           as("List") %>%
#           extractList(gr$comp, .) %>%
#           sapply(function(y) {
#             ifelse(sum(y == 'A') > sum(y == 'B'), 'A', 'B')
#           }) %>%
#           setNames(nds)
#         ints <-  o[[p]][[s]] %>%
#           mutate(i = sprintf('%s:%d-%d', V1, V2, V3),
#                  j = sprintf('%s:%d-%d', V4, V5, V6)) %>%
#           select(i, j)
#         cs <- ints %>%
#           as.matrix() %>%
#           graph_from_edgelist() %>%
#           max_cliques() %>%
#           lapply(function(y) {
#             tibble(node = names(y)) %>%
#               add_count()
#           }) %>%
#           bind_rows(.id = 'idx') %>%
#           filter(n > 2) %>% 
#           separate(node, c('V1', 'V2', 'V3'), '[:-]', F) %>%
#           arrange(idx, V1, V2) %>%
#           select(idx, node) %>%
#           setDT() %>%
#           .[, {tmp <- combn(node, 2); .(i = tmp[1,], j = tmp[2,])}, idx] %>%
#           select(-idx) %>%
#           merge(ints)
#         freq <- c(cs$i, cs$j) %>%
#           unique() %>%
#           tibble(nd = .) %>%
#           mutate(cmp = ab[nd]) %>%
#           count(cmp) %>%
#           mutate(n = n / sum(n)) %>%
#           deframe()
#         obs <- cs %>% 
#           mutate(i = ab[i], j = ab[j]) %>% 
#           as.data.frame()
#         lapply(c('A-A', 'B-B', 'Het'), function(x) {
#           fx <- if (x == 'A-A') {
#             function(d, i) {
#               d2 <- d[i, ]
#               return(sum(d2$i == 'A' & d2$j == 'A'))
#             }
#           } else if (x == 'B-B') {
#             function(d, i) {
#               d2 <- d[i, ]
#               return(sum(d2$i == 'B' & d2$j == 'B'))
#             }
#           } else {
#             function(d, i) {
#               d2 <- d[i, ]
#               return(sum(d2$i != d2$j))
#             }
#           }
#           set.seed(42)
#           b <- boot(obs, fx, 50000)
#           boot.ci(b, type = 'perc') %>%
#             .$percentile %>%
#             as.data.frame() %>%
#             `colnames<-`(c('level', 'i1', 'i2', 'v1', 'v2')) %>% 
#             mutate(exp = case_when(x == 'A-A' ~ freq['A']^2,
#                                    x == 'B-B' ~ freq['B']^2,
#                                    T ~ 2 * freq['A'] * freq['B']) * nrow(obs),
#                    k = x,
#                    obs = b$t0)
#         }) %>% bind_rows()
#       }) %>% bind_rows(.id = 'samp') -> rr
#   }, ., names(.)) %>% bind_rows(.id = 'proj')

load('../data/cliques/clique.ab.bootstrap.rda')

res %>%
  mutate(x = case_when(samp == 'Ba' ~ 'B\u03B1',
                       samp == 'ESC' & study == 'germ' ~ 'ESC (2i)',
                       samp == 'ESC' & study == 'neural' ~ 'ESC (serum)',
                       T ~ samp) %>%
           sub('^D', 'Day ', .) %>%
           rnm() %>%
           factor(samps),
         study = c(germ = 'Germline', neural = 'Neural',
                   stadhouders = 'Reprogram',
                   cardio = 'Cardiac')[study] %>%
           factor(c('Germline', 'Neural', 'Reprogram', 'Cardiac')),
         k = ifelse(k == 'Het', 'A-B & B-A', k)) %>%
  #filter(k == 'A-B & B-A') %>%
  na.omit() %>%
  ggplot(aes(x = x, y = obs / exp, color = x)) +
  geom_hline(yintercept = 1, color = 'grey70') +
  geom_line(aes(group = study), color = 'grey50', alpha = .5, size = 1) +
  geom_linerange(aes(ymin = v1 / exp, ymax = v2 / exp)) +
  geom_point(size = 3) +
  facet_grid(k~study, scales = 'free', space = 'free_x') +
  scale_color_manual(values = clrs) +
  ylab('O/E # of clique interactions') +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())  -> p
  ggsave('sf2_c.pdf', height = 4.84, width = 5.8, device = cairo_pdf, bg = 'transparent')
