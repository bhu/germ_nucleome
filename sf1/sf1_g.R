library(data.table)
library(tidyverse)
library(GenomicRanges)
library(ggpubr)

load('../data/tads/insulscore.rda')

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}

# dat <- list(germ = c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC'),
#             neural = c('ESC', 'NPC', 'CN'),
#             stadhouders = c('D6', 'D4', 'D2', 'Ba'),
#             cardio = c('CM', 'CPC', 'PCM', 'VCM')) %>%
#   lapply(function(x) setNames(x, x)) %>%
#   Map(function(ss, p) {
#     grs <- tads[[p]][ss] %>%
#       lapply(function(x) {
#         x %>%
#           filter(boundary_strength_100000 > .2 & !is_bad_bin) %>%
#           mutate(start = start - 1e4 + 1,
#                  end = start + 1e4) %>%
#           makeGRangesFromDataFrame()
#       })
#     u <- Reduce(c, grs) %>% reduce()
#     obs <- lapply(grs, function(x) overlapsAny(u, x)) %>%
#       bind_cols() %>%
#       rowSums() %>%
#       tibble(num = .) %>%
#       mutate(grp = case_when(
#         num == 1 ~ '1',
#         num == 2 ~ '2',
#         num > 2 ~ '3+'
#       )) %>%
#       count(grp) %>%
#       deframe()
#     inp <- tads[[p]][ss] %>%
#       bind_rows(.id = 'samp') %>%
#       filter(boundary_strength_100000 > .2 & !is_bad_bin) %>%
#       mutate(start = start - 1e4 + 1,
#              end = start + 1e4)
#     lapply(1:100000, function(x) {
#       grs <- inp %>%
#         mutate(samp = sample(samp)) %>%
#         split(., .$samp) %>%
#         lapply(makeGRangesFromDataFrame) %>%
#         lapply(function(y) overlapsAny(u, y)) %>%
#         bind_cols() %>%
#         rowSums() %>%
#         tibble(num = .) %>%
#         mutate(grp = case_when(
#           num == 1 ~ '1',
#           num == 2 ~ '2',
#           num > 2 ~ '3+'
#         )) %>%
#         count(grp)
#     }) %>%
#       bind_rows(.id = 'idx') %>%
#       split(., .$grp) %>%
#       Map(function(x, grp) {
#         tibble(n = obs[grp], p = (sum(x$n > obs[grp])+1) / (nrow(x)+1))
#       }, ., names(.)) %>%
#       bind_rows(.id = 'grp') %>%
#       mutate(perc = 100 * n / sum(n))
#   }, ., names(.)) %>%
#   bind_rows(.id = 'study') %>%
#   mutate(study = c(germ = 'Germline',
#                    neural = 'Neural',
#                    stadhouders = 'Reprogramming',
#                    cardio = 'Cardiac')[study] %>%
#            fct_inorder(),
#          p.signif = as.character(signif.num(p)),
#          ypos = n + 250,
#          group1 = grp,
#          group2 = grp)

load('../data/tads/overlap.permutation.rda')
m <- 1e2
b <- 0

sclr <- 'firebrick3'
dat <- dat %>%
  mutate(study = as.character(study) %>%
           sub('Germ', 'Germline', .) %>%
           sub('ming', '', .) %>%
           fct_inorder()) 
ggplot(dat, aes(x = grp, y = n)) +
  geom_col(aes(fill = grp), alpha = .7) +
  geom_point(aes(y = m * perc + b), color = sclr) +
  stat_pvalue_manual(dat, x = 'grp', y.position = 'ypos', label = 'p.signif') +
  facet_wrap(~ study, nrow = 2) +
  scale_fill_manual(values = blues9[c(4, 7, 9)]) +
  scale_y_continuous('# of boundaries', expand = expansion(c(0, .05)), 
                     sec.axis = sec_axis(~ (. - b) / m, '% of all boundaries')) +
  xlab('# of samples boundary found in') +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.y.right = element_text(color = sclr),
        axis.title.y.right = element_text(color = sclr, vjust = 1.5),
        panel.grid = element_blank()) -> p
  ggsave('sf1_g.pdf', p, height = 4.6, width = 3.1)

