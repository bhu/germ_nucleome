library(data.table)
library(tidyverse)
library(CAGEfightR)
library(pbapply)
library(InteractionSet)
library(pals)
library(ggpubr)
library(ggdist)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# register(MulticoreParam(workers=44))
# ds <- list.files('../data/netcage', pattern = '.bw', full.names = T) %>%
#   tibble(f = .) %>%
#   mutate(b = basename(f)) %>%
#   separate(b, c('Name', 'sign', NA), '\\.') %>%
#   mutate(sign = c(plus = 'BigWigPlus', minus = 'BigWigMinus')[sign]) %>%
#   pivot_wider(names_from = "sign", values_from = "f") %>%
#   mutate(nm = Name) %>%
#   column_to_rownames("nm") %>%
#   DataFrame()
# 
# bws <- list(plus = ds$BigWigPlus,
#             minus = ds$BigWigMinus) %>%
#   lapply(function(x) {
#     BigWigFileList(x) %>%
#       `names<-`(ds$Name)
#   })
# gn <- fread('../data/resources/mm10.chrom.sizes') %>%
#   {Seqinfo(seqnames = .$V1, seqlengths = .$V2, genome = "mm10")}
# 
# 
# ctss <- quantifyCTSSs(plusStrand = bws$plus,
#                       minusStrand = bws$minus,
#                       design = ds,
#                       genome = gn)
# 
# ctss <- calcTPM(ctss) %>%
#   calcPooled(inputAssay = "TPM") %>%
#   calcSupport(inputAssay = "counts",
#               outputColumn = "support",
#               unexpressed = 0)
# sctss <- subset(ctss, support > 0)
# sctss <- calcTPM(sctss)
# 
# tcs <- clusterUnidirectionally(sctss,
#                                pooledCutoff = 0.1,
#                                mergeDist = 20)
# 
# bcs <- clusterBidirectionally(ctss, balanceThreshold = 0.8) %>%
#   calcBidirectionality(samples = ctss)
# enhs <- subset(bcs, bidirectionality > 0)
# 
# tcs <- quantifyClusters(ctss,
#                         clusters = tcs,
#                         inputAssay="counts") %>%
#   calcTPM(totalTags = "totalTags") %>%
#   calcSupport(inputAssay = "TPM",
#               outputColumn = "support",
#               unexpressed = 1)
# 
# enhs <- quantifyClusters(ctss,
#                          clusters = enhs,
#                          inputAssay="counts") %>%
#   calcTPM(totalTags = "totalTags") %>%
#   calcSupport(inputAssay = "TPM",
#               outputColumn = "support",
#               unexpressed = 1)
# 
# tcs2 <- subset(tcs, support > 0)
# enhs2 <- subset(enhs, support > 0)
# 
# tcs2 <- assignTxType(tcs2, txModels = txdb)
# enhs2 <- assignTxType(enhs2, txModels = txdb)
# 
# enhs2 <- subset(enhs2, txType %in% c("intergenic", "intron"))
# 
# rowRanges(tcs2)$clusterType <- "TSS"
# rowRanges(enhs2)$clusterType <- "enhancer"
# 
# colData(enhs2) <- colData(tcs2)
# SE <- combineClusters(object1 = tcs2, 
#                       object2 = enhs2,
#                       removeIfOverlapping = "object1")
# 
# rowRanges(SE)$clusterType <- factor(rowRanges(SE)$clusterType,
#                                     levels = c("TSS", "enhancer"))
# 
# SE <- calcTPM(SE, totalTags="totalTags")
# 
# TCBC_pairs <- findLinks(SE, 
#                         inputAssay = "TPM",
#                         maxDist = 1000000, 
#                         directional = "clusterType",
#                         method = "kendall")
load('../data/netcage/cagefightr.rda')

# lthres <- .5
# ethres <- 1
# 
# ses <- rowRanges(SE) %>% split(., .$clusterType)
# o <- colnames(SE) %>%
#   setNames(., .) %>%
#   lapply(function(s) {
#     gr <- rowRanges(SE) %>% .[assay(SE, 'TPM')[,s] > ethres] %>% split(., .$clusterType)
#     l <- TCBC_pairs[TCBC_pairs$estimate > lthres]
#     grl <- sapply(gr, length)
#     obs <- (overlapsAny(anchors(l)$first, gr$TSS) & 
#               overlapsAny(anchors(l)$second, gr$enhancer)) |
#       (overlapsAny(anchors(l)$second, gr$TSS) & 
#          overlapsAny(anchors(l)$first, gr$enhancer))
#     
#     oo <- pblapply(1:100, function(ii) {
#       l <- TCBC_pairs[TCBC_pairs$estimate > lthres] %>%
#         sample(replace = T) 
#       obs <- (overlapsAny(anchors(l)$first, gr$TSS) & 
#                 overlapsAny(anchors(l)$second, gr$enhancer)) |
#         (overlapsAny(anchors(l)$second, gr$TSS) & 
#            overlapsAny(anchors(l)$first, gr$enhancer))
#       sum(obs)
#     })
#     
#     obs <- sum(obs)
#     rand <- pblapply(1:100000, function(i) {
#       gr <- list(TSS = sample(ses$TSS, grl['TSS']),
#                  enhancer = sample(ses$enhancer, grl['enhancer']))
#       k <- (overlapsAny(anchors(l)$first, gr$TSS) & 
#               overlapsAny(anchors(l)$second, gr$enhancer)) |
#         (overlapsAny(anchors(l)$second, gr$TSS) & 
#            overlapsAny(anchors(l)$first, gr$enhancer))
#       sum(k)
#     })
#     list(obs, rand)
#   })

load('../data/netcage/permutation.rda')
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}
samps <- c('ESC', 'EpiLC', 'd4c7PGCLC','GSC') 

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}


clrs <- dark <- setNames(tableau20(11)[seq(1,10,2)[-3]], samps)
light <- setNames(tableau20(11)[seq(2,10,2)[-3]], samps)
pd <- lapply(o, function(x) {
  unlist(x$rand) %>%
    tibble(v = .) %>%
    mutate(i = 1:n(),
           obs = x$obs)
}) %>% bind_rows(.id = 'samp') %>%
  mutate(type = sub('_.*', '', samp) %>%
           factor(rev(samps)),
         rep = sub('.*_', '', samp)) %>%
  na.omit() 

ggplot() +
  stat_slab(aes(x = type, y = v, group = samp, fill = type),  side = 'left',
               scale = .5, data = pd[pd$rep == '1',], position = position_nudge(x = -.1), alpha = .8) +
  stat_slab(aes(x = type, y = v, group = samp, fill = type), side = 'right',
               scale = .5, data = pd[pd$rep == '2',], position = position_nudge(x = +.1), alpha = .8) +
  coord_flip(xlim = c(1,4.1)) +
  scale_color_manual(values = dark) +
  scale_fill_manual(values = light) +
  geom_point(aes(x = type, y = obs, color = type, group = samp), shape = 18, size = 3,
             data = distinct(pd, samp, type, obs), position = position_dodge(.4)) +
  labs(y = '# of co-transcribed E-P pairs') +
  facet_grid(.~'NET-CAGE co-expression') +
  scale_x_discrete(labels= rnm) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.y = element_blank()) -> p
ggsave('sf5_f_left.pdf',p, height = 4.9, width = 3.2)


pd2 <- lapply(o, function(x) {
  tibble(obs = x$obs,
         exp = mean(unlist(x$rand)),
         std = sd(unlist(x$rand)),
         p_hi = (sum(x$obs < unlist(x$rand)) + 1) / (length(x$rand)+1),
         p_lo = (sum(x$obs > unlist(x$rand)) + 1) / (length(x$rand)+1),
         n = length(x$rand))
}) %>% bind_rows(.id = 'samp') %>%
  mutate(type = factor(sub('_.*', '', samp), rev(samps))) %>%
  na.omit() %>%
  rowwise() %>%
  mutate(p = min(2 * p_hi, 2 * p_lo),
         z = (obs - exp) / std,
         y.position = (obs / exp) * 1.05,
         p.signif = as.character(signif.num(p))) %>%
  ungroup() %>%
  separate(samp, c('group1', 'group2'), '_', F)


ggplot(pd2, aes(x = type, y = obs/exp, group = group2)) +
  geom_hline(yintercept = 1, color = 'grey70') +
  geom_point(aes(color = type), position = position_dodge(.7)) +
  stat_pvalue_manual(pd2, label = 'p.signif', x = 'type',
                     position = position_dodge(.7), coord.flip = T) +
  scale_color_manual(values = clrs) +
  labs(y = 'O/E # of co-expressed E-P pairs') +
  facet_grid(.~'NET-CAGE associations') +
  scale_y_log10(limits = c(.8, 1.8), breaks = c(.8,1,1.25, 1.5)) +
  scale_x_discrete(labels = rnm) +
  coord_flip() +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        strip.clip = 'off',
        axis.title.y = element_blank()) -> p
ggsave('sf5_f_right.pdf',p, height = 4.9, width = 3.6)

