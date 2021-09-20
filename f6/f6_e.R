library(tidyverse)
library(pals)
library(rtracklayer)
library(DMwR)
library(pals)
library(caret)
library(pROC)
library(scales)

load('../data/tracks/inpnorm.50kb.rda')

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

pmd <- d50 %>% mutate(start = start + 1) %>% 
  distinct(chr, start, end) %>%
  makeGRangesFromDataFrame() %>% 
  overlapsAny(import.bed('../data/PMD/GSC.bed'))

# o <- d50 %>%
#   split(., .$samp) %>%
#   lapply(function(x) {
#     tmp <- x %>%
#       dplyr::select(-c(chr, start, end, samp, PC1,ATAC)) %>%
#       mutate(targ = factor(ifelse(pmd, 'PMD', 'nonPMD')))
#     set.seed(42)
#     idx <- createDataPartition(tmp$targ, p = .7, list = F)
#     d.train <- tmp[idx,]
#     d.test <- tmp[-idx,]
#     mod <- train(targ ~ ., data = d.train, method = 'gbm',
#                  trControl = trainControl(method = 'repeatedcv',
#                                           number = 10, repeats = 10,
#                                           verboseIter = T, sampling = 'smote'))
#     list(mod = mod, dat = tmp, d.train = d.train, d.test = d.test, idx = idx)
#   })
# 
# perf <- lapply(o, function(x) {
#   roc(x$d.test$targ, predict(x$mod, x$d.test, type = 'prob')[,'PMD']) %>%
#     ci.auc() %>%
#     as.numeric() %>%
#     {tibble(lower = .[1], auc = .[2], upper = .[3])}
# })

load('../data/PMD/auc.rda')

perf %>%
  bind_rows(.id = 'samp') %>%
  mutate(samp = factor(samp, samps)) %>%
  na.omit() %>%
  ggplot(aes(x = samp, y = auc, ymin = lower, ymax = upper, color= samp)) +
  geom_point() +
  geom_linerange() +
  facet_grid(.~'Predicting GSC PMDs') +
  ylab('Model AUROC') +
  scale_color_manual(values = clrs) +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  theme(axis.text.x = element_text(angle = 30, hjust= 1),
        legend.position = 'none',
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#ffffff66", color = NA),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank(),
        strip.clip = 'off',
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank()) +
  ggsave('f6_e.pdf', height = 2, width = 2.5)


