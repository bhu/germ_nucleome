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

# perf <- d50 %>%
#   split(., .$samp) %>%
#   lapply(function(x) {
#     tmp <- x %>%
#       dplyr::select(-c(chr, start, end, samp, PC1,ATAC)) %>%
#       mutate(targ = factor(ifelse(pmd, 'PMD', 'nonPMD')))
#     set.seed(42)
#     idx <- createDataPartition(tmp$targ, p = .7, list = F)
#     d.train <- tmp[idx,]
#     d.test <- tmp[-idx,]
#     
#     c('glmnet', 'gbm', 'nnet') %>%
#       lapply(function(m) {
#         mod <- train(targ ~ ., data = d.train, method = m,
#                      trControl = trainControl(method = 'repeatedcv',
#                                               number = 10, repeats = 10,
#                                               verboseIter = T, sampling = 'smote'))
#         roc(d.test$targ, predict(mod, d.test, type = 'prob')[,'PMD']) %>%
#           ci.auc() %>%
#           as.numeric() %>%
#           {tibble(lower = .[1], auc = .[2], upper = .[3],
#                   method = m)} 
#       }) %>%
#       bind_rows()
#   }) %>%
#   bind_rows(.id = 'samp')

load('../data/PMD/auc.rda')

perf %>%
  mutate(samp = factor(samp, samps)) %>%
  na.omit() %>%
  ggplot(aes(x = samp, y = auc, ymin = lower, ymax = upper, color= samp)) +
  geom_line(aes(group = method), color = 'grey50', alpha = .5) +
  geom_point() +
  geom_linerange() +
  facet_grid(.~'Predicting GSC PMDs') +
  ylab('Model AUROC') +
  scale_color_manual(values = clrs) +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  scale_x_discrete(labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  facet_grid(method~'Predicting GSC PMDs') +
  theme(axis.text.x = element_text(angle = 45, hjust= 1),
        legend.position = 'none',
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "#ffffff66", color = NA),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.line.x = element_line(color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank(),
        strip.clip = 'off',
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y = element_blank())  -> p
ggsave('f5_l.pdf', p, height = 4.1, width = 2.9, device = cairo_pdf, bg = 'transparent')



