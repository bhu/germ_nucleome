library(data.table)
library(tidyverse)
library(pals)
library(ggh4x)
library(scales)
library(readxl)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC') 
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)
marks <- c('K36me2', 'K9me2', 'K9me3', 'Laminb1', 'DNAme')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
anchors <- c('ESC', 'EpiLC','d4c7PGCLC', 'GSC', 'GSCLC')
types <- c('ESC', 'EpiLC', 'd2PGCLC','d4c7PGCLC', 'GSC', 'MEF') 

d <- c('old', 'new') %>%
  setNames(., .) %>%
  lapply(function(s) {
    read_excel('../data/histone_ratios.xlsx', sheet = s) %>%
      {.[,colSums(is.na(.)) != nrow(.)]} %>% 
      na.omit() %>% 
      {.[,c(1, which(.[1,] == 'Area'))]} %>%
      tail(-1) %>%
      rename(p = 1) %>%
      `names<-`(., sub('\\..*', '', names(.))) %>%
      separate(p, c('pep', 'mod'), '\\ ') %>%
      mutate(pep = sub('H33', 'H3', pep)) %>%
      group_by(pep) %>%
      mutate(across(-mod, function(x) as.numeric(x) / sum(as.numeric(x)))) %>%
      ungroup() %>% 
      filter(grepl('^H[34]', pep) & mod != 'unmod') %>%
      mutate(mod = strsplit(mod, '(?<=.)(?=K)', perl = T)) %>% 
      unnest(mod) %>%
      mutate(mod = paste(substr(pep, 1, 2), mod)) %>%
      dplyr::select(-pep) %>%
      group_by(mod) %>%
      summarise(across(everything(), sum), .groups = 'drop') %>%
      pivot_longer(-mod, names_to = 'samp', values_to = 'a') %>%
      separate(samp, c('type', 'rep'), '_', F) 
  })

nrm <- d %>%
  lapply(function(x) {
    x %>%
      filter(type %in% anchors) %>%
      group_by(mod, type) %>%
      summarise(a = median(a), .groups = 'drop')
  }) %>%
  bind_rows(.id = 'batch') %>%
  pivot_wider(names_from = 'batch', values_from = 'a') %>%
  split(., .$mod) %>% 
  lapply(function(x) {
    fit <- tryCatch(
      rlm(old ~ new, x)$coefficients,
      error = function(msg) {
        lm(old ~ new, x)$coefficients
      },
      warning = function(msg) {
        lm(old ~ new, x)$coefficients
      }
    )
    d$new %>%
      filter(mod == x$mod[1]) %>%
      mutate(a = fit[1] + a * fit[2])
  }) %>%
  bind_rows() %>%
  mutate(batch = 'new') %>%
  rbind(mutate(d$old, batch = 'old')) %>%
  separate(mod, c('his', 'mod'), '\\ ') %>%
  filter(his == 'H3' & mod %in% marks & type %in% samps) %>%
  group_by(mod, type) %>%
  summarise(a = mean(a), .groups = 'drop_last') %>%
  mutate(a = a / mean(a)) %>%
  ungroup()

load('../data/aggregate/GSC.PMD.rda')

a %>%
  separate(track, c('samp', 'mark'), '_') %>%
  merge(nrm, all.x = T, by.x = c('samp', 'mark'), by.y = c('type', 'mod')) %>%
  mutate(a = replace_na(a, 1),
         mu = mu * a,
         sd = sd * a,
         se = sd / sqrt(num),
         mark = factor(mark, marks),
         samp = factor(samp, samps)) %>%
  filter(reg == 'PMD' & mark %in% marks) %>%
  arrange(samp, mark) %>%
  mutate(mark = as.character(mark) %>%
           sub('^K', 'H3K', .) %>%
           sub('b1', ' B1', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x, y = mu, ymin = mu - se, ymax = mu + se)) +
  annotate('rect', xmin = 20, xmax = 40, ymin = -Inf, ymax = Inf, fill = 'tan', alpha = .2) +
  geom_ribbon(aes(fill = samp), alpha = .3, color = NA) +
  geom_line(aes(color = samp)) +
  facet_wrap(.~ mark, nrow = 1, scales = "free") +
  ylab("Enrichment") +
  scale_color_manual(values = clrs, labels = rnm) +
  scale_fill_manual(values = clrs, labels = rnm) +
  facetted_pos_scales(x = list(
    scale_x_continuous('PMD',
                       breaks = c(1, 20, 40, 60),
                       labels = c("-1mb", "", "", "+1mb")),
    scale_x_continuous(breaks = c(1, 20, 40, 60),
                        labels = c("","","","")),
    scale_x_continuous(breaks = c(1, 20, 40, 60),
                       labels = c("","","","")),
    scale_x_continuous(breaks = c(1, 20, 40, 60),
                       labels = c("","","","")),
    scale_x_continuous(breaks = c(1, 20, 40, 60),
                       labels = c("","","",""))
  )) +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.text = element_text(color = "black", size = 11),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color = "grey70", linetype = "dashed"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = "black", size = 13, face = 'bold'),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11),
        legend.margin = margin(t = -10),
        legend.title = element_blank()) -> p
ggsave('sf12_g.pdf', height = 3.87, width = 7.5)

