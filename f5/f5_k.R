library(data.table)
library(tidyverse)
library(readxl)
library(pals)
library(cowplot)
library(patchwork)
library(scales)

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC') 
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)
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
  filter(his == 'H3' & mod == 'K9me3' & type %in% samps) %>%
  group_by(mod, type) %>%
  summarise(a = mean(a), .groups = 'drop_last') %>%
  mutate(a = a / mean(a)) %>%
  ungroup()

load('../data/aggregate/TSS.K9.rda')
a %>%
  separate(track, c('samp', 'mark'), '_') %>%
  merge(nrm, all.x = T, by.x = c('samp', 'mark'), by.y = c('type', 'mod')) %>%
  mutate(a = replace_na(a, 1),
         mu = mu * a,
         sd = sd * a,
         se = sd / sqrt(num),
         mark = sub('K', 'H3K', mark),
         samp = factor(samp, samps)) %>%
  filter(reg != 'other') %>%
  ggplot(aes(x = x, y = mu, ymin = mu - se, ymax = mu + se, color = samp, fill = samp)) +
  geom_vline(xintercept = 30.5, color = 'grey', size = .5) +
  geom_ribbon(alpha = .3, color = NA) +
  geom_line() + 
  facet_grid(mark ~ 'Setdb1-repressed germline genes', scales = 'free_y') +
  scale_color_manual(values = clrs, guide = guide_legend(nrow = 1), labels = rnm) +
  scale_fill_manual(values = clrs, labels = rnm) +
  scale_y_continuous('Normalized coverage', breaks = pretty_breaks(3), expand = expansion(.1)) +
  scale_x_continuous(breaks = c(11, 30.5, 50), expand = expansion(0),
                     labels = c("-2kb", "TSS", "+2kb")) +
  theme(panel.background = element_rect(color = 'black', size = 1, fill = NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text = element_text(color = "black", size = 11),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = "black", size = 13, face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'bottom',
        strip.clip = 'off',
        legend.title = element_blank(),
        legend.text = element_text(size = 11))  +
  ggsave('f5_k.pdf', height = 4.07,width=3.5)
  