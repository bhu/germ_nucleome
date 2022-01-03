library(tidyverse)
library(readxl)
library(pals)
library(MASS)
library(patchwork)

anchors <- c('ESC', 'EpiLC','d4c7PGCLC', 'GSC', 'GSCLC')
types <- c('ESC', 'EpiLC', 'd2PGCLC','d4c7PGCLC', 'GSC', 'MEF') 
clrs <- setNames(tableau20(20)[c(1, 3, 5, 7, 9, 17)], types)
mods <- c("K4me1", "K4me3", "K9me2", "K9me3", "K27me3", "K36me2", "K36me3", "K9ac","K14ac", "K18ac", "K23ac", "K27ac")
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}

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

pd <- lapply(d, function(x) {
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
  filter(his == 'H3') %>%
  mutate(type = factor(type, types),
         mod = factor(mod, mods)) %>%
  na.omit() %>%
  arrange(mod, type) %>%
  mutate(mod = paste0('H3', mod) %>% fct_inorder())

p1 <- ggplot(pd, aes(x = type, y = a * 100)) +
  geom_hline(yintercept = c(-Inf,Inf), color = 'black', size = 1) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), 
               geom = 'pointrange', mapping = aes(color = type)) +
  facet_grid(.~mod, scales = 'free_x') +
  scale_color_manual(values = clrs) +
  coord_cartesian(ylim = c(0,60)) +
  scale_y_continuous(expand = expansion(c(0,.05))) +
  scale_x_discrete(labels = rnm) +
  ylab('Abundance (%)') +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.clip = "off",
        panel.spacing = unit(1.5, "lines"),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        panel.grid.major.x = element_line(color = 'grey90')) 

p2 <- pd %>%
  filter(mod == 'H3K4me3') %>%
  ggplot(aes(x = type, y = a * 100)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
               geom = 'pointrange', mapping = aes(color = type)) +
  scale_color_manual(values = clrs) +
  coord_cartesian(ylim = c(.1,.3)) +
  scale_y_continuous(expand = expansion(c(.15,.08)), breaks = c(.1,.2,.3)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'white', size = 13, face = 'bold'),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.title = element_blank(),
        #panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        #panel.grid.major.x = element_line(color = 'grey90'),
        panel.grid = element_blank()) 

p3 <- pd %>%
  filter(mod == 'H3K27ac') %>%
  ggplot(aes(x = type, y = a * 100)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
               geom = 'pointrange', mapping = aes(color = type)) +
  scale_color_manual(values = clrs) +
  coord_cartesian(ylim = c(.1,.5)) +
  scale_y_continuous(expand = expansion(c(.15,.08)), breaks = c(.1,.3,.5)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.title = element_blank(),
        #panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        #panel.grid.major.x = element_line(color = 'grey90'),
        panel.grid = element_blank()) 

p4 <- pd %>%
  filter(mod == 'H3K9ac') %>%
  ggplot(aes(x = type, y = a * 100)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
               geom = 'pointrange', mapping = aes(color = type)) +
  scale_color_manual(values = clrs) +
  coord_cartesian(ylim = c(.2,1.4)) +
  scale_y_continuous(expand = expansion(c(.15,.08)), breaks = c(.2,.8,1.4)) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.title = element_blank(),
        #panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        #panel.grid.major.x = element_line(color = 'grey90'),
        panel.grid = element_blank())

p1 + 
  inset_element(p2, .057,.135,.153,1.06) + 
  inset_element(p3, .9101, .135, 1.006, 1.06) +
  inset_element(p4, .5689, .135, .665, 1.06) &
  theme(plot.background = element_blank()) &
  ggsave('f2_a.pdf', height = 3, width = 13.3)

