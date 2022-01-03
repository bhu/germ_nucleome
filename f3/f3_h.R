library(data.table)
library(tidyverse)
library(readxl)
library(pals)
library(patchwork)
library(scales)

samps <- c("d4c7PGCLC", "GSC")
clrs <- setNames(tableau20()[c(7, 9)], samps)
marks <- c("CTCF", "methyl", "K9me2", "K9me3", "K36me2", "K36me3", "insul")
load('../data/aggregate/d4c7.GSC.CTCF.rda')

read_excel('../data/histone_ratios.xlsx') %>%
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
  filter(grepl('^H3', pep) & mod != 'unmod') %>%
  mutate(mod = strsplit(mod, '(?<=.)(?=K)', perl = T)) %>% 
  unnest(mod) %>%
  dplyr::select(-pep) %>%
  group_by(mod) %>%
  summarise(across(everything(), sum), .groups = 'drop') %>%
  pivot_longer(-mod, names_to = 'samp', values_to = 'a') %>%
  separate(samp, c('type', 'rep'), '_', F) %>%
  filter(mod %in% marks & type %in% samps) %>%
  group_by(mod, type) %>%
  summarise(a = mean(a), .groups = 'drop_last') %>%
  mutate(a = n() * a / sum(a)) %>%
  ungroup() %>%
  dplyr::rename(mark = mod) %>%
  merge(separate(a, track, c('type', 'mark')), all.y = T) %>%
  filter(mark %in% marks) %>%
  mutate(a = replace_na(a, 1),
         se = sd / sqrt(num) * a,
         mu = mu * a,
         mark = factor(mark, marks)) %>%
  arrange(mark) %>%
  mutate(mark = as.character(mark) %>% 
           sub('^K', 'H3K', .) %>%
           sub('methyl', 'DNAme', .) %>%
           fct_inorder()) %>%
  filter(reg == 'd4c7_up') %>%
  split(., .$mark) %>%
  lapply(function(d) {
    p <- if (d$mark[1] != 'insul') {
      ggplot(d, aes(x = x, y = mu, ymin = mu - se, ymax = mu + se, color = type, fill = type)) +
        geom_vline(xintercept = 30.5, color = "grey", size = 0.5) +
        geom_ribbon(alpha = .5, color = NA) +
        geom_line() +
        scale_x_continuous(breaks = c(1, 30.5, 60),
                           labels = c("-3kb", "", "+3kb")) +
        scale_color_manual(values = clrs) +
        scale_fill_manual(values = clrs) +
        scale_y_continuous(breaks = pretty_breaks(n = 3)) +
        facet_grid(.~mark) +
        ylab('Normalized CPM') +
        theme(legend.position = "none",
              panel.background = element_blank(),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.line.y = element_blank(),
              axis.line.x = element_line(size = 0.5, color = "black"),
              axis.text = element_text(color = "black", size = 11),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_line(color = "grey", linetype = "dashed"),
              strip.background = element_rect(fill = NA),
              strip.clip = 'off',
              strip.text = element_text(color = "black", size = 13, face = 'bold'),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.key = element_blank(),
              legend.title = element_blank(),
              plot.margin = margin(0,0,0,5))
    } else {
      
      ggplot(d, aes(x = x, y = mu, ymin = mu - se, ymax = mu + se, color = type, fill = type)) +
        geom_vline(xintercept = 40.5, color = "grey", size = 0.5) +
        geom_ribbon(alpha = .5, color = NA) +
        geom_line() +
        scale_x_continuous(breaks = c(1, 40.5, 80),
                           labels = c("-1mb", "", "+1mb")) +
        scale_color_manual(values = clrs) +
        scale_fill_manual(values = clrs) +
        facet_grid('GSC-low' ~ 'Insulation') +
        ylab(expression(log[2]~'IS')) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.line.y = element_blank(),
              axis.line.x = element_line(size = 0.5, color = "black"),
              axis.text = element_text(color = "black", size = 11),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_line(color = "grey", linetype = "dashed"),
              strip.background = element_rect(fill = NA),
              strip.text = element_text(color = "black", size = 13, face = 'bold'),
              legend.box.background = element_blank(),
              legend.background = element_blank(),
              legend.key = element_blank(),
              strip.clip = 'off',
              legend.title = element_blank(),
              plot.margin = margin(0,0,0,5)) 
    }
    if (d$mark[1] != 'CTCF' & d$mark[1] != 'insul') {
      p <- p +
        theme(axis.text.x = element_blank(),
              axis.title.y = element_blank())
    }
    p
  }) %>%
  {wrap_plots(., nrow = 1) &
  theme(plot.background = element_blank())} %>%
  ggsave('f3_h.pdf', ., height = 1.5, width = 8)

a %>%
  separate(track, c('type', 'mark')) %>%
  ggplot(aes(x = x, y = mu, color = type)) + geom_line(size = 1) + 
  scale_color_manual(values = clrs, labels = function(x) sub('PGC', ' mPGC', x)) +
  theme(legend.position = 'bottom',
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank()) -> p
  ggsave('f3_h.leg.pdf', p, width = 8.2, height = 2)
