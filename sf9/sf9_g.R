library(tidyverse)
library(readxl)
library(ggh4x)
library(pals)
library(scales)

bsz <- 5e5
dec <- 100
load('../data/tracks/inpnorm.50kb.rda')

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
marks <- c('Laminb1', 'K9me3', 'K9me2')

rnm <- function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))
rnm2 <- function(x) sub('Laminb1', 'LMNB1', x)
  
sclrs <- setNames(tableau20(9)[seq(1, 9, 2)], rnm(samps))
mclrs <- setNames(c('forestgreen', 'darkorange3', 'darkorchid3'), rnm2(marks)) %>%
  c(bg = 'grey70')

cap <- d50 %>% 
  filter(chr == 'chr3' & samp %in% samps) %>%
  dplyr::select(Laminb1, K9me3, K9me2) %>% 
  sapply(function(x) unname(quantile(x, .999))) 

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

limits <- d %>%
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
  mutate(a = -log10(a / mean(a))) %>%
  ungroup() %>%
  rbind(tibble(mod = 'Laminb1', type = samps, a = 0)) %>%
  mutate(samp = factor(rnm(type), rnm(samps)),
         v = cap[mod] + a,
         v = case_when(mod == 'K9me2' ~ v * .7, 
                       mod == 'K9me3' ~ v * .7,
                       T ~ v),
         l = floor(v * dec) / dec,
         mark = factor(rnm2(mod), rnm2(marks))) %>%
  arrange(samp, mark) %>%
  rowwise() %>%
  mutate(l = list(scale_y_continuous(limits = c(0, v), breaks = c( 0, l), labels = c(0,l),
                                     expand = expansion(c(0, .1))))) %>%
  pull(l)

p <- d50 %>% 
  filter(chr == 'chr3' & samp %in% samps) %>%
  dplyr::select(all_of(marks), samp, x = start) %>% 
  mutate(x = round(x / bsz) * bsz + bsz / 2) %>%
  group_by(samp, x) %>% 
  summarise_all(mean) %>% 
  ungroup() %>%
  pivot_longer(-c(samp, x), names_to = 'mark', values_to = 'y') %>%
  mutate(mark = factor(rnm2(mark), rnm2(marks)),
         samp = factor(rnm(samp), rnm(samps)),
         pos = y > 0,
         up = ifelse(pos, y, 0),
         dn = ifelse(pos, 0, y),
         reg = ifelse(pos, as.character(mark), 'bg')) %>%
  ggplot(aes(x = x, y = y)) +
  #geom_hline(yintercept = 0, color = 'grey70') +
  #geom_line(aes(alpha = pos, group = 1)) +
  scale_alpha_manual(values = c(.1, .5)) +
  geom_ribbon(aes(ymin = 0, ymax = up, fill = mark), alpha = .5) +
  geom_ribbon(aes(ymin = 0, ymax = dn), fill = 'grey70', alpha = .25) +
  #geom_rug(aes(color = reg), sides = "b") +
  geom_vline(xintercept = Inf, color = 'black', size = 1) +
  scale_color_manual(values = mclrs) +
  scale_fill_manual(values = mclrs) +
  facet_nested(samp + mark ~ ., scales = 'free_y', nest_line = T) +
  ylab('log2 enrichment over input') +
  facetted_pos_scales(y = limits) +
  scale_x_continuous('chr3', breaks = c(0, 5e7, 1e8, 1.5e8),
                     limits = c(0, 1.6e8),
                     labels = c('0', '50mb', '100mb', '150mb'),
                     expand = expansion(0)) +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'))

g <- ggplot_gtable(ggplot_build(p))
strip <- which(grepl('strip-r', g$layout$name))[1:5]
fills <- sclrs
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k + 1
}

ggsave('sf7_g.pdf', g, height = 11.75, width = 6, bg = 'transparent')

