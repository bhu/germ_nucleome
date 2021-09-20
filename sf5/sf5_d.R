library(data.table)
library(tidyverse)
library(pals)
library(gghalves)
library(rtracklayer)

clrs <- c(blue = '#0040FF', red = '#EB2224', yellow = '#EEE000', 'chr3\ntrans.' = '#8AC43F')

d <- c('EpiLC','GSC') %>%
  lapply(function(x) {
    cons <- fread(sprintf("../data/chrY/%s_trimmed_bismark_bt2.CX_report_all_CG.bedGraph", x), 
                  col.names = c("chr", "start", "end", "strand", "M", "U")) %>%
      mutate(chr = sub('yelow', 'yellow', chr),
             up = (1:n()) %% 2 == 1, 
             end = ifelse(up, end + 1, end), 
             start = ifelse(up, start, start - 1)) %>%
      group_by(chr, start, end) %>%
      summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
      mutate(depth = M + U,
             perc = 100 * M / depth,
             reg = 'all',
             chr = sub('tr_chr3', 'chr3\ntrans.', chr) %>% factor(names(clrs)),
             type = 'chrY amplicons\n from Soh et al., 2014') %>%
      select(-end) %>%
      arrange(chr) 
    trad <- fread(sprintf("../data/chrY/%s_allCG_list.txt.gz", x), 
                         col.names = c("chr", "start", "M", "U", "ratio")) %>%
      filter(chr %in% c("chr1", "chrY")) %>%
      mutate(up = (1:n()) %% 2 == 1, 
             start = ifelse(up, start, start - 1)) %>%
      group_by(chr, start) %>%
      summarise(across(c(M, U), sum), .groups = 'drop') %>%
      mutate(depth = M + U,
             perc = 100 * M / depth,
             end = start + 1) %>%
      {mutate(., reg = overlapsAny(makeGRangesFromDataFrame(.), import.bed('../data/PMD/GSC.bed')))} %>%
      mutate(reg = ifelse(reg, 'PMD', 'non-PMD'),
             type = 'Reference assembly') %>%
      select(-end)  %>%
      arrange(chr)
    rbind(trad, cons) %>%
      mutate(samp = x, chr = fct_inorder(chr))
  }) %>%
  bind_rows() %>%
  filter(depth > 5) %>% 
  select(-depth, -M, -U) %>% 
  pivot_wider(names_from = 'samp', values_from = 'perc') %>% 
  na.omit() %>% 
  mutate(dif = GSC - EpiLC, 
         type = factor(type, c('Reference assembly', 'chrY amplicons\n from Soh et al., 2014')),
         r = case_when(type == 'Reference assembly' ~ sprintf('%s (%s)', chr, reg),
                       T ~ sub('yelow', 'yellow', chr)) %>% fct_inorder()) 

clrs <- c(clrs, setNames(c(stepped2()[c(17,19)], stepped3()[c(5,7)]), 
                         c('chr1 (non-PMD)', 'chr1 (PMD)', 'chrY (non-PMD)', 'chrY (PMD)')))


ggplot(d, aes(x = r, y = dif, color = r, fill = r)) + 
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_half_violin(color = NA, side = 'l', alpha = .5) +
  geom_boxplot(outlier.color = NA, width = .3, notch = T, position = position_nudge(x = .25)) + 
  stat_summary(geom = "crossbar", width = 0.15, fatten = 0, color = "white", position = position_nudge(x = .25),
               fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  facet_grid(.~type, scales = 'free_x') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_y_continuous(breaks = c(-50,0,50), expand = expansion(0)) +
  coord_cartesian(ylim = c(-100,100)) +
  labs(x = 'Region', y = expression(Delta~'mCG% (GSC - EpiLC)')) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.line.x = element_line(color = 'black'),
        legend.position = 'none',
        strip.clip = 'off',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        plot.margin = margin(5,5,5,20),
        panel.grid = element_blank()) +
  ggsave('sf5_d.pdf', height = 4.37, width = 5.8, device = cairo_pdf, bg = 'transparent')
