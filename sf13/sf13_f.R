library(data.table)
library(tidyverse)
library(pals)
library(rtracklayer)
library(patchwork)
library(scales)
library(ggh4x)


clrs <- c(blue = '#0040FF', red = '#EB2224', yellow = '#EEE000', 'chr3\ntrans.' = '#8AC43F')

cons <- fread("../data/chrY/GSC_trimmed_bismark_bt2.CX_report_all_CG.bedGraph", 
      col.names = c("chr", "start", "end", "strand", "M", "U")) %>%
  mutate(chr = sub('yelow', 'yellow', chr),
         up = (1:n()) %% 2 == 1, 
         end = ifelse(up, end + 1, end), 
         start = ifelse(up, start, start - 1)) %>%
  group_by(chr, start, end) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
  mutate(depth = M + U,
         perc = 100 * M / depth,
         chr = sub('tr_chr3', 'chr3\ntrans.', chr) %>% factor(names(clrs)))

trad <- fread("../data/chrY/GSC_allCG_list.txt.gz", 
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
  mutate(reg = ifelse(reg, 'PMD', 'non-PMD'))

deep <- trad %>%
  group_by(chr, reg) %>%
  summarise(r = 100 * sum(depth >= 5) / n())

p1 <- trad %>% 
  count(chr, reg, depth) %>%
  ggplot(aes(x = depth, y = n)) +
  geom_vline(xintercept = 5, color = 'firebrick', alpha = .5) +
  geom_col(width = .02, alpha = .5, position = 'identity') +
  geom_label(aes(label = sprintf("%s%% \u22655x", round(r,2)), x = Inf, y = Inf), 
             hjust = 1.05, vjust = 1.2, data = deep, color = 'firebrick', alpha = .5) + 
  scale_x_continuous('Coverage', trans = 'log1p',
                     breaks = c(0, 10, 100)) +
  coord_cartesian(xlim = c(0,100)) +
  ylab('# of CpGs') +
  facet_nested(chr + reg ~ 'Reference assembly', scales = "free_y", nest_line = T) +
  scale_y_continuous(breaks = pretty_breaks(2), labels = comma, expand = expansion(c(0, .05))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'))

deep2 <- cons %>%
  group_by(chr) %>%
  summarise(r = 100 * sum(depth >= 5) / n())

p2 <- cons %>% 
  count(chr, depth) %>%
  ggplot(aes(x = depth, y = n)) +
  geom_vline(xintercept = 5, color = 'firebrick', alpha = .5) +
  geom_col(width = .02, alpha = .5, position = 'identity') +
  geom_label(aes(label = sprintf("%s%% \u22655x", round(r, 2)), x = Inf, y = Inf), 
             hjust = 1.05, vjust = 1.2, data = deep2, color = 'firebrick', alpha = .5) + 
  scale_x_continuous('Depth', trans = 'log1p',
                     breaks = c(0, 10, 100)) +
  coord_cartesian(xlim = c(0,100)) +
  scale_y_continuous('# of CpGs', expand = expansion(c(0, .05)),
                     breaks = pretty_breaks(2)) +
  facet_grid(chr ~ 'chrY amplicons\n from Soh et al., 2014', scales = "free_y") +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        strip.clip = 'off',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.title = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'))


g2 <- ggplot_gtable(ggplot_build(p2))
strip <- which(grepl('strip-r', g2$layout$name))
k <- 1
for (i in strip) {
  j <- which(grepl('title', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- clrs[k]
  k <- k+1
}


{ wrap_plots(p1, wrap_ggplot_grob(g2), nrow = 1) &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf13_f.pdf', ., height = 4, width = 7.3, device = cairo_pdf, bg = 'transparent')

