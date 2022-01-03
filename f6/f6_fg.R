library(tidyverse)
library(pals)
library(GenomicRanges)
library(broom)
library(ggpubr)
library(patchwork)
library(boot)

load('../data/tracks/inpnorm.50kb.rda')
load('../data/calder.rda')

samps <- c("GSC",  'GSCLC') 
clrs <- setNames(c( '#9467bd', '#499894'), samps)
mrks <- c('PC1','ATAC','CTCF','Rad21',
          'Ring1b','H2Aub','K27me3', 'K27ac', 'K4me1','K4me3',
          'K36me2','K36me3','K9me2','K9me3', 'Laminb1')

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}


sc <- lapply(o$`50000`, makeGRangesFromDataFrame, keep.extra.columns = T)

bins <- distinct(d50, chr, start, end) %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame()

lbs <- o$`50000` %>%
  lapply(function(x) {
    makeGRangesFromDataFrame(x) %>%
      findOverlaps(bins, .) %>%
      as("List") %>%
      extractList(x$lab, .) %>%
      sapply(`[`, 1)
  }) %>%
  bind_cols()

d <- unique(o[[1]][[1]]$lab) %>%
  setNames(.,.) %>%
  lapply(function(x) {
    apply(lbs, 2, function(y) c(sum(y == x, na.rm = T), sum(!is.na(y)))) %>%
      t() %>%
      .[samps,] %>%
      prop.test() %>%
      tidy()
  }) %>%
  bind_rows(.id = 'sc') %>%
  mutate(q = p.adjust(p.value, method = 'fdr'),
         fc = log2(estimate2 / estimate1),
         symb = as.character(signif.num(q)),
         x = as.numeric(as.factor(sc))) %>%
  rowwise() %>%
  mutate(y = max(estimate1, estimate2) * 100 + 1.5) %>%
  ungroup() 

b <- 40
a <- 8

c2 <- 'firebrick2'

p1 <- d %>%
  select(sc, GSC = estimate1, GSCLC = estimate2) %>%
  pivot_longer(-sc, names_to = 'samp', values_to = 'p') %>%
  ggplot(aes(x = sc, y = p * 100, fill = samp)) +
  geom_hline(yintercept = 8, color = 'grey70') +
  geom_col(position = 'dodge', alpha = .7) +
  geom_point(aes(x = sc, y = fc * b + a), data = distinct(d, sc, fc), 
             inherit.aes = F, color = c2) +
  geom_signif(aes(xmin = x - .22, xmax = x+ .22, y_position = y, annotations = symb, group = sc), 
              manual = T, data = d, inherit.aes = F, tip_length = 0, vjust = -0) +
  scale_fill_manual(values = clrs) +
  scale_y_continuous(expand = expansion(0), limits = c(0,20),
                     breaks = seq(0, 16, 4),
                     sec.axis = sec_axis(~(.-a)/b, name = 'log2(GSCLC/GSC)',
                                         breaks = seq(-.2,.2,.1))) +
  labs(x = 'Subcompartment', y = '% of genome') +
  facet_grid(.~'Propotion differences') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y.left = element_line(color= 'black'),
        axis.line.y.right = element_line(color = c2),
        axis.text = element_text(color = 'black', size = 11),
        legend.title = element_blank(),
        legend.position = c(.5,1),
        legend.justification = c(.5,1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks.y.right = element_line(color = c2),
        axis.text.y.right = element_text(color = c2),
        axis.title.y.right = element_text(color = c2),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.line.x = element_line(color = 'black'))

fx <- function(dat, idx) {
  i <- which(dat[idx,1])
  j <- which(dat[idx,2])
  return(length(intersect(i, j)) / (length(i) + length(j) - length(intersect(i, j))))
}

jacc <- unique(lbs$GSC) %>%
  na.omit() %>%
  sort() %>%
  setNames(.,.) %>%
  lapply(function(x) {
    a <- replace_na(lbs$GSC == x, F)
    b <- replace_na(lbs$GSCLC == x, F)
    bt <- boot(cbind(a, b), fx, 1e3)
    boot.ci(bt, type = 'perc') %>%
      .$percent %>%
      as.data.frame() %>%
      `colnames<-`(c('level', 'i1', 'i2', 'v1', 'v2')) %>%
      mutate(ji = bt$t0)
  }) %>%
  bind_rows(.id = 'sc')
p4 <- jacc %>%
  ggplot(aes(x = sc, y='',fill = ji)) + 
  geom_tile() +
  scale_fill_distiller('Jaccard index', palette = 'Purples', direction = 1, breaks = c(.3,.6),
                       guide = guide_colorbar(title.position = 'left', barheight = 5, barwidth = .5, 
                                              title.hjust = 1)) +
  coord_cartesian(expand = F) +
  facet_grid(~'Subcompartment concordance') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 13, color = 'black', face = 'bold'),
        axis.line.x = element_line(color = 'black'),
        axis.title = element_blank(),
        #axis.text = element_text(size = 11, color = 'black')
        axis.text = element_blank(),
        legend.title = element_text(angle = 90),
        legend.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        legend.justification = 'top')

p3 <- lbs[,c('GSC', 'GSCLC')] %>%
  dplyr::count(GSCLC, GSC) %>% 
  na.omit() %>%  mutate(n = 100 * n / sum(n)) %>% 
  ggplot(aes(x = GSC, y = GSCLC)) +
  geom_tile(aes(fill = n)) + 
  scale_fill_gradientn('% of genome', breaks = c(0, 5, 10),
                       colors = plasma(25)) +
  coord_cartesian(expand = F) +
  guides(fill = guide_colorbar(title.position = 'left',
                               barheight = 5,
                               barwidth = .5)) +
  facet_grid(. ~'Subcompartment transitions') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        legend.title = element_text(angle = 90),
        #plot.margin = margin(b = 40),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.direction = 'vertical',
        strip.text = element_text(color = 'black', size = 13, face = 'bold'))

pd5 <- lbs[,c('GSC', 'GSCLC')] %>%
  filter(GSCLC != GSC) %>% 
  dplyr::count(GSCLC, GSC) %>% 
  na.omit() %>% 
  mutate(n = 100 * n / sum(n),
         across(c(GSC, GSCLC), ~as.numeric(factor(.)), .names = '{.col}.idx'),
         higher = ifelse(GSC.idx < GSCLC.idx, 'GSC', 'GSCLC')) %>%
  rowwise() %>% 
  mutate(nm = paste(sort(c(GSC, GSCLC)), collapse=' ')) %>%
  ungroup() %>%
  dplyr::select(nm, higher, n) %>%
  pivot_wider(names_from = 'higher', values_from = 'n') %>%
  mutate(lfc = log2(GSCLC / GSC),
         x = 'Lower / upper')

anns <- pd5 %>% slice_max(lfc, n = 1) %>%
  mutate(p = wilcox.test(pd5$lfc, alternative = 'less')$p.value %>%
           signif.num() %>%
           as.character())

p5 <- ggplot(pd5, aes(x = 'Lower / upper', y = lfc)) +
  geom_boxplot() +
  geom_text(aes(label = p), data = anns, vjust =  0, size = 5) +
  geom_hline(yintercept = 0, color = 'grey70', alpha = .5) +
  scale_y_continuous(position = 'right') +
  labs(y = 'log2(# more active in GSCLC / GSC)') +
  coord_cartesian(clip = 'off') +
  theme(axis.text = element_text(size = 11, color = 'black'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        #axis.text.x = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.title.x = element_blank(),
        axis.ticks = element_blank())

{ wrap_plots(p4, p3 + theme(axis.title.x = element_text(margin = margin(t = -20, unit = "pt"))),
           nrow = 2 ,heights = c(1, 8)) %>%
  wrap_plots(p5, p1 + theme(axis.title.x = element_text(margin = margin(t = -20, unit = "pt"))),
             nrow = 1, widths = c(10, 1,10)) &
  theme(plot.background = element_blank())  } %>%
  ggsave('f6_fg.pdf', ., height = 4.54, width = 8.5)
  
ggsave('f6_leg_bot.pdf', cowplot::get_legend(p3 + theme(legend.position = 'left', legend.justification = c(.5,.5))),height=4.54)
ggsave('f6_leg_top.pdf', cowplot::get_legend(p4 + theme(legend.position = 'left', legend.justification = c(.5,.5))),height=4.54)

  
