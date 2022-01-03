library(tidyverse)
library(pals)
library(ppcor)

load('../data/tracks/inpnorm.50kb.rda')

samps <- c('ESC', 'EpiLC', 'd2PGCLC','d4c7PGCLC', 'GSC')
marks<- c('CTCF', 'Rad21', 'K4me1', 'K4me3', 'K27ac',
          'K27me3', 'Ring1b', 'H2Aub', 'K36me2', 'K36me3', 'K9me2', 'K9me3')
rnm <- function(x) {
  case_when(x == 'ESC' ~ 'mESC',
            x == 'd2PGCLC' ~ 'd2 mPGCLC',
            x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
            T ~ x)
}
  
setdiff(samps, 'd4c7PGCLC') %>%
  lapply(function(x) {
    dm <- list(d50[d50$samp == x,],
         d50[d50$samp == 'd4c7PGCLC',]) %>%
      lapply(function(y) {
        y %>% dplyr::select(-c(samp, chr, start, end, PC1, Laminb1))
      }) %>%
      {.[[1]] - .[[2]]} 
    cm <- pcor(dm)$estimate
    colnames(cm) <- rownames(cm) <- colnames(dm)
    cm[,'ATAC'] %>% 
      data.frame(v = .) %>%
      rownames_to_column('mark') %>%
      mutate(samp = x)
  }) %>%
  bind_rows() %>%
  filter(mark != 'ATAC') %>%
  mutate(samp = factor(rnm(samp), rnm(samps)),
         mark = factor(mark, marks)) %>%
  na.omit() %>%
  arrange(samp, mark) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           sub('^', '\u0394', .) %>%
           fct_inorder(),
         samp = paste0(samp, 's') %>%
           fct_inorder()) %>%
  ggplot(aes(x = samp, y = mark, fill = z)) +
  geom_tile(aes(fill = v)) +
  coord_cartesian(expand = F) +
  scale_fill_gradientn(expression('Partial Pearson correlation vs'~Delta*'ATAC-seq'),
                       colors = coolwarm(25), limits = c(-1,1), breaks = -1:1,
                       guide = guide_colorbar(barheight = 2, barwidth = .5, title.position = 'bottom', title.vjust = 0)) +
  facet_grid(.~'Determinants of openness') +
  xlab('Cell type compared\nagainst d4c7PGCLCs') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(angle = -90),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13,face = 'bold'))  +
  ggsave('f2_e.pdf', height = 4.96, width = 3.5, bg = 'transparent', device = cairo_pdf)


