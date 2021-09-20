library(data.table)
library(tidyverse)
library(geometry)
library(pals)
library(ggh4x)
library(ggpubr)

load('../data/compscore/100kb.rda')
samps <- c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC','GSCLC')
sampss <- c('ESC (2i)','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
            'ESC (serum)', 'NPC','CN',
            'Day 6','Day 4','Day 2','B\u03B1',
            'CM','CPC','PCM','VCM')
rnm <- function(x) {
  #x %>%
  #  sub('C ', 'Cs ', .) %>%
  #  sub('^(.*[CMN])$', '\\1s', .) 
  x <- sub('ESC', 'mESC', x)
  x <- case_when(x == 'ESC' ~ 'mESC',
                 x == 'd2PGCLC' ~ 'd2 mPGCLC',
                 x == 'd4c7PGCLC' ~ 'd4c7 mPGCLC',
                 T ~ x)
}
sampss <- rnm(sampss)
clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(12)[-4]), sampss)
light <- paste0(clrs, 'aa')
dd <- list.files('../data/3dmods', pattern = 'xyz', full.names = T) %>% 
  setNames(., sub('.xyz', '', basename(.))) %>% 
  tibble(f = .) %>%
  mutate(s = sub('.xyz', '', basename(f))) %>%
  separate(s, c('x', 'chr'), '_') %>% 
  mutate(study = case_when(x %in% as.character(1:6) ~ 'Nagano',
                           x %in% c('ES','CN','NPC') ~ 'Bonev2017',
                           x %in% c(paste0('D', c(2,4,6,8)), 'PSC', 'B', 'Ba') ~ 'Stadhouders2018',
                           T ~ 'Zhang2019'),
         x = case_when(x %in% as.character(1:6) ~ samps[as.integer(x)],
                       T ~ sub('^ES$', 'ESC', x)),
         chr = case_when(study == 'Zhang2019' ~ sub('23', 'X', chr),
                         T ~ sub('20', 'X', chr)),
         gn = case_when(study == 'Zhang2019' ~ 'hg38',
                        T ~ 'mm10')) %>%
  split(.,.$study) %>%
  lapply(function(x) {
    split(x, x$x) %>%
      lapply(function(y) {
        split(y, y$chr) %>%
          lapply(function(z) {
            dat <- fread(z$f)
            k <- mutate(bins[[z$gn]], s = mats[[z$gn]][[z$study]][[z$x]]) %>%
              dplyr::filter(chrom == paste0('chr', z$chr)) %>%
              na.omit() %>%
              pull(start) %>%
              {c(., as.integer(. + 5e4))}
            m <- dat[dat$V1 %in% k, 2:4] %>% as.matrix()
            d <- dist(m) %>% as.matrix()
            d[upper.tri(d, diag = T)] <- NA
            rg2.raw <- sum(d^2,na.rm = T)*2 / (2 * (nrow(d)^2))
            rh_1.raw <- sum(1/d,na.rm = T)*2/nrow(d)^2
            ratio.raw <- sqrt(rg2.raw)*rh_1.raw
            l <- sapply(1:(nrow(m) - 1), function(i) {
              sqrt(sum((m[i,] - m[i + 1,])^2))
            }) %>% sum()
            d <- d / l
            rg2 <- sum(d^2,na.rm = T)*2 / (2 * (nrow(d)^2))
            rh_1 <- sum(1/d,na.rm = T)*2/nrow(d)^2
            ratio <- sqrt(rg2)*rh_1
            tibble(v = convhulln(m / l, output.options = 'FA')$vol,
                   v.raw = convhulln(m, output.options = 'FA')$vol,
                   rg.raw = sqrt(rg2.raw),
                   rh.raw = 1/rh_1.raw,
                   rg = sqrt(rg2),
                   rh = 1/rh_1,
                   ratio = ratio,
                   ratio.raw = ratio.raw)
          }) %>%
          bind_rows(.id = 'chr')
      }) %>%
      bind_rows(.id = 'x')
  }) %>%
  bind_rows(.id = 'study') 

ddd <- dd %>%
  mutate(samp = case_when(x == 'Ba' ~ 'B\u03B1',
                          x == 'ESC' & study == 'Nagano' ~ 'ESC (2i)',
                          x == 'ESC' & study == 'Bonev2017' ~ 'ESC (serum)',
                          T ~ sub('^D', 'Day ', x)) %>%
           rnm() %>%
           factor(sampss),
         study = c(Nagano = 'Germline', Bonev2017 = 'Neural',
                   Stadhouders2018 = 'Reprogram',
                   Zhang2019 = 'Cardiac')[study] %>%
           factor(c('Germline', 'Neural', 'Reprogram', 'Cardiac'))) %>%
  arrange(samp, chr) %>%
  na.omit()


anns <- split(ddd, ddd$study) %>%
  lapply(function(d) {
    compare_means(v ~ samp, paired = T, data = d, method = 'wilcox.test') %>%
      mutate(x1 = as.numeric(factor(group1, sampss)),
             x2 = as.numeric(factor(group2, sampss))) %>%
      dplyr::filter(x2 - x1 == 1) %>%
      mutate(idx = 1:n(),
             y = case_when(d$study[1] == 'Germline' ~ 10^(-4.1 + idx * .1),
                           d$study[1] == 'Neural' ~ 10^(-4.3 + idx * .1),
                           d$study[1] == 'Reprogram' ~ 10^(-3.7 - idx * .1),
                           d$study[1] == 'Cardiac' ~ 10^(-5.55 - idx * .1)))
  }) %>%
  bind_rows(.id = 'study') %>%
  mutate(study = fct_inorder(study))

ggplot(ddd, aes(x = samp, y = v)) +
  geom_line(aes(group = chr), color = 'grey70') +
  geom_boxplot(outlier.colour = NA, aes(color = samp, fill = samp)) +
  stat_summary(geom = "crossbar", width = 0.7, fatten = 0, color = "white", 
               fun.data = function(x) return(c(y = median(x), ymin = median(x), ymax = median(x)))) +
  geom_segment(aes(x = group1, xend = group2, y = y, yend = y),
               data = anns, inherit.aes = F) +
  geom_text(aes(x = idx + .5, y = y, label = p.signif),
            data = anns, inherit.aes = F, vjust = -.3) +
  scale_color_manual(values = paste0(clrs, '99')) +
  scale_fill_manual(values = paste0(clrs, 'cc')) +
  #geom_signif(aes(xmin = x1, xmax = x2, annotations = p.signif), y = .000001,
  #            manual = T, data = anns, inherit.aes = F, tip_length = 0, vjust = -0) +
  #stat_pvalue_manual(anns, label = 'p.signif', tip.length = 0) +
  labs(y = 'Volume of chromosome territories (a.u.)') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
               labels = scales::trans_format("log10", scales::math_format(10^.x)),
               expand = expansion(c(.01,.01))) +
  #coord_cartesian(ylim = c(5e-6, 1.5e-4)) +
  facet_grid(.~study, scales = 'free', space = 'free') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        #axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = "bold")) +
  ggsave('f1_j.pdf', height = 3, width = 4.85, device = cairo_pdf, bg = 'transparent')

