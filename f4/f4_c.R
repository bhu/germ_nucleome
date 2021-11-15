library(tidyverse)
library(ggh4x)
library(ggrastr)
library(pals)

load('../data/tracks/inpnorm.50kb.rda')

p <- d50 %>% 
  filter(samp %in% c('EpiLC','GSC')) %>%
  filter(chr == 'chr3') %>%
  split(., .$samp) %>%
  lapply(function(x) {
    lapply(x[,c('Laminb1', 'PC1')], function(y) {
      smooth.spline((x$start + x$end) / 2, y) %>% 
        predict() %>% 
        {tibble(x = .$x, y = .$y)}
    }) %>% bind_rows(.id = 'mark')
  }) %>% bind_rows(.id = 'samp') %>%
  mutate(mark = c(PC1 = 'Compartment score', Laminb1 = 'Lamin B1 ChIP-seq')[mark] %>%
           factor(),
         pos = ifelse(y > 0, y, 0),
         neg = ifelse(y < 0, y, 0)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_line(alpha = .5) +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'forestgreen', alpha = .5,
              data = ~ subset(., mark == 'Lamin B1 ChIP-seq')) +
  geom_ribbon(aes(ymin = 0, ymax = pos), fill = 'firebrick', alpha = .5,
              data = ~ subset(., mark == 'Compartment score')) +
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'grey70', alpha = .5,
              data = ~ subset(., mark == 'Lamin B1 ChIP-seq')) +
  geom_ribbon(aes(ymin = 0, ymax = neg), fill = 'steelblue', alpha = .5,
              data = ~ subset(., mark == 'Compartment score')) +
  geom_rug(aes(color = ifelse(y > 0, "A", "B")), sides = "b",
           data = ~ subset(., mark == 'Compartment score')) +
  geom_rug(aes(color = ifelse(y > 0, "C", "D")), sides = "b",
          data = ~ subset(., mark == 'Lamin B1 ChIP-seq')) +
  geom_vline(xintercept = Inf, color = 'black', size = 1) +
  scale_color_manual(values = c(A = 'firebrick', B = 'steelblue', C = 'forestgreen', D = 'grey70')) +
  facet_nested(mark + samp ~., scales = 'free_y', nest_line = T) +
  ylab('Compartment score or log2 enrichment over input') +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5)),
    scale_y_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5)),
    scale_y_continuous(limits = c(-.25, .25), breaks = c(-.2, 0, .2)),
    scale_y_continuous(limits = c(-.25, .25), breaks = c(-.2, 0, .2))
  )) +
  scale_x_continuous('chr3', breaks = c(0, 5e7, 1e8, 1.5e8),
                     limits = c(0, 1.6e8),
                     labels = c('0', '50mb', '100mb', '150mb'),
                     expand = expansion(0)) +
  theme(plot.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'))

g <- ggplot_gtable(ggplot_build(p))
strip <- which(grepl('strip-r', g$layout$name))[3:6]
fills <- tableau20()[c(3,9,3,9)]
k <- 1
for (i in strip) {
  #j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- fills[k]
  k <- k+1
}

ggsave('f4_c.pdf', g, height = 4.6, width = 7)

