library(tidyverse)
library(matrixStats)
library(ggrastr)
library(ggnewscale)
library(ggforce)
library(pals)
library(ggrepel)

load("../data/peaks/ATAC/public.fpkm.rda")
rnm <- function(x) {
  sub('PGCLC', ' mPGCLC', sub('ESC', 'mESC', x))
}
samps <- c("ESC", "EpiLC", "d2PGCLC", "d4PGCLC", "d4c7PGCLC", "GSC", "MEF") %>% rnm()
clrs <- setNames(tableau20(20)[c(1, 3, 5, 6, 7, 9, 17)], samps)
clrs[rnm('d4PGCLC')] <- '#81642A'

res <- atac.my.pub.pks %>% 
  select(matches('BDF121|MEF_Gia|GSC_AAG')) %>%
  slice_max(rowVars(as.matrix(.)), n = 1e4) %>%
  mutate(idx = 1:n()) %>%
  pivot_longer(-idx, names_to = 'samp', values_to = 'v') %>%
  mutate(samp = sub('_[0-9]$', '', samp)) %>%
  group_by(idx, samp) %>%
  summarise(v = mean(v), .groups = 'drop') %>%
  pivot_wider(names_from = 'samp', values_from = 'v') %>%
  select(-idx) %>%
  t() %>%
  prcomp(scale. = T, center = T)

comps <- summary(res)$importance[2,] * 100

axes <- c(1, 2) %>%
  setNames(paste0('PC', .))

pd <- res$x %>%
  data.frame() %>%
  rownames_to_column("samp") %>% 
  mutate(samp = sub('_.*', '', samp) %>%
           rnm() %>%
           factor(samps)) %>% 
  arrange(samp) %>%
  rename(!!"ax1" := names(axes)[1], !!"ax2" := names(axes)[2]) %>%
  mutate(xend = lead(ax1), yend = lead(ax2))

p <- ggplot(pd, aes(x = ax1, y = ax2, xend = xend, yend = yend))

for (i in 1:5) {
  if (i < 5) {
    p <- p + geom_link(aes(color = stat(index)), data = pd[i,], 
                        show.legend = F, alpha = 1) +
      scale_color_gradient(low = clrs[i], high = clrs[i + 1]) +
      new_scale_color()
  } else {
    p <- p + geom_link(aes(color = stat(index)), data = pd[i,], 
                        show.legend = F, alpha = .5) +
      scale_color_gradient(low = clrs[i], high = clrs[i + 1]) +
      new_scale_color()
  }
}

p +
  geom_point(aes(color = samp), size = 2) +
  scale_color_manual(values = clrs) +
  labs(x = sprintf('%s (%.1f%%)', names(axes)[1], comps[axes[1]]),
       y = sprintf('%s (%.1f%%)', names(axes)[2], comps[axes[2]])) +
  facet_grid(.~'ATAC-seq PCA') +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 11),
        strip.background = element_rect(fill = NA),
        legend.text = element_text(color = 'black', size = 11),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size =13, face = 'bold'),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed')) -> p
ggsave('sf5_b.pdf', p, width = 4.2, height = 3.82)

