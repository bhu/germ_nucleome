library(data.table)
library(tidyverse)
library(ggtext)
library(rasterly)
library(cowplot)
library(patchwork)
library(pals)
library(broom)
library(ggnewscale)

odr <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
clrs <- dark <- setNames(tableau20(12)[seq(1, 10, 2)], odr)
light <- setNames(tableau20(12)[seq(2, 10, 2)], odr)
cclrs <- kelly(22)[c(8,10,11,12,15,14,20)]
load('../data/clust/promoter.rda')

d <- d %>%
  filter(cl != -1 & type != 'GSCLC') %>%
  {.[.$gene %in% Reduce(intersect, split(.$gene, .$type)),]} 

clabs <- c('Repressed', 'Polycomb', 'Active', 'Bivalent', 'Cohesin-CTCF',
           'Readthrough', 'Enhancer-like')

pd <- odr %>%
  lapply(function(x) {
    n1 <- d$type == x
    unique(d$cl) %>%
      lapply(function(y) {
        n2 <- d$cl == y
        m <- tibble(a = sum(n1 & n2),
                    b = sum(n1) - a,
                    c = sum(n2) - a,
                    d = nrow(d) - a - b - c)
        fisher.test(matrix(unlist(m), 2,2), alternative = 'two.sided') %>% 
          tidy() %>% .[,1:4] %>% 
          mutate(p = fisher.test(matrix(unlist(m), 2,2), alternative = 'greater')$p.value) %>% 
          cbind(m) %>%
          mutate(x = x, y = y)
      }) %>%
      bind_rows()
  }) %>% bind_rows() %>%
  mutate(lor = log2(estimate),
         samp = factor(x, rev(odr)),
         lo = log2(conf.low),
         hi = log2(conf.high)) %>%
  arrange(samp, -y) %>%
  mutate(xl = fct_inorder(sprintf("<span style='color:%s'>%s</span>", cclrs[y], clabs[y])))

rcts <- distinct(pd, y, xl) %>% 
  mutate(alt = as.character(1:n() %% 2))

slope <- 1e4
const <- 2e4


  sec.axis = sec_axis(~., 'Odds ratio for label vs cell type', breaks = (1:3) * 1e4,
                      labels = c('0.5', '1', '2'))
  
p1 <- ggplot() +
  scale_x_discrete(expand = expansion(add = .5)) +
  coord_flip(ylim = c(0,3e4)) +
  geom_rect(aes(xmin = y - .5, xmax = y + .5, ymin = -Inf, ymax = Inf, fill = alt), 
            data = rcts, show.legend = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  new_scale_fill() +
  geom_col(aes(x = xl, y = a, fill = samp), position = 'dodge', data = pd, alpha = .8) +
  scale_y_continuous(expand = expansion(c(0,.05))) +
  scale_fill_manual(values = light, name = 'Cell type') +
  scale_color_manual(values = dark) +
  #facet_grid('Proportion of promoters in each cluster' ~ .) +
  labs(y = '# of promoters in cluster',
       x = 'Promoter cluster') +
  theme(axis.text.y = element_markdown(),
        plot.background = element_blank(),
        axis.line.y = element_line(color = 'black'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.text = element_text(color = 'black',size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.position = 'none',
        legend.title = element_blank(),
        axis.title.y = element_text(vjust = 2),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'))


p2 <- ggplot() +
  scale_x_discrete(expand = expansion(add = .5)) +
  coord_flip() +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_rect(aes(xmin = y - .5, xmax = y + .5, ymin = -Inf, ymax = Inf, fill = alt), 
            data = rcts, show.legend = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  new_scale_fill() +
  geom_point(aes(x = xl, y = lor , color = samp, group = samp),
             data = pd, show.legend = F, position = position_dodge(.9)) +
  geom_linerange(aes(x = xl, ymin = lo , ymax = hi , color = samp, group = samp),
                 data = pd, show.legend = F, position = position_dodge(.9)) +
  scale_y_continuous(expand = expansion(c(0,.05)), breaks = c(-1, 0, 1), labels = c('0.5', '1', '2')) +
  scale_fill_manual(values = light, name = 'Cell type') +
  scale_color_manual(values = dark) +
  #facet_grid('Proportion of promoters in each cluster' ~ .) +
  labs(y = 'Odds ratio for cluster vs cell type') +
  theme(axis.text.y = element_blank(),
        plot.background = element_blank(),
        axis.line.y = element_line(color = 'black'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.text = element_text(color = 'black',size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.position = 'none',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed')) 

{wrap_plots(p1,p2,nrow=1) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf7_d.pdf', ., height = 4.77, width = 8)
 
