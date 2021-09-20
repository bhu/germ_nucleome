library(data.table)
library(tidyverse)
library(plotly)
library(ggtext)
library(broom)
library(cowplot)
library(patchwork)
library(pals)
library(ggnewscale)
odr <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC')
odrr <- c('repressed', 'CTCF', 'enhancer',  'bivalent', 'promoterCTCF', 'promoter')
cclrs <- setNames(kelly(22)[2:7], odrr)
clrs <- dark <- setNames(tableau20(12)[seq(1, 10, 2)], odr)
light <- setNames(tableau20(12)[seq(2, 10, 2)], odr)

load('../data/clust/open.rda')
d <- d %>%
  filter(cl != -1 & type %in% odr) %>% 
  mutate(V1 = cl)

p1 <- d %>%
  mutate(cl = case_when(V1 %in% c(1,10,11) ~ 'promoterCTCF',
                        V1 %in% c(7,8,9,12) ~ 'bivalent',
                        V1 %in% c(6,13) ~ 'promoter',
                        V1 %in% c(4,5,14,17) ~ 'enhancer',
                        V1 %in% c(15,16,18) ~ 'repressed',
                        V1 %in% c(0,2,3) ~ 'CTCF',
                        T ~ NA_character_) %>%
           factor(odrr),
         V1 = V1 + 1) %>% 
  group_by(cl, V1) %>%
  summarise_if(is.numeric, median) %>%
  select(lab = cl, clu = V1, CTCF, Rad21, ATAC, K4me1, K4me3, K27ac,
         K27me3, Ring1b, H2Aub, K36me2, K36me3, K9me2, K9me3) %>%
  group_by(lab) %>% 
  mutate(mi = min(clu)) %>% 
  ungroup() %>%
  #arrange(mi, lab, clu) %>%
  arrange(lab, clu) %>%
  select(-mi) %>%
  pivot_longer(-c('lab','clu'), names_to = 'mark', values_to = 'med') %>%
  mutate(mark = case_when(grepl('^K', mark)~ paste0('H3', mark), T ~ mark),
         mark = fct_inorder(mark),
         clu = fct_inorder(as.character(clu)),
         xl = fct_inorder(sprintf("<span style='color:%s'>%s</span>", cclrs[lab], clu))) %>%
  ggplot(aes(x = xl, y = mark)) +
  geom_tile(aes(fill = med)) +
  scale_fill_gradientn(colors = rev(brewer.brbg(25)), limits = c(-4,4), breaks = c(-4, 0, 4),
                       name = 'Enrichment', oob = scales::squish,
                       guide = guide_colorbar(barheight = .5, barwidth = 4,
                                              title.position = 'top')) +
  scale_x_discrete(expand = expansion(add = .5)) +
  scale_y_discrete(expand = expansion(0)) +
  labs(y = 'Mark') +
  theme(axis.text.x = element_markdown(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        legend.text = element_text(size = 11),
        legend.direction = 'horizontal',
        axis.text = element_text(size = 11),
        axis.text.y = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_blank())

or <- odr %>%
  lapply(function(x) {
    n1 <- d$type == x
    unique(d$V1) %>%
      lapply(function(y) {
        n2 <- d$V1 == y
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
  mutate(lor = log10(estimate))

slope <- 1e4
const <- 3e4

pd <- or %>%
  mutate(V1 = y,
         lab =case_when(V1 %in% c(1,10,11) ~ 'promoterCTCF',
                        V1 %in% c(7,8,9,12) ~ 'bivalent',
                        V1 %in% c(6,13) ~ 'promoter',
                        V1 %in% c(4,5,14,17) ~ 'enhancer',
                        V1 %in% c(15,16,18) ~ 'repressed',
                        V1 %in% c(0,2,3) ~ 'CTCF',
                        T ~ NA_character_) %>%
           factor(odrr),
         V1 = V1 + 1) %>%
  group_by(lab) %>% 
  mutate(mi = min(V1)) %>%
  ungroup() %>%
  #arrange(mi, lab, V1) %>%
  arrange(lab, V1) %>%
  select(-mi) %>%
  mutate(clu = fct_inorder(as.character(V1)),
         xl = fct_inorder(sprintf("<span style='color:%s'>%s</span>", cclrs[lab], clu)),
         samp = factor(x, odr))

rcts <- distinct(pd, clu, xl) %>% 
  mutate(alt = as.character(1:n() %% 2),
         x = as.numeric(clu))

p2 <- ggplot() +
  scale_x_discrete(expand = expansion(add = .5)) +
  geom_hline(yintercept = 3e4, color = 'grey70') +
  geom_rect(aes(xmin = x - .5, xmax = x + .5, ymin = -Inf, ymax = Inf, fill = alt), 
            data = rcts, show.legend = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  new_scale_fill() +
  geom_col(aes(x = xl, y = a, fill = samp), position = 'dodge', data = pd) +
  geom_point(aes(x = xl, y = lor * slope + const, color = samp), data = pd, show.legend = F) +
  scale_y_continuous(expand = expansion(c(0,.05)),
                     sec.axis = sec_axis(~., 'Odds ratio for label vs cell type', breaks = (1:5) * 1e4,
                                         labels = c('0.01', '0.1', '1', '10', '100'))) +
  scale_fill_manual(values = light, name = 'Cell type', labels = function(x) sub('ESC', 'mESC', sub('PGC', ' mPGC', x))) +
  scale_color_manual(values = dark) +
  facet_grid(.~'Characterization of open site clusters') +
  labs(y = '# of open sites in cluster') +
  theme(axis.text.x = element_markdown(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'))
leg <- data.frame(c = cclrs) %>% 
  rownames_to_column('state') %>% 
  filter(state != 'noise') %>%
  mutate(state = fct_inorder(state)) %>%
  ggplot(aes(x = c, y =state, color = state, group = 1)) + geom_line(size = 1) + 
  scale_color_manual(values = cclrs) +
  guides(color = guide_legend(nrow = 1, byrow = T)) +
  theme(legend.position = 'bottom',
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank())
                   
{wrap_plots(p2, p1, get_legend(leg), ncol = 1, heights = c(1,1,.2)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('f3_b.pdf', ., height = 6.5, width = 9)

