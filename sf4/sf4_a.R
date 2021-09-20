library(tidyverse)
library(pals)
library(LOLA)
library(GenomicRanges)
regionDB <- loadRegionDB('../data/resources/reps')
regs <- regionDB$regionGRL %>% setNames(sub('.bed', '', regionDB$regionAnno$filename)) 
load('../data/tracks/inpnorm.50kb.rda')

gr <- distinct(d50, chr, start, end) %>% 
  mutate(start = start + 1) %>% 
  makeGRangesFromDataFrame()

d <- lapply(regs, function(x) {
  findOverlaps(gr, x) %>%
    as("List") %>%
    sapply(length) %>%
    tibble(v = .) %>%
    cbind(filter(d50, samp == 'GSC')) %>%
    {cor.test(.$Laminb1, .$v, method = 'pearson')} %>%
    {tibble(estimate = .$estimate, lower = .$conf.int[1], 
            upper = .$conf.int[2], p = .$p.value)} 
}) %>%
  bind_rows(.id = 'fam') %>%
  arrange(estimate) %>% 
  mutate(fam = fct_inorder(fam))

rcts <- d %>%
  mutate(x = as.numeric(fam),
         alt = (x %% 2) == 0)

ggplot(d, aes(x = fam, y = estimate, ymin = lower, 
             ymax = upper, color = estimate)) +
  scale_x_discrete(expand = expansion(0)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_rect(aes(xmin = x - .5, xmax = x + .5, ymin = -Inf, ymax = Inf, fill = alt),
            data = rcts, inherit.aes = F, show.legend = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  geom_point() +
  geom_linerange() +
  scale_color_gradientn('Pearson\'s r', colors = coolwarm(25), limits = c(-.3,.3),
                        breaks = c(-.2, 0, .2)) +
  guides(color = guide_colorbar(title.position = 'top',
                                title.hjust = .5,
                                barheight = .5)) +
  labs(x = 'Repeat family', y = 'Correlation') +
  facet_grid(.~'Correlation between repeat density and Lamin B1 enrichment in GSCs') +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.text = element_text(size = 11),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.line.x = element_line(color = 'black'),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.key = element_blank()) +
  ggsave('sf4_a.pdf', height = 3.5, width = 13)

