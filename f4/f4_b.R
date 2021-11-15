library(tidyverse)
library(pals)


load('../data/tracks/inpnorm.50kb.rda')

d <- d50 %>% 
  select(-c(chr, start, end)) %>%
  {split(select(., -samp), .$samp)}

smps <- c('ESC', 'EpiLC','d2PGCLC','d4c7PGCLC')
mrks <- c('ATAC','CTCF','Rad21',
          'Ring1b','H2Aub','K27me3', 'K27ac', 'K4me1','K4me3',
          'K36me2','K36me3','K9me2','K9me3', 'Laminb1')

lapply(d, function(x) {
  {x - d$GSC} %>%
    cor(method = 'spearman') %>%
    .['PC1',] %>%
    data.frame(y = .) %>%
    rownames_to_column('mark') 
}) %>%
  bind_rows(.id = 'samp') %>%
  filter(!(mark %in% c('Stag1', 'Stag2'))) %>%
  filter(samp %in% smps & mark != 'PC1') -> dd

dd %>%
  mutate(samp = factor(samp, smps),
         mark = factor(mark, mrks)) %>%
  na.omit() %>%
  arrange(samp, mark) %>%
  mutate(mark = sub('^K', 'H3K', mark) %>%
           sub('^', '\u0394', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x = samp, y = mark, fill = y)) +
  geom_tile() +
  scale_fill_gradientn(expression('Spearman\'s' ~ rho), breaks = c(-.5, 0, .5),
                       colors = coolwarm(25), limits = c(-.62, .62),
                       guide = guide_colorbar(title.position = 'top', title.hjust = .5,
                                              barheight = .5, barwidth = 5)) +
  scale_x_discrete(labels = function(x) sub('^', 'GSC - ', sub('ESC', 'mESC', sub('PGC', ' mPGC', x)))) +
  facet_grid(.~'bold(Corr.~vs~Delta*comp.~score)', labeller = label_parsed) +
  coord_cartesian(expand = F) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        legend.justification = c(1,2),
        plot.margin = margin(5,5,35,5),
        legend.direction = 'horizontal',
        legend.position = c(1,0)) -> p 
ggsave('f4_b.pdf', p, height = 4.67, width = 3, device = cairo_pdf, bg = 'transparent')
