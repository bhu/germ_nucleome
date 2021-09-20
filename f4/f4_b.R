library(tidyverse)
library(pals)


load('../data/tracks/inpnorm.50kb.rda')

d50 %>% 
  filter(samp %in% c('ESC', 'GSC')) %>% 
  select(-c(chr, start, end)) %>%
  {split(select(., -samp), .$samp)} %>% 
  {.$ESC - .$GSC} %>%
  cor(method = 'spearman') %>%
  .['PC1',] %>%
  data.frame(y = .) %>%
  rownames_to_column('x') %>%
  arrange(y) %>%
  filter(!(x %in% c('PC1','Stag1','Stag2'))) %>%
  mutate(x = sub('^K', 'H3K', x) %>%
           fct_inorder(),
         idx = as.numeric(x),
         alt = as.character(idx %% 2)) %>%
  ggplot(aes(x = x, y = y)) +
  scale_x_discrete(expand = expansion(0)) +
  geom_rect(aes(xmin = idx - .5, xmax = idx + .5, fill = alt, 
                ymin = -Inf, ymax = Inf), show.legend = F) +
  scale_fill_manual(values = c('#ffffff00', '#11111111')) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_linerange(aes(ymin = 0, ymax = y, color = y)) +
  geom_point(aes(color = y)) +
  scale_color_gradientn(colors = coolwarm(25), limits = c(-.65,.65)) +
  coord_flip() +
  facet_grid(.~'bold(Corr.~vs~Delta*comp.~score)', labeller = label_parsed) +
  labs(y = expression('Spearman\'s'~rho), x = 'GSC - ESC') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.clip = 'off',
        strip.text = element_text(color = 'black', size = 13),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
        axis.line.x = element_line(color = 'black')) +
  ggsave('f4_b.pdf', height = 4.6, width = 2.8, device = cairo_pdf, bg = 'transparent')
