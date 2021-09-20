library(tidyverse)
library(readxl)
library(pals)
library(patchwork)

samps <- c("d4c7PGCLC","GSC",  'GSCLC') 

clrs <- setNames(c('#d62728', '#9467bd', '#499894'), samps)

ms <- read_excel('../data/histone_ratios.xlsx') %>%
  {.[,colSums(is.na(.)) != nrow(.)]} %>% 
  na.omit() %>% 
  {.[,c(1, which(.[1,] == 'Area'))]} %>%
  tail(-1) %>%
  rename(p = 1) %>%
  `names<-`(., sub('\\..*', '', names(.))) %>%
  separate(p, c('pep', 'mod'), '\\ ') %>%
  mutate(pep = sub('H33', 'H3', pep)) %>%
  group_by(pep) %>%
  mutate(across(-mod, function(x) as.numeric(x) / sum(as.numeric(x)))) %>%
  ungroup() %>% 
  filter(grepl('^H3', pep)) %>%
  mutate(mod = strsplit(mod, '(?<=.)(?=K)', perl = T)) %>% 
  unnest(mod) %>%
  dplyr::select(-pep) %>%
  group_by(mod) %>%
  summarise(across(everything(), sum), .groups = 'drop') %>%
  pivot_longer(-mod, names_to = 'samp', values_to = 'a') %>%
  separate(samp, c('type', 'rep'), '_', F)


mods <- c("K4me3" , "K4me1"  , "K36me2", "K36me3",  "K9me2" , "K9me3" , "K27ac" ,"K27me3" )
p1 <- ms %>%
  filter(type %in% c('GSC', 'GSCLC') & mod %in% mods) %>%
  group_by(mod, type) %>% 
  summarise(std = sd(a),
            se = std / n(),
            mu = mean(a), 
            .groups = 'drop') %>%
  pivot_wider(names_from = 'type', values_from = c('mu', 'std', 'se')) %>%
  mutate(fc = mu_GSCLC / mu_GSC,
         se = sqrt((se_GSCLC/mu_GSCLC)^2 + (se_GSC/mu_GSC)^2)) %>% 
  arrange(fc) %>%
  ggplot(aes(x = fct_inorder(mod), y = fc)) +
  geom_hline(yintercept = 1, color = 'grey70') +
  geom_point() +
  geom_linerange(aes(ymin = fc - se, ymax = fc + se)) +
  scale_y_continuous(breaks = c(1,1.4)) +
  labs(x = 'H3 modification', y = 'GSCLC / GSC') +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major.y = element_line(color = 'grey90'),
        panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed')) 

p2 <- ms %>%
  filter(type %in% samps & grepl('K27me3|K9me2|K9me3', mod)) %>%
  mutate(type = factor(type, samps),
         mod = paste0('H3', mod)) %>%
  ggplot(aes(x = type, y = a * 100, color = type)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), size = .7, fatten = 1,
               geom = 'pointrange', mapping = aes(color = type)) +
  facet_grid(.~mod, scales = 'free_y') +
  scale_color_manual(values = clrs) +
  scale_x_discrete(labels = function(x) sub('PGC', ' mPGC', x)) +
  ylab('Abundance (%)') +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black',size = 13, face = 'bold'),
        axis.title.x = element_blank(),
        strip.clip = 'off',
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank()) +
  ylab('Abundance (%)') 

wrap_plots(p1 + theme(axis.title.x = element_text(margin = margin(t = -45))),
           p2, nrow = 1, widths = c(1,3)) &
  theme(plot.background = element_blank()) &
  ggsave('f7_h.pdf', height = 2.6, width = 5.4)

