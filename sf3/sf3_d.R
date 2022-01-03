library(tidyverse)
library(pals)

load('../data/compscore/100kb.rda')

samps <- c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC',
           'ESC (serum)', 'Neural progenitors', 'Cortical neurons',
           'Day 6','Day 4','Day 2','B\u03B1',
           'Cardiac mesoderm', 'Cardiac progenitors',
           'Primitive cardiomyocytes', 'Ventricular cardiomyocytes')

rnm <- function(x) {
  sub('ESC', 'mESC', sub('PGCLC', ' mPGCLC', x))
}

samps <- rnm(samps)

clrs <- setNames(c(tableau20(12)[seq(1, 9, 2)],
                   stepped3(16)[-4]), samps)

bsz <- 100000

ac <- mats %>%
  Map(function(x, gn) {
    lapply(x, function(y) {
      o <- y %>%
        cbind(bins[[gn]]) %>%
        pivot_longer(-c('chrom','start','end'), 
                    names_to = 'samp', values_to = 'score') %>% 
        dplyr::rename(chr = chrom) %>% 
        filter(!grepl('_[0-9]$', samp)) %>% 
        na.omit() %>% 
        add_count(chr, start) %>% 
        filter(n == max(n)) %>% 
        dplyr::select(-n) 
      
      dm <- o %>%
        mutate(idx = as.integer(start / bsz)) %>%
        arrange(samp, chr, idx) 
      
      dl <- distinct(dm, chr, idx) %>%
        split(., .$chr) %>%
        lapply(function(x) {
          group_by(x, gw = cumsum(c(1, diff(idx) != 1))) %>%
            summarise(n = n(), i = idx[1], j = idx[n]) %>%
            ungroup()
        }) %>%
        bind_rows(.id = 'chr') %>%
        group_by(chr) %>%
        top_n(1, n) %>%
        ungroup() %>%
        mutate(l = n * bsz / 1e6)
      
      split(dm, dm$samp) %>%
        lapply(function(x) {
          split(x, x$chr) %>%
            lapply(function(y) {
              y %>%
                dplyr::filter(y$idx > dl[dl$chr == y$chr[1],]$i &
                                y$idx < dl[dl$chr == y$chr[1],]$j) %>%
                pull(score) %>%
                acf(na.action = na.pass, plot = F,
                    lag.max = dl[dl$chr == y$chr[1],]$n) %>%
                {tibble(lag = .$lag[,,1], c = .$acf[,,1])} %>%
                dplyr::filter(lag >= 0) %>%
                mutate(lag = lag * bsz,
                       c = abs(c))
            }) %>%
            bind_rows(.id = 'chr') %>%
            dplyr::filter(!(chr %in% c('chrY'))) %>%
            group_by(lag) %>%
            summarise(s = sd(c, na.rm = T) / sqrt(20),
                      c = mean(c, na.rm = T),
                      lo = c - s, hi = c + s) 
        }) %>%
        bind_rows(.id = 'x') %>%
        na.omit() %>%
        dplyr::filter(lag != 0) 
      
    }) %>%
      bind_rows(.id = 'study')
  }, ., names(.)) %>%
  bind_rows(.id = 'gn')


ac %>%
  mutate(x = case_when(x == 'Ba' ~ 'B\u03B1',
                       x == 'CM' ~ 'Cardiac mesoderm',
                       x == 'CPC' ~ 'Cardiac progenitors',
                       x == 'PCM' ~ 'Primitive cardiomyocytes',
                       x == 'VCM' ~ 'Ventricular cardiomyocytes',
                       x == 'ESC' & study == 'Nagano' ~ 'ESC',
                       x == 'ESC' & study == 'Bonev2017' ~ 'ESC (serum)',
                       x == 'NPC' ~ 'Neural progenitors',
                       x == 'CN' ~ 'Cortical neurons',
                       T ~ x) %>%
           sub('^D', 'Day ', .) %>%
           rnm() %>%
           factor(samps),
         study = c(Nagano = 'Germline', Bonev2017 = 'Neural',
                   Stadhouders2018 = 'Reprogram',
                   Zhang2019 = 'Cardiac')[study] %>%
           factor(c('Germline', 'Neural', 'Reprogram', 'Cardiac'))) %>%
  na.omit() %>%
  filter(study == 'Germline') %>%
  ggplot(aes(x = lag, y = c, color = x)) + 
  geom_line(aes(color = x)) +
  coord_cartesian(xlim = c(bsz*2, 2e6)) +
  scale_color_manual(values = clrs[1:5]) +
  labs(x = 'Shift', y = 'Compartment score\nPearson auto-correlation') +
  facet_grid(.~'Compartment broadness') +
  scale_x_continuous(breaks = c(0,5e5, 1e6,1.5e6,2e6,3e6,4e6,5e6),
                     labels = c('0','500kb', '1mb', '1.5mb','2mb', '3mb', '4mb', '5mb')) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        legend.title = element_blank(),
        #legend.position = c(1,1),
        #legend.justification = c(1,1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'),
        axis.ticks.y = element_blank()) -> p
ggsave('sf3_d.pdf', p, height = 4.93, width = 6.5)

