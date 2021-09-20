library(tidyverse)
library(mixtools)
library(patchwork)
library(pals)

dat <- readRDS('../data/peaks/ATAC/fetal_comparison.rds') %>%
  select(-c(chr, start, end, GSCLC))

samps <- c('ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC', 'd4PGCLC') 
clrs <- setNames(c(tableau20(9)[seq(1,9,2)], '#81642A'), samps)
clrs2 <- setdiff(names(dat), samps) %>% 
  {setNames(rep_along(., 'black'), .)}

fits <- dat %>%
  lapply(function(x) {
    normalmixEM2comp(x[x > 2], mu = c(4, 8), sigsqrd = c(.5, .5), lambda = .5)
  })

s <- 'ESC'
h <- hist(dat[[s]], 100, plot = F)
w <- median(diff(h$breaks))
h <- tibble(x = h$mids, y = h$counts) %>%
  mutate(grp = case_when(x > fits[[s]]$mu[2] ~ 'More open',
                         x < fits[[s]]$mu[1] ~ 'Less open',
                         T ~ 'Intermediate'))
d <- lapply(1:2, function(i) {
  tibble(x = h$x, 
         y = fits[[s]]$lambda[i] * dnorm(h$x, mean = fits[[s]]$mu[i], sd = fits[[s]]$sigma[i]),
         cmp = i)
}) %>% 
  bind_rows() %>%
  mutate(cmp = c('Less open', 'More open')[cmp] %>% factor(c('More open', 'Less open')))

m <- 2e4
b <- 0

eclrs <- c('More open' = 'firebrick3',
           'Less open' = 'steelblue2',
           'Intermediate' = 'grey50')
p1 <-  ggplot(mapping = aes(x = x, y = y)) +
  geom_col(aes(fill = grp), data = h, width = w, alpha = .6, show.legend = F) +
  geom_line(aes(y = m * y + b, color = cmp), data = d, size = 2, alpha = .5) +
  scale_color_manual(values = eclrs) +
  scale_fill_manual(values = eclrs) +
  scale_y_continuous(expand = expansion(c(0, .05)),
                     sec.axis = sec_axis(~ (. - b)/m, name = 'Fitted density')) +
  scale_x_continuous(expand = expansion(0)) +
  labs(x = 'log2(FPKM + 1)', y = '# of sites') +
  facet_grid(.~ 'Example fit for mESC') +
  theme(legend.position = 0:1,
        legend.justification = 0:1,
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_rect(fill = '#ffffff99', color = NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'), 
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(color = 'black', size = 11),
        axis.title.y.right = element_text(vjust = 1.5),
        panel.grid.major.y = element_line(color = 'grey70', linetype = 'dashed'))
  
p2 <- fits %>%
  Map(function(x, samp) {
    dat[[samp]] %>%
      {tibble(low = sum(. < x$mu[1]),
              high = sum(. > x$mu[2]))}
  }, ., names(.)) %>%
  bind_rows(.id = 'samp') %>%
  mutate(r = high / low) %>%
  arrange(r) %>%
  mutate(clr = c(clrs, clrs2)[samp],
         samp = case_when(samp == 'facialprominence' ~ 'craniofacial prominence',
                          samp == 'neuraltube' ~ 'neural tube',
                          T ~ samp) %>%
           sub('ESC', 'mESC', .) %>%
           sub('PGC', ' mPGC', .) %>%
           fct_inorder()) %>%
  ggplot(aes(x = samp, y = r)) +
  geom_hline(yintercept = 1) +
  geom_linerange(aes(ymin = r, ymax = 1, color = clr)) +
  geom_point(aes(color = clr)) +
  scale_color_identity() +
  ylab('# of more open sites / less') +
  facet_grid(.~ 'Comparison to mouse fetal samples') +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(color = 'black', size = 13, face = 'bold'), 
        axis.text = element_text(color = 'black', size = 11),
        legend.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        axis.ticks = element_blank(),
        axis.line.x = element_line(color = 'black'),
        strip.background = element_rect(fill = NA))

{ wrap_plots(p1, p2, widths = c(1,1.2)) &
  theme(plot.background = element_blank()) } %>%
  ggsave('sf2_e.pdf', ., height = 5, width = 9) 
