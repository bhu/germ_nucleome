library(tidyverse)
library(GenomicRanges)
library(pals)
library(patchwork)

load('../data/peaks/CTCF/counts.rda')
m <- mcols(db)

samps <- c('ESC','EpiLC','d2PGCLC','d4c7PGCLC','GSC')
rnm <- function(x) sub('ESC', 'mESC', sub('PGC', '\nmPGC', x))
cm <- m %>% 
  as_tibble() %>%
  `+`(1) %>%
  log10() %>%
  as.list() %>%
  split(., sub('(.*)_(.*)_(.*)_.*', '\\1', names(.))) %>%
  lapply(function(x) {
    bind_cols(x) %>% rowMeans()
  }) %>%
  bind_cols() %>%
  select(all_of(samps))
samps <- rnm(samps)
colnames(cm) <- samps
clrs <- setNames(tableau20(9)[seq(1,9,2)], samps)

l <- c(0, 3.6)

ps <- lapply(seq_along(samps), function(i) {
  lapply(seq_along(samps), function(j) {
    p <- if (i > j) {
      tibble(y = cm[[i]], x = cm[[j]]) %>%
        ggplot(aes(x = x, y = y)) +
        geom_density_2d_filled(aes(color = ..level..), bins = 50) +
        scale_fill_viridis_d(option = 'A') +
        scale_color_viridis_d(option = 'A') +
        coord_cartesian(xlim = l, ylim = l, expand = F) +
        scale_x_continuous(breaks = 0:3) +
        scale_y_continuous(breaks = 0:3) +
        theme(plot.background = element_blank(),
              panel.background = element_rect(fill = 'black', color = 'black', size = 1),
              axis.text = element_blank(),
              legend.position = 'none',
              panel.grid = element_blank(),
              axis.title = element_blank())
    } else if (i == j) {
      tibble(x = cm[[i]]) %>%
        ggplot(aes(x = x)) +
        geom_histogram(fill = clrs[i], alpha = .5, bins = 50) +
        scale_y_continuous(expand = expansion(c(0, .05)), breaks = c(0, 2500, 5000)) +
        scale_x_continuous(expand = expansion(0), breaks = 0:3) +
        coord_cartesian(xlim = l, ylim = c(0, 6000)) +
        theme(plot.background = element_blank(),
              panel.background = element_rect(fill = NA, color = 'black', size = 1),
              axis.text = element_blank(),
              panel.grid = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.x = element_line(color = 'grey70', linetype = 'dashed'),
              axis.title = element_blank())
    } else {
      p <- ggplot() +
        annotate("text", x = 1, y = 1,
                 label = sprintf('Pearson\'s r\n%.2f', cor(cm[[i]], y = cm[[j]], method = 'pearson'))) + 
        theme(plot.background = element_blank(),
              panel.background = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank())
    }
    
    p <- if (i == 1 & j == length(samps)) {
      p + facet_grid(eval(parse(text = sprintf('"%s"', samps[i]))) ~
                       eval(parse(text = sprintf('"%s"', samps[j])))) +
        theme(strip.background.x = element_rect(fill = NA),
              strip.text.x = element_text(color = clrs[j], size = 13, face = 'bold'),
              strip.background.y = element_rect(fill = NA),
              strip.text.y = element_text(color = clrs[i], size = 13, face = 'bold'))
    } else if (i == 1) {
      p + facet_grid(.~eval(parse(text = sprintf('"%s"', samps[j])))) +
        theme(strip.background.x = element_rect(fill = NA),
              strip.text.x = element_text(color = clrs[j], size = 13, face = 'bold'))
    } else if (j == length(samps)) {
      p + facet_grid(eval(parse(text = sprintf('"%s"', samps[i]))) ~ .) +
        theme(strip.background.y = element_rect(fill = NA),
              strip.text.y = element_text(color = clrs[i], size = 13, face = 'bold'))
    } else {
      p
    }
    
    p <- if (i == length(samps) & j == 1) {
      p + theme(axis.text = element_text(color = 'black', size = 11))
    } else if (i == length(samps)) {
      p + theme(axis.text.x = element_text(color = 'black', size = 11))
    } else if (j == 1) {
      p + theme(axis.text.y = element_text(color = 'black', size = 11))
    } else {
      p
    }
    p
  })
})

{wrap_plots(unlist(ps, recursive = F), nrow = 5) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf3_f.pdf', ., height = 6.45, width = 6.8)

