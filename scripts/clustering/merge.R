library(data.table)
library(tidyverse)

d <- system("find . -name 'dat_*.csv'", intern=T) %>%
  tibble(f = .) %>% 
  mutate(s = sub('^\\./' , '', f) %>% 
           sub('.csv', '', .) %>% 
           sub('dat_', '', .) %>%
           gsub('/', '_', .))

i <- str_count(d$s, '_') %>% 
  max() %>% {1:.} %>% 
  lapply(function(x) {
    sub(sprintf("^((?:[^_]*_){%d}).*", x), "\\1", d$s) %>%
      sub('_$', '', .) %>%
      {.[duplicated(.)]} %>% 
      unique()
  }) %>%
  unlist()

d %>% 
  filter(!(s %in% i)) %>%
  split(.,.$s) %>%
  lapply(function(x) {  
    fread(x$f, select = c('type', 'crd')) 
  }) %>% 
  rbindlist(idcol = 'clust')  %>%
  distinct(clust, crd, type) %>% 
  merge(fread('dat.csv'), ., all.x = T) %>%
  mutate(cl = {as.numeric(factor(clust))-1} %>% 
           replace_na(-1)) %>%
  fwrite('lab_dat.csv', row.names = F)

