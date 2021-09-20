# For DAPI analysis
# import library
suppressMessages(library(RImageJROI))
suppressMessages(library(rgeos))
suppressMessages(library(sp))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(tiff))
suppressMessages(library(ggsignif))
# set working directory
setwd("/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis")

# Firstly, extract the information of nuclear mask, nucleus, nuclear rim, and DAPI dense regions from czi file with .ijm script running on ImageJ
acceptable_cells <- list(list("ESC", c(1:18)), list("EpiLC", c(1:23)),list("d2PGCLC", c(1:19)), list("d4c7PGCLC", c(1:27)), list("GSC", c(1:22)), list("GSCLC", c(1:24)))
# Note that unsorted GSCs are excluded in the downsteam analysis.
# Then, calculate some stats of DAPI dense regions in the nucleus (i.e. distance from nucleus rim). Note that we have to upzip the output above (*_nuclear_periphery_roiset.zip) 
acceptable_cells %>%
  lapply(., function(p){
    #cell <- "EpiLC"
    #num <- "1"
    #dir <- "/Users/masahiro/Desktop/Experiments/IF/dapi/210821/EpiLC"
    cell <- p[[1]]
    num <- p[[2]]
    dir <- sprintf("/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/%s", cell)
    setwd(dir)
    num %>%
      lapply(., function(n){
        #n <- 1
        # read nuclear mask
        readTIFF(sprintf("%s_nucleus_mask.tif", n), all = T) -> a
        # read nuclear(rim)
        readTIFF(sprintf("%s_nucleus_rim_mask.tif", n), all = T) -> r
        # read image
        readTIFF(sprintf("%s.tif", n), all = T) -> im
        
        # make nucleus mask area data in each slice
        c(1:length(a)) %>% lapply(., function(l){
          a[[l]]%>% data.frame() %>% 
            mutate(row_ind = row_number()) %>%
            gather(-row_ind, key=col_ind, value=value) %>%
            mutate(col_ind=parse_number(col_ind)) %>% filter(value==1) %>% nrow(.) -> tmp
          c(l, tmp)
        }) %>% do.call(rbind, .) %>% `colnames<-`(c("slice", "area")) %>%
          data.frame() -> nucleus.area
        
        # decide +- 5 (i) slices from the max slice
        i <- 5
        nucleus.area %>% mutate(m=(area + lag(area)+ lead(area))/3) %>% na.omit() %>% {which.max(.$m)} -> max.slice
        usable.slices <- c((max.slice-i):(max.slice+i))
        
        # Next calculate DAPI dense regions and nuclear rim
        usable.slices %>%
          lapply(., function(l){
            # define nuclear rim area
            list.files(sprintf("%s_nuclear_periphery_roiset/", n), full.names = T) %>%
              # take the slice of the number
              .[str_detect(., sprintf("/%04d", l))] %>%
              # read the roi
              read.ijroi(.) %>%
              # change to sp object
              {SpatialPoints(.$coords)} -> boundary.sp
            
            # read information of DAPI dense regions
            fread(sprintf("%s_DAPI_dense_regions.csv", n)) %>% 
              dplyr::rename(id="V1") %>%
              # exclude too big regions
              filter(Area < 5000) %>%
              # consider pixel information
              transmute(id, col=round(X/0.0353844), row=round(Y/0.0353844), slice=sub(".*:z:([0-9]*)/.*", "\\1",Label), area=Area/(0.0353844^2)) -> dapi.dense
            
            # define dapi.dense regions position
            dapi.dense %>%
              filter(slice==l) %>%
              transmute(x=col, y=row) %>% SpatialPoints(.) -> dense.sp
            
            # measure distance in each pair of position and take the minimum and also measure the distance from nearly centroid position to exclude the points outside of the cell
            gDistance(dense.sp, boundary.sp , byid=T) %>%
              data.frame() %>%
              lapply(., function(x){
                data.frame(x) %>% mutate(row=row_number()) %>%
                  # choose the nearest boundary roi
                  .[.$x==min(.$x),] %>%
                  # if several mminimum points exit, choose one
                  .[1,]
              }) -> tmp
            
            c(1:length(tmp)) %>%
              lapply(., function(j){
                boundary.sp[tmp[[j]]$row,] -> nearest.bdr
                g <- c(x=mean(boundary.sp$x %>% range()),y=mean(boundary.sp$y %>% range())) %>% data.frame() %>% t() %>% SpatialPoints(.)
                gDistance(g, nearest.bdr , byid=T) %>% as.numeric() -> g_bdr
                gDistance(g, dense.sp[j] , byid=T) %>% as.numeric() -> g_p
                tmp[[j]] %>% `colnames<-`(c("dist", "bdr_row")) %>%
                  transmute(dist, gTobdr=g_bdr, gTopoint=g_p)
              }) %>%
              do.call(rbind, .) %>%
              bind_cols(dapi.dense %>% filter(slice==l),.) 
          }) %>%
          bind_rows() %>% fwrite(sprintf("%s_%s_DAPI_dense_stats.tsv",cell, n), sep="\t", col.names = T, row.names = F, scipen = 999)
      })
    
  })



# summarize
cells <- c("ESC", 
           "EpiLC",
           "d2PGCLC", 
           "d4c7PGCLC", 
           "GSC",
           "GSCLC")

cells %>% lapply(., function(c){
  list.files(sprintf("/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/%s", c), pattern = "DAPI_dense_stats.tsv", full.names = T) %>%
    lapply(., function(x){
      id2 <- x %>% sub(".*/.*_([0-9]*)_DAPI_dense_stats.tsv", "\\1", .) %>% parse_number(x)
      cell_n <- sub(".*/(.*)_[0-9].*", "\\1",x)
      fread(x) %>% mutate(number=id2, cell=cell_n)
    }) %>% bind_rows()
}) %>% bind_rows() %>%
  mutate(cell=factor(cell, levels=cells)) %>%
  filter(gTobdr > gTopoint) -> dense.dat

acceptable_cells %>%
  lapply(., function(p){
    dense.dat %>% filter(cell==p[[1]] & number %in% p[[2]])
  })  %>% bind_rows() -> dense.dat

save(file = "/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/dense.dat.rda", dense.dat)
####################################################################################################################################################################################################################################################################
## For DAPI variance analysis in the nucleus, we calculate DAPI normalized variance
acceptable_cells %>% lapply(., function(p){
  p[[1]] -> cell
  p[[2]] -> n
  setwd(sprintf("/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/%s", cell))
  n %>%
    lapply(., function(n){
      # read nuclear mask
      readTIFF(sprintf("%s_nucleus_mask.tif", n), all = T) -> a
      # read nuclear(rim)
      readTIFF(sprintf("%s_nucleus_rim_mask.tif", n), all = T) -> r
      # read image
      readTIFF(sprintf("%s.tif", n), all = T) -> im
      
      # make nucleus mask area data in each slice
      c(1:length(a)) %>% lapply(., function(l){
        a[[l]]%>% data.frame() %>% 
          mutate(row_ind = row_number()) %>%
          gather(-row_ind, key=col_ind, value=value) %>%
          mutate(col_ind=parse_number(col_ind)) %>% filter(value==1) %>% nrow(.) -> tmp
        c(l, tmp)
      }) %>% do.call(rbind, .) %>% `colnames<-`(c("slice", "area")) %>%
        data.frame() -> nucleus.area
      
      # decide the max slice
      nucleus.area %>% mutate(m=(area + lag(area)+ lead(area))/3) %>% na.omit() %>% {which.max(.$m)} -> l
      # max slice +- 2 slices
      
      ls <- c((l-2):(l+2))
      ls %>% lapply(., function(l){
        a[[l]] %>% {which(. > 0, arr.ind = T)} %>%
          data.frame() %>% 
          {mapply(function(row, col){
            im[[l]] %>% .[row, col] -> val
            c(row, col, val)  
          }, .$row, .$col, SIMPLIFY = F)} %>% do.call(rbind, .) %>% 
          data.frame() %>%
          `colnames<-`(c("row", "col", "sig")) %>% transmute(sig, number=n, celltype=cell) -> tmp
        r <- (max(tmp$sig) - min(tmp$sig))
        tmp %>% mutate(norm_range=sig/r, norm_mean=sig/mean(tmp$sig), norm_median=sig/median(tmp$sig), slice=l) %>%
          group_by(celltype, number, slice) %>%  summarise(Var1=var(norm_range), Var2=var(norm_mean), Var3=var(norm_median)) %>%
          ungroup() 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) -> dapi.varinace_all
save(file = "/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/DAPI_variance_all.rda", dapi.varinace_all)

##############################################################################################################################
##############################################################################################################################
# With Robject generated above, the folllowing commands are used for figure generation.
samps <- c("ESC", "EpiLC", "d2PGCLC", "d4PGCLC", "d4c7PGCLC", "GSC", "GSCLC")
clrs <- setNames(pals::tableau20(20)[c(1, 3, 5, 6, 7, 9, 11)], samps)

# cells used for analysis (We used sorted GSC and GSCLC for consistency.)
cells <- c("ESC", "EpiLC","d2PGCLC", "d4c7PGCLC", "GSC", "GSCLC")

setwd("/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis")
load("dense.dat.rda")

dense.dat %>% filter(area <3000 & area > 500) %>% 
  mutate(cell= factor(cell, levels=c("ESC", "EpiLC", "d2PGCLC", "d4c7PGCLC", "GSC", "GSCLC"))) -> dense.dat

dense.dat %>%
  filter(cell!="GSCLC") %>%
  ggplot(aes(x=cell, y=dist*0.035, fill=cell)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("ESC", "EpiLC"), c("EpiLC", "d2PGCLC"), c("d2PGCLC", "d4c7PGCLC"), c("d4c7PGCLC", "GSC")), 
                        test = wilcox.test,
                        step_increase = 0.075,
                        map_signif_level = TRUE) +
  theme_bw() +
  scale_fill_manual(values=clrs) +
  labs(x="", y="Distance of DAPI dense regions from nuclear periphery (um)") +
  theme(legend.position = "none") -> p
p
ggsave("dapi_dense_regions_distance_ESC_GSC.png", width = 5, height = 6)


dense.dat %>%
  filter(cell %in% c("d4c7PGCLC", "GSC", "GSCLC")) %>%
  ggplot(aes(x=cell, y=dist*0.035, fill=cell)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("d4c7PGCLC", "GSC"), c("GSC", "GSCLC"), c("d4c7PGCLC", "GSCLC")), 
                        test = wilcox.test,
                        step_increase = 0.075,
                        map_signif_level = TRUE) +
  theme_bw() +
  scale_fill_manual(values=clrs) +
  labs(x="", y="Distance of DAPI dense regions from nuclear periphery (um)") +
  theme(legend.position = "none") -> p
p
ggsave("dapi_dense_regions_distance_d4c7_GSCLC.png", width = 5, height = 6)

dense.dat %>% filter(area <3000 & area > 500) %>% 
  filter(cell!="GSCLC") %>%
  ggplot(aes(x=cell, y=area*(0.035^2), fill=cell)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("ESC", "EpiLC"), c("EpiLC", "d2PGCLC"), c("d2PGCLC", "d4c7PGCLC"), c("d4c7PGCLC", "GSC")), 
                        test = wilcox.test,
                        step_increase = 0.075,
                        map_signif_level = TRUE) +
  theme_bw() +
  scale_fill_manual(values=clrs) +
  labs(x="", y="Size of DAPI dense regions (um^2)") +
  theme(legend.position = "none") -> p
p
ggsave("dapi_dense_regions_size_ESC_GSC.png", width = 5, height = 6)

dense.dat %>% filter(area <3000 & area > 500) %>% 
  filter(cell %in% c("d4c7PGCLC", "GSC", "GSCLC")) %>%
  ggplot(aes(x=cell, y=area*(0.035^2), fill=cell)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("d4c7PGCLC", "GSC"), c("GSC", "GSCLC"), c("d4c7PGCLC", "GSCLC")),
                        test = wilcox.test,
                        step_increase = 0.075,
                        map_signif_level = TRUE) +
  theme_bw() +
  scale_fill_manual(values=clrs) +
  labs(x="", y="Size of DAPI dense regions (um^2)") +
  theme(legend.position = "none") -> p
p
ggsave("dapi_dense_regions_size_d4c7_GSCLC.png", width = 5, height = 6)


dense.dat %>% filter(area <3000 & area > 500) %>% 
  filter(cell %in% c("d4c7PGCLC", "GSC", "GSCLC")) %>%
  group_by(number, cell, slice) %>%
  summarise(total=n()) %>%
  ggplot(aes(x=cell, y=total, fill=cell)) +
  geom_boxplot() +
  ggsignif::geom_signif(comparisons = list(c("d4c7PGCLC", "GSC"), c("GSC", "GSCLC"), c("d4c7PGCLC", "GSCLC")), 
                        test = wilcox.test,
                        step_increase = 0.075,
                        map_signif_level = TRUE) +
  theme_bw() +
  scale_fill_manual(values=clrs) +
  labs(x="", y="Number of DAPI dense regions per slice")+
  theme(legend.position = "none") -> p
p
ggsave("dapi_dense_regions_number_d4c7_GSCLC.png", width = 5, height = 6)

## Analyze variance in the nucleus
load("DAPI_variance_all.rda")
cells <- c("ESC", "EpiLC","d2PGCLC", "d4c7PGCLC", "GSC_sort", "GSCLC_sort")

dapi.varinace_all %>%
  mutate(celltype= factor(celltype, levels=c("ESC", "EpiLC","d2PGCLC", "d4c7PGCLC", "GSC", "GSCLC"))) -> dapi.varinace_all

dapi.varinace_all %>%
  gather(-c("celltype", "number", "slice"), key=type, value=val) %>%
  mutate(celltype=factor(celltype, levels = c("ESC", "EpiLC","d2PGCLC", "d4c7PGCLC", "GSC", "GSCLC"))) -> var.dat

var.dat %>%
  filter(type=="Var2") %>%
  filter(celltype!="GSCLC") %>%
  ggplot(aes(x=celltype, y= val, fill=celltype)) +
  geom_boxplot(alpha=.4) +
  geom_violin(alpha=.4) +
  ggsignif::geom_signif(comparisons = list(c("ESC", "EpiLC"), c("EpiLC", "d2PGCLC"), c("d2PGCLC", "d4c7PGCLC"), c("d4c7PGCLC", "GSC")), 
                        test = wilcox.test,
                        step_increase = 0.075,
                        map_signif_level = TRUE) +
  theme_bw() +
  scale_fill_manual(values=clrs) +
  labs(x="", y="DAPI Variance within cell") +
  theme(legend.position = "none") -> p
p
ggsave("DAPI_uniformity_ESC_GSC.png", width=5, height=6, p)
