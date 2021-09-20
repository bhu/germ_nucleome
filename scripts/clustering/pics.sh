#!/bin/bash

ncpus=20

ls -A1 *UMAP.sav | sed 's/.UMAP.sav//' | parallel -j $ncpus python plot_clust.py {}
mkdir -p pics/clu
mv  clu*corr*png pics/clu/
mv -t pics/ *marks.png *types.png
tar -zcf pics.tgz pics 
