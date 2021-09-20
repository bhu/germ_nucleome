1. `umap_hdbscan.py` is ran using various parameters
2. `pics.sh` is ran to produce a set of diagnostic plots used to identify the optimal settings
3. `sep.py` is ran to split the input table
4. Individual subsets of the original input table is symlinked into new subdirectories (e.g., `0/dat.csv ` is linked to `./dat_0.csv`)
5. Go back to step 1 in each of the subdirectories and repeat until data points no longer form epigenetically distinct clusters
6. Finally, from the top directory, `merge.R` is ran to aggregate results

