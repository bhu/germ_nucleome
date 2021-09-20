import umap
import umap.plot
import pandas as pd
import joblib
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_scatter_density
import hdbscan

d = pd.read_csv('dat.csv', index_col=['type', 'crd'])

um = umap.UMAP(verbose=True,
               n_epochs = 200,
               n_neighbors=int(sys.argv[1]),
               min_dist=float(sys.argv[2]),
               spread=1.0,
               metric=sys.argv[3],
               n_components=int(sys.argv[4]),
               random_state=42).fit(d.values)

fnm =  '%s_%s_%s_%s' % (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
joblib.dump(um, fnm + '.UMAP.sav')



if sys.argv[4] == '2':
  umap.plot.points(um, labels=d.index.get_level_values('type').values)
  plt.savefig(fnm + '.types.png', bbox_inches="tight")
  plt.close()
  fig, axs = plt.subplots(4, 4, figsize=(18, 10), subplot_kw=dict(projection='scatter_density'))
  faxs = axs.flatten()
  for i in range(d.shape[1]):
    mark = d.columns[i]
    c = d[mark].values
    dens = faxs[i].scatter_density(x=um.embedding_[:,0], y=um.embedding_[:,1],
                                   c=c, vmin=-2, vmax=2, cmap='coolwarm')
    faxs[i].set_facecolor('black')
    faxs[i].set_title(mark, loc='left', x=0.05, y=0.85, color='white')
    faxs[i].set_xticks([])
    faxs[i].set_yticks([])
    div = make_axes_locatable(faxs[i])
    cax = div.append_axes("bottom", size="5%", pad=0.05)
    fig.colorbar(dens, cax=cax, orientation="horizontal")
  plt.savefig(fnm + '.marks.png', bbox_inches="tight")

fnm =  '%s_%s_%s_%s' % (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
vs = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
for i in range(len(vs)):
    mem = joblib.Memory(location='cache/%d' % i)
    x = vs[i]
    for j in range(i, len(vs)):
        y = vs[j]
        hd = hdbscan.HDBSCAN(min_cluster_size=y,
                             min_samples=x,
                             core_dist_n_jobs=4,
                             memory=mem) \
                    .fit(um.embedding_)
        joblib.dump(hd, 'clu_%s_%d_%d.sav' % (fnm, y, x))
