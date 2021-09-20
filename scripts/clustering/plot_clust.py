import sys
from textwrap import wrap
import umap
import umap.plot
import numpy as np
import hdbscan
from joblib import load

vs = [50, 100, 200, 500, 1000, 2000, 5000]
um = load('%s.UMAP.sav' % re.sub(r'_[0-9]*$', '_2', sys.argv[1]))
for i in range(len(vs)):
    x = vs[i]
    for j in range(i, len(vs)):
        y = vs[j]
        try:
            hd = load('clu_%s_%d_%d.sav' % (sys.argv[1], y, x))
            umap.plot.points(um, width=300, height=200, labels=hd.labels_)
            ax = plt.gca()
            uniq, cts = np.unique(hd.labels_, return_counts=True)
            ttl = 'minpts = %s, minsamps = %s ' % (y, x) \
                   + ', '.join(['%d:%d' % (a, b) for a, b in zip(uniq, cts)])
            ax.set_title("\n".join(wrap(ttl, 40)))
            ax.get_legend().remove()
            plt.savefig('clu_%s_%d_%d.png' % (sys.argv[1], y, x), dpi=300, bbox_inches="tight")
            plt.close()
        except:
            print('%s - %s is bad' % (y, x))

