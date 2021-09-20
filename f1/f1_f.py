import matplotlib.pyplot as plt
from numpy.ma import masked_array
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.offsetbox import AnchoredText
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'sans-serif']
def tri_mask(A, side = 'upper', k = 1):
    m = A.shape[0]
    r = np.arange(m)
    if side == 'upper':
        mask = r[:,None] > (r - k)
    else:
        mask = r[:,None] < (r - k)
    return mask
samps = ['ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC' ]
fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=False, figsize=(15, 1.8))
norm = LogNorm(vmin = .5, vmax = 2)
for ax, s in zip(axes.flat, samps):
    scis = pd.read_csv('../data/saddle/' + s + '.cis.digitized.tsv', sep='\t')
    strans = pd.read_csv('../data/saddle/' + s + '.trans.digitized.tsv', sep='\t')
    dcis = np.load('../data/saddle/' + s + '.cis.saddledump.npz')
    dtrans = np.load('../data/saddle/' + s + '.trans.saddledump.npz')

    vcis = scis.groupby(["E1.d"])['E1'].median().values[1:-1]
    vtrans = strans.groupby(["E1.d"])['E1'].median().values[1:-1]
    xs = np.arange(49)

    upmap = dcis['saddledata'][1:-1,1:-1]
    lomap = dtrans['saddledata'][1:-1,1:-1]
    upmap = masked_array(upmap, tri_mask(upmap, 'upper', 0))
    lomap = masked_array(lomap, tri_mask(lomap, 'lower', 0))

    idx = [i for i in range(len(upmap)) if dcis['binedges'][i] < 0 and dcis['binedges'][i+1] > 0][0]
    ix = (idx + .5) / len(upmap)
    iy = 1 - ix
    upi = ax.imshow(upmap, norm=norm, cmap='RdYlBu_r')
    loi = ax.imshow(lomap, norm=norm, cmap='RdYlBu_r')
    ax.plot([1, 0], [0, 1], transform=ax.transAxes, color='white', linewidth=5)
    ax.axhline(idx, color='white', linewidth=5, alpha=.5)
    ax.axvline(idx, color='white', linewidth=5, alpha=.5)

    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    div = make_axes_locatable(ax)
    cax = div.append_axes("left", size="5%", pad="3%")
    cbar = plt.colorbar(loi, ax=ax, cax=cax, extend='both', aspect=10, ticks=[.5, 1, 2])
    cbar.ax.set_yticklabels(['.5', '1', '2'])
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.ax.minorticks_off()
    cbar.ax.get_yaxis().labelpad = -35
    cbar.ax.set_ylabel('Observed / expected', rotation=90)
    cbar.solids.set_rasterized(True)
    if s != 'ESC':
        fig.delaxes(cax)
    else:
        ts =[]
        ts.append(AnchoredText('Cis', loc='upper right', frameon=True, pad=0.25))
        ts.append(AnchoredText('Trans', loc='lower left', frameon=True, pad=0.25))
        ts.append(AnchoredText('B-B', loc='lower right', bbox_to_anchor=(ix, iy),
                               bbox_transform=ax.transAxes, frameon=True, pad=0.25))
        ts.append(AnchoredText('A-A', loc='upper left', bbox_to_anchor=(ix, iy),
                               bbox_transform=ax.transAxes, frameon=True, pad=0.25))
        ts.append(AnchoredText('B-A', loc='upper right', bbox_to_anchor=(ix, iy),
                               bbox_transform=ax.transAxes, frameon=True, pad=0.25))
        ts.append(AnchoredText('A-B', loc='lower left', bbox_to_anchor=(ix, iy),
                               bbox_transform=ax.transAxes, frameon=True, pad=0.25))
        for t in ts:
            t.patch.set(alpha=0.7, edgecolor=(1,1,1,0))
            ax.add_artist(t)
plt.savefig('f1_f.pdf', transparent=True, dpi=300)

