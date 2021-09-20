import matplotlib
import matplotlib.pyplot as plt
from numpy.ma import masked_array
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.offsetbox import AnchoredText
import numpy as np
import cooler
import cooltools.lib.plotting
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan, observed_over_expected
import pandas as pd
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'sans-serif']
def oe(A):
    A = np.array(A)
    A[~np.isfinite(A)] = 0
    mask = A.sum(axis=0) > 0
    OE, _, _, _ = observed_over_expected(A, mask)
    OE = np.clip(OE, 0, np.percentile(OE[mask, :][:, mask], 99.9))
    return OE
def tri_mask(A, side = 'upper', k = 1):
    m = A.shape[0]
    r = np.arange(m)
    if side == 'upper':
        mask = r[:,None] > (r - k)
    else:
        mask = r[:,None] < (r - k)
    return mask

samps = ['ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC']
res = 250000
chrom = 'chr1'
fig, axes = plt.subplots(nrows=1, ncols=5, constrained_layout=False, figsize=(15, 1.8))
norm = LogNorm(vmin=1e-6, vmax=5e-1)
for ax, s in zip(axes.flat, samps):
    clr = cooler.Cooler('../data/mcool/%s.mcool::/resolutions/%d' % (s, res))
    pmap = adaptive_coarsegrain(clr.matrix(balance=True).fetch(chrom),
                                clr.matrix(balance=False).fetch(chrom),
                                cutoff=3, max_levels=8)
    pmap = interp_nan(pmap)
    oemap = oe(pmap)
    cmap = np.corrcoef(oemap)
    upmap = masked_array(pmap, tri_mask(pmap, 'upper', 1))
    lomap = masked_array(cmap, tri_mask(cmap, 'lower', -1))
    upi = ax.imshow(upmap, norm=norm, cmap='fall')
    loi = ax.imshow(lomap, cmap = 'coolwarm', vmax=.7, vmin=-.7)
    upi.set_rasterized(True)
    loi.set_rasterized(True)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    div = make_axes_locatable(ax)
    locax = div.append_axes("left", size="5%", pad="3%")
    lobar = plt.colorbar(loi, ax=ax, cax=locax, extend='both', aspect=10)
    lobar.ax.yaxis.set_ticks_position('left')
    lobar.ax.get_yaxis().labelpad = -40
    lobar.ax.set_ylabel('Pearson\'s r', rotation=90)
    lobar.solids.set_rasterized(True)
    lobar.set_ticks([-.5,0,.5])
    upcax = div.append_axes("right", size="5%", pad="3%")
    upbar = plt.colorbar(upi, ax=ax, cax=upcax, extend='max', aspect=10)
    upbar.ax.get_yaxis().labelpad = 15
    upbar.set_ticks([1e-1, 1e-3, 1e-5])
    upbar.ax.set_ylabel('Contact probability', rotation=270)
    upbar.solids.set_rasterized(True)
    if s != 'ESC':
        fig.delaxes(locax)
    if s != 'GSC':
        fig.delaxes(upcax)
plt.savefig('f1_e.pdf', transparent=True, dpi=300)

