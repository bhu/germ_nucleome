import sys
sys.path.append('../scripts')
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
import diverging_map
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
res = 250000
ms = []
chrom = 'chr1'
for s in ['GSC', 'GSCLC']:
    clr = cooler.Cooler('../data/mcool/%s.mcool::/resolutions/%d' % (s, res))
    pmap = adaptive_coarsegrain(clr.matrix(balance=True).fetch(chrom),
                                clr.matrix(balance=False).fetch(chrom),
                                cutoff=3, max_levels=8)
    oemap = oe(pmap)
    cmap = np.corrcoef(oemap)
    ms.append(interp_nan(pmap))
dmap = ms[1] / ms[0]
norm = LogNorm(vmin=1e-6, vmax=5e-1)
upmap = masked_array(ms[0], tri_mask(pmap, 'upper', 1))
lomap = masked_array(ms[1], tri_mask(cmap, 'lower', -1))
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True, figsize=(3, 1.5))
upi = ax.imshow(upmap, norm=norm, cmap='fall')
loi = ax.imshow(lomap, norm=norm, cmap='fall')
upi.set_rasterized(True)
loi.set_rasterized(True)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
div = make_axes_locatable(ax)
cax = div.append_axes("left", size="3%", pad="2%")
cbar = plt.colorbar(upi, ax=ax, cax=cax, extend='both', aspect=10)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.get_yaxis().labelpad = -50
cbar.ax.set_ylabel('Contact probability', rotation=90)
cbar.solids.set_rasterized(True)
ts = []
ts.append(AnchoredText('GSC', loc='upper center', frameon=True, pad=0.25))
ts.append(AnchoredText('GSCLC', loc='lower center', frameon=True, pad=0.25))
for t in ts:
    t.patch.set(alpha=0.7, color='white')
    ax.add_artist(t)
plt.savefig('f6d_bot.pdf', transparent=True, dpi=300)  
RGB2 = np.array([73,152,148])
RGB1 = np.array([148, 103, 189])
nc = 99
colormap = diverging_map.ColorMapCreator(RGB1, RGB2, numColors=nc) 
clrs = colors.ListedColormap(np.c_[colormap.generateColorMap(RGB1, RGB2, 255), np.ones(nc)])
norm2 = LogNorm(vmin=.25, vmax=4)
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True, figsize=(3, 1.1))
im = ax.imshow(dmap, norm=norm2, cmap=clrs)
im.set_rasterized(True)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
div = make_axes_locatable(ax)
cax = div.append_axes("right", size="5%", pad="4%")
cbar = plt.colorbar(im, ax=ax, cax=cax, extend='both', aspect=10)
cbar.ax.get_yaxis().labelpad = 10
#cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.set_ylabel('GSCLC/GSC', rotation=270)
cbar.solids.set_rasterized(True)
cbar.ax.yaxis.set_ticks([4,1,1/4])
cbar.ax.yaxis.set_ticklabels(['4','1','0.25'])
plt.savefig('f6d_top.pdf', transparent=True, dpi=300)  
