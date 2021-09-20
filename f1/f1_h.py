import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools.lib.plotting
from matplotlib.colors import LogNorm
from matplotlib.ticker import EngFormatter
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
samps = ['ESC', 'EpiLC','d2PGCLC', 'd4c7PGCLC', 'GSC']
clrs = {'ESC': 'tab:blue', 'EpiLC': 'tab:orange', 'd2PGCLC': 'tab:green',
       'd4c7PGCLC': 'tab:red', 'GSC': 'tab:purple'}
res = 25000

region = ('chr3', 5_500_000, 12_000_000)
norm = LogNorm(vmax=0.1, vmin=5e-4)
fig, axes = plt.subplots(nrows=5, ncols=1, constrained_layout=True, figsize=(3, 5.5))
sclr = '#7AFAC2'
tclr = '#52A883'
for ax, s in zip(axes.flat, samps):
    clr = cooler.Cooler('../data/mcool/%s.mcool::/resolutions/%d' % (s, res))
    data = clr.matrix(balance=True).fetch(region)
    im = pcolormesh_45deg(ax, data, start=region[1], resolution=res, norm=norm, cmap='fall')
    ax.set_ylim(0, 2800000)
    format_ticks(ax, rotate=False)
    ax.set_xticks(np.arange(7,12) * 1e6)
    ax.set_xlim(region[1]+1e6, region[2]-1e6)
    ax.xaxis.tick_top()
    ax.set_yticks([0, 1e6, 2e6])
    if s != 'ESC':
        ax.xaxis.set_visible(False)
        ax.set_yticks([])
    
    div = make_axes_locatable(ax)
    cax = div.append_axes("left", size="2%", pad="3%")
    cbar = plt.colorbar(im, ax=ax, cax=cax, extend='max', aspect=5)
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.ax.get_yaxis().labelpad = -50
    cbar.ax.set_ylabel("")
    if s != 'GSC':
        fig.delaxes(cax)
plt.savefig('f1_h.pdf', transparent=True, dpi=300, bbox_inches='tight')

