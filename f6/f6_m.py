import sys
sys.path.append('../scripts')

import diverging_map

import pandas as pd
import numpy as np

import fanc
import fanc.plotting as fancplot

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import LogLocator
from matplotlib.patches import Arc, Circle, Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

RGB1 = np.array([73,152,148])
RGB2 = np.array([148, 103, 189])
nc = 99
colormap = diverging_map.ColorMapCreator(RGB1, RGB2, numColors=nc) 
clrs = colors.ListedColormap(np.c_[colormap.generateColorMap(RGB1, RGB2, 255), np.ones(nc)])

a = 112652e3
b = 112927e3
c = 112674e3
x1 = 112625483
x2 = 112897580


hic = fanc.load('../data/hicsr/GSC.GSCLC.10kb.fanc')
reg = 'chr13:112450000-113400000'
fig, axes = plt.subplots(5, 1, sharex=True, figsize=[3.4,4], constrained_layout=True,
                        gridspec_kw={'height_ratios': [4, 1,1,1,1]})
cl = 1.2
p1 = fancplot.TriangularMatrixPlot(hic, ax=axes[0], show_colorbar=False,
                                   vmin=-cl, vmax=cl, oe=False, log=True, 
                                   colormap=clrs, max_dist='950kb')
p1.plot(reg)

axes[0].add_patch(Circle(((x1 + x2) / 2, (.5)**.5 * (x2 - x1) * .8), 30000,
                             fill=False, lw=2, alpha=.5, edgecolor='purple'))
axes[0].add_patch(Rectangle((a, 0), 30000, 2.75e5, 315,
                            fill=False, lw=2, alpha=.25, edgecolor='blue'))

div = make_axes_locatable(axes[0])
cax = div.append_axes("right", size="2%", pad="1%")
p1.collection.set_clim(vmin=-cl, vmax=cl)
cbar = plt.colorbar(p1.collection, ax=axes[0], cax=cax, extend='both',
                    aspect=10)
cbar.ax.get_yaxis().labelpad = 10
cbar.ax.set_ylabel('log(GSC/GSCLC)', rotation=270)


p2 = fancplot.LinePlot('../data/tracks/Ddx4/K27ac.log2.bw', ax=axes[1], colors=('#6AA56E'),
                       n_yticks=2, fill=True, bin_size=500, ylim=[-.8,.8])

p3 = fancplot.LinePlot('../data/tracks/Ddx4/CTCF.log2.bw', ax=axes[2], colors=('#5372AB'), 
                       n_yticks=2, fill=True, bin_size=500, ylim=[-1,1])

p4 = fancplot.LinePlot('../data/tracks/Ddx4/K36me3.log2.bw', ax=axes[3], colors=('#B65655'),
                       n_yticks=2, fill=True, bin_size=2000, ylim=[-.301,.301])
                       
for p, ax in zip([p2, p3, p4], axes[1:4]):
    p.plot(reg)
    p.lines[0].set_linewidth(1)
    p.fills[0].set_alpha(.5)
    ax.axvline(a, color='#00A1FF', alpha=.2, lw=3)
    ax.axvline(b, color='#00A1FF', alpha=.2, lw=3)
    ax.axvline(c, color='red', alpha=.2, lw=3)
    
p5 = fancplot.GenePlot('../data/resources/gencode.vM25.basic.protein_coding.annotation.gff3.gz',
                       ax=axes[4], group_by="gene_id", squash=True,
                       label_field='gene_name', arrow_size=4, line_width=2, box_height=0.5)
p5.plot(reg)

axes[4].set_xticks([112.5e6, 113e6])

for t in p5.texts:
    if t.get_text() == 'Ddx4':
        t.set_color('darkturquoise')
    else:
        t.set_visible(False)

for ax, lab in zip(axes, ['Hi-C', 'H3K27ac', 'CTCF', 'H3K36me3', '']):
    if lab != 'Hi-C':
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="2%", pad="1%")
        fig.delaxes(cax)
    if lab != '':
        at = AnchoredText(lab, loc='upper right', frameon=False, pad=0, borderpad=0)
        ax.add_artist(at)

plt.savefig('f6_m.pdf', transparent=True, dpi=600, bbox_inches='tight')

