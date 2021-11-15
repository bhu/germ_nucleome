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
a = 112652e3
b = 112927e3
c = 112674e3
x1 = 112625483
x2 = 112897580
for s in ['d4c7PGCLC', 'GSC']:
    hic = fanc.load('../data/hicsr/' + s + '.10kb.fanc')
    reg = 'chr13:112450000-113400000'
    fig, axes = plt.subplots(5, 1, sharex=True, figsize=[2.9,3.5], constrained_layout=True,
                            gridspec_kw={'height_ratios': [4, 1,1,1,1]})
    cl = 3
    p1 = fancplot.TriangularMatrixPlot(hic, ax=axes[0], show_colorbar=False,
                                       vmin=-cl, vmax=cl, oe=True, log=True, 
                                       colormap='twilight_shifted', max_dist='950kb')
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
    cbar.ax.get_yaxis().labelpad = 5
    cbar.ax.set_ylabel('log(O/E)', rotation=270)

    p2 = fancplot.LinePlot('../data/tracks/Ddx4/' + s + '_K27ac.bw', ax=axes[1], colors=('#6AA56E'),
                           n_yticks=2, fill=True, bin_size=500, ylim=[0,5.1])

    p3 = fancplot.LinePlot('../data/tracks/Ddx4/' + s + '_CTCF.bw', ax=axes[2], colors=('#5372AB'), 
                           n_yticks=2, fill=True, bin_size=500, ylim=[0,5.1])

    p4 = fancplot.LinePlot('../data/tracks/Ddx4/' + s + '_K36me3.bw', ax=axes[3], colors=('#B65655'),
                           n_yticks=2, fill=True, bin_size=2000, ylim=[0,1.02])

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

    plt.savefig('f3_f_' + s + '.pdf', transparent=True, dpi=600, bbox_inches='tight')
    plt.close()

for s in ['d4c7PGCLC','GSC']:
    reg = 'chr13:112927200-112928000'
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=[3.08,.65], constrained_layout=True,
                            gridspec_kw={'height_ratios': [1]})
    axes.axhline(0, c='black', lw=.5)
    axes.axvline(112927660, color='#00A1FF', alpha=.2, lw=30)
    p1 = fancplot.LinePlot(['../data/tracks/Ddx4/%s.CTSS.normalized.%s.bw' % (s, std) for std in ['plus', 'minus']],
                           ax=axes, colors=('#FF00FF'),n_yticks=2, fill=True, ylim=[-10,10])

    p1.plot(reg)
    
    p1.lines[0].set_linewidth(1)
    p1.lines[1].set_linewidth(1)

    axes.add_artist(AnchoredText(' +', loc='upper left', frameon=False, pad=0, borderpad=0)) 
    axes.add_artist(AnchoredText('NET-CAGE', loc='upper right', frameon=False, pad=0, borderpad=0)) 
    axes.add_artist(AnchoredText(' –', loc='lower left', frameon=False, pad=0, borderpad=0)) 

    div = make_axes_locatable(axes)
    cax = div.append_axes("right", size="2%", pad="1%")
    fig.delaxes(cax)

    axes.set_xticks([112927500, 112928000])
    axes.set_xticklabels(['chr13:112927500', '11928000'])
    plt.savefig(s + '_bidir.pdf', transparent=True, dpi=600, bbox_inches='tight')
    plt.close()

