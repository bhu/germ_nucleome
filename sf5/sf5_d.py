import pandas as pd
import numpy as np
import fanc
import fanc.plotting as fancplot
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['text.color'] = 'black'
rcParams['axes.labelcolor'] = 'black'
rcParams['xtick.color'] = 'black'
rcParams['ytick.color'] = 'black'
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import LogLocator
from matplotlib.patches import Arc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
abc = pd.read_csv('../data/abc/ESC.csv.gz')
hic = fanc.load('../data/hicsr/ESC.10kb.fanc')
ep = abc.loc[(abc['class'] != 'promoter') & 
             (abc['ABC.Score'] > 0.02) & 
             (abc['TargetGene'] == 'Nanog') &
             (abc['start'] > 122.3e6) &
             (abc['end'] < 123e6)]\
        .sort_values('ABC.Score', ascending=False)
reg = 'chr6:122.3mb-123mb'
fig, axes = plt.subplots(6, 1, sharex=True, figsize=[3.5,4], constrained_layout=True,
                        gridspec_kw={'height_ratios': [3,1,1,1,1,1]})
cl = 3
p1 = fancplot.TriangularMatrixPlot(hic, ax=axes[0], show_colorbar=False,
                                   vmin=-cl, vmax=cl, oe=True, log=True, 
                                   colormap='twilight_shifted', max_dist='700kb')
p1.plot(reg)

div = make_axes_locatable(axes[0])
cax = div.append_axes("right", size="2%", pad="1%")
p1.collection.set_clim(vmin=-cl, vmax=cl)
cbar = plt.colorbar(p1.collection, ax=axes[0], cax=cax, extend='both',
                    aspect=10)
cbar.ax.get_yaxis().labelpad = 5
cbar.ax.set_ylabel('log(O/E)', rotation=270)

for index, row in ep.iterrows():
    g = row['TargetGeneTSS']
    e = (row['start'] + row['end']) / 2
    a = Arc(((g + e) / 2, 0),  abs(g - e), 1, 
            theta1=0, theta2=180, edgecolor='#1f77b4', lw=1)
    axes[1].add_patch(a)

axes[1].set_ylim([0,.52])
axes[1].invert_yaxis()
axes[1].set_axis_off()

p2 = fancplot.HicSlicePlot(hic, '%s:%d' % (row['chr'], row['TargetGeneTSS']), ax=axes[2],
                           colors=('#ff7f0e'), n_yticks=2, fill=False)
p3 = fancplot.LinePlot('../data/tracks/sf2/ESC_ATAC.bw', ax=axes[3],
                       colors=('#9467bd'), n_yticks=2, fill=True, bin_size=200)
p4 = fancplot.LinePlot('../data/tracks/sf2/ESC_K27ac.bw', ax=axes[4],
                       colors=('#6AA56E'), n_yticks=2, fill=True, bin_size=200)
                       
for p in [p2, p3, p4]:
    p.plot(reg)
    p.lines[0].set_linewidth(1)
    if isinstance(p, fanc.plotting.hic_plotter.HicSlicePlot):
        continue
    else:
        p.fills[0].set_alpha(.5)
p5 = fancplot.GenePlot('../data/resources/gencode.vM25.basic.protein_coding.annotation.gff3.gz',
                       ax=axes[5], group_by="gene_id", squash=True,
                       label_field='gene_name')
p5.plot(reg)
axes[5].set_xticks([122.3e6, 122.6e6, 122.9e6])
for t in p5.texts:
    if t.get_text() == 'Nanog':
        t.set_color('orangered')
    else:
        t.set_visible(False)
for ax, lab in zip(axes, ['Hi-C', 'ABC', 'Virtual 4-C', 'ATAC-seq', 'H3K27ac', '']):
    if lab != 'Hi-C':
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="2%", pad="1%")
        fig.delaxes(cax)
    if lab != '':
        at = AnchoredText(lab, loc='upper right', frameon=False, pad=0)
        ax.add_artist(at)
plt.savefig('sf4_f.pdf', transparent=True, dpi=600, bbox_inches='tight')

