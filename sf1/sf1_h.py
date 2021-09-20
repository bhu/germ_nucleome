from pymol import cmd
for i, s in enumerate(['ESC', 'EpiLC', 'd2PGCLC', 'd4c7PGCLC', 'GSC']):
    cmd.reinitialize()
    cmd.load('../data/3dmods/%d_16.xyz' % (i + 1))
    cmd.spectrum()
    cmd.orient()
    cmd.png('sf1_h.%s.png' % s, ray=1)

import matplotlib.pyplot as plt
import matplotlib as mpl
fig, ax = plt.subplots(figsize=(.7, .5))
fig.subplots_adjust(bottom=0.5)
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                  cax=ax, orientation='horizontal')
cb.set_ticks([0, 1])
cb.set_ticklabels(['pter', 'qter'])
cb.outline.set_visible(False)
cb.ax.tick_params(width=0)
plt.savefig('sf1_h.leg.pdf', bbox_inches='tight')

