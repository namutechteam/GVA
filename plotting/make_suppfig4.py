"""Supplementary Fig. 4 | Measurement groups and shot counts (redrawn with matplotlib from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: fewest measurement groups over the three colouring heuristics for QWC and GC partitions.
b: shots per energy evaluation (sigma = 1.6 mHa) without grouping and for the cheapest QWC and GC partition;
   labels give the GC/QWC ratio (sheet SuppFig4_Measurement).
Usage: python make_suppfig4.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import numpy as np
import openpyxl

from gva_style import COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
HEUR = ('largest_first', 'dsatur', 'sorted_insertion')

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
rows = list(wb['SuppFig4_Measurement'].iter_rows(values_only=True))
meas = [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


def best(r, kind, key):
    """Fewest groups or fewest shots over the colouring heuristics."""
    return min(r[f'{kind} {h} {key}'] for h in HEUR)


# grouping strategies are not ansaetze: greys from light (no grouping) to dark (GC)
NONE_, QWC, GC = '#CFCFCF', '#8C8C8C', '#333333'
labels = ['H$_2$', 'H$_4$', 'LiH', 'BeH$_2$', 'Cisplatin\nCAS(10e,10o)']
xx = np.arange(len(meas))

fig, axs = plt.subplots(1, 2, figsize=(COL2 * 0.82, 60 / 25.4), gridspec_kw={'wspace': 0.30})
fig.subplots_adjust(left=0.09, right=0.98, top=0.90, bottom=0.20)

ax = axs[0]
w = 0.38
q = [best(r, 'QWC', 'groups') for r in meas]
g = [best(r, 'GC', 'groups') for r in meas]
ax.bar(xx - w / 2, q, w, color=QWC, label='Qubit-wise commuting', lw=0)
ax.bar(xx + w / 2, g, w, color=GC, label='General commuting', lw=0)
for i in range(len(xx)):
    ax.text(xx[i] - w / 2, q[i] * 1.12, str(q[i]), ha='center', va='bottom', fontsize=5.6)
    ax.text(xx[i] + w / 2, g[i] * 1.12, str(g[i]), ha='center', va='bottom', fontsize=5.6)
ax.set_yscale('log')
ax.set_ylim(1, 1e4)
ax.set_ylabel('Measurement groups')
ax.legend(loc='upper left', handlelength=1.0, handleheight=0.8)

ax = axs[1]
u = [r['Shots, no grouping'] for r in meas]
qs = [best(r, 'QWC', 'shots') for r in meas]
gs_ = [best(r, 'GC', 'shots') for r in meas]
w = 0.27
ax.bar(xx - w, u, w, color=NONE_, label='No grouping', lw=0)
ax.bar(xx, qs, w, color=QWC, label='QWC', lw=0)
ax.bar(xx + w, gs_, w, color=GC, label='GC', lw=0)
for i in range(len(xx)):
    ax.text(xx[i] + w, max(qs[i], gs_[i]) * 1.25,   # clear of the neighbouring QWC bar
            f'{gs_[i] / qs[i]:.2f}×', ha='center', va='bottom', fontsize=5.6)
ax.set_yscale('log')
ax.set_ylabel('Shots per energy evaluation\n(σ = 1.6 mHa)')
ax.legend(loc='upper left', handlelength=1.0, handleheight=0.8)

for a in axs:
    a.set_xticks(xx)
    a.set_xticklabels(labels, linespacing=1.1)
    a.tick_params(axis='x', length=0, pad=3)
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=20))
    a.tick_params(axis='y', which='minor', labelleft=False)

fig.canvas.draw()
for a, lab in zip(axs, 'ab'):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'SuppFig4')
print('groups QWC', q, 'GC', g)
print('shots none', [f'{v:.3g}' for v in u], 'QWC', [f'{v:.3g}' for v in qs], 'GC', [f'{v:.3g}' for v in gs_],
      'ratio', [round(b / a, 2) for a, b in zip(qs, gs_)])
