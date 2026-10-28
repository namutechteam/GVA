"""Supplementary Fig. 1 | Spin contamination (redrawn with matplotlib from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
<S^2> of the final state against |dE| vs FCI for every single- and two-layer UCCSD and GVA run on the
small molecules (sheet Runs_SmallMolecules). |dE| below 1e-6 mHa and <S^2> below 1e-12 drawn at those floors.
Usage: python make_suppfig1.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator
import numpy as np
import openpyxl

from gva_style import C, COL1, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
FLOOR, S2_FLOOR, CHEM = 1e-6, 1e-12, 1.6

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


runs = [r for r in sheet('Runs_SmallMolecules') if r['Ansatz'] in ('UCCSD', 'GVA')]

fig, ax = plt.subplots(figsize=(COL1, 62 / 25.4))
for r in runs:
    col = C[r['Ansatz']]
    ax.plot(max(abs(r['dE vs FCI (mHa)']), FLOOR), max(r['<S^2>'], S2_FLOOR),
            'o' if r['Layers'] == 1 else 's', ms=3, mec=col,
            mfc='white' if r['Ansatz'] == 'UCCSD' else col, mew=0.7, alpha=0.85, zorder=3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.axvline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7, zorder=1)
ax.set_xlabel('|ΔE| relative to FCI (mHa)')
ax.set_ylabel(r'$\langle\mathrm{S}^2\rangle$ of final state')
handles = [Line2D([], [], marker='o', ls='', ms=3.4, mec=C['UCCSD'], mfc='white', mew=0.7, label='UCCSD 1L'),
           Line2D([], [], marker='s', ls='', ms=3.4, mec=C['UCCSD'], mfc='white', mew=0.7, label='UCCSD 2L'),
           Line2D([], [], marker='o', ls='', ms=3.4, color=C['GVA'], label='GVA 1L'),
           Line2D([], [], marker='s', ls='', ms=3.4, color=C['GVA'], label='GVA 2L')]
ax.legend(handles=handles, loc='upper left', handletextpad=0.3)

# minor ticks only at unlabelled decades
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
ax.tick_params(which='minor', labelleft=False, labelbottom=False)

save(fig, OUT, 'SuppFig1')
s2 = np.array([r['<S^2>'] for r in runs], dtype=float)
print('runs', len(runs), 'max <S^2>', s2.max(), 'n <S^2> < 1e-12:', int((s2 < S2_FLOOR).sum()))
