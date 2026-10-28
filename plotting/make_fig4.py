"""Fig. 4 | Savings from symmetry against the order of the symmetry group (redrawn from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a, b: GVA/UCCSD parameter and CNOT ratios against h (sheet Fig4_RatioVsGroupOrder):
      1D chain 1+1 bands (N_k >= 2), 3D mesh 2+2 and 8+8 bands (N_k >= 8), small molecules,
      cisplatin CAS with >= 10 active orbitals; dashed lines 1/(2h) and 1/h.
c: parameters and CNOT gates per layer for 3D meshes, 2+2 bands per cell (sheet TableS13_PeriodicScaling).
Usage: python make_fig4.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import numpy as np
import openpyxl

from gva_style import C, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


ratio = sheet('Fig4_RatioVsGroupOrder')
scal = [r for r in sheet('TableS13_PeriodicScaling') if r['Set'] == 'three_d' and r['o per cell'] == 2]


def pick(typ, hmin=1, cas_min=None):
    out = [r for r in ratio if r['Type'] == typ and r['Group order h'] >= hmin]
    if cas_min is not None:     # cisplatin CAS(ne,no): keep no >= cas_min
        out = [r for r in out if int(r['System / mesh'].split(',')[1].rstrip('o)')) >= cas_min]
    return out


# crystals in greys (all lie on the line), molecules in purple (above the line)
MOL = C['SPIN']
series = [
    (pick('1D mesh, 1+1 bands', 2), dict(marker='o', ls='-', color='#7A7A7A', ms=2.5, lw=0.7),
     '1D chain, 1+1 bands'),
    (pick('3D mesh, 2+2 bands', 8), dict(marker='s', ls='-', color='#222222', ms=3, lw=0.7),
     '3D mesh, 2+2 bands'),
    (pick('3D mesh, 8+8 bands, 8+8 bands', 8), dict(marker='D', ls='', color='#222222', ms=3.2, mfc='white',
                                                     mew=0.8), '3D mesh, 8+8 bands'),
    (pick('molecule'), dict(marker='^', ls='', color=MOL, ms=3.8, mew=0), 'small molecules'),
    (pick('cisplatin CAS (C2v)', cas_min=10), dict(marker='v', ls='', color=MOL, ms=3.8, mfc='white', mew=0.8),
     'cisplatin CAS'),
]

fig, axs = plt.subplots(1, 3, figsize=(COL2, 60 / 25.4), gridspec_kw={'wspace': 0.48})
fig.subplots_adjust(left=0.07, right=0.985, top=0.90, bottom=0.18)
hh = np.logspace(0, np.log10(160), 50)
for k, (ax, col, law, lab, txt) in enumerate((
        (axs[0], 'Param. ratio', 0.5 / hh, 'parameters', '$1/(2h)$'),
        (axs[1], 'CNOT ratio', 1 / hh, 'CNOT gates', '$1/h$'))):
    ax.plot(hh, law, color='#999999', lw=0.7, ls=(0, (3.2, 1.6)), zorder=1)
    for rows, kw, name in series:
        ax.plot([r['Group order h'] for r in rows], [r[col] for r in rows], label=name, zorder=3, **kw)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('symmetry group order $h$ ($h = N_k$ for crystals)')
    ax.set_ylabel(f'GVA / UCCSD, {lab}')
    ax.set_xlim(1.5, 160)
    ax.text(0.97, 0.93, txt, transform=ax.transAxes, ha='right', va='top', color='#777777')
axs[1].legend(loc='lower left', fontsize=5.5, handlelength=1.3, handletextpad=0.4, labelspacing=0.35, borderaxespad=0.2)   # empty corner; shared by a and b

# ---------------- c ----------------
ax = axs[2]
q = [r['Qubits'] for r in scal]
ax.plot(q, [r['UCCSD params'] for r in scal], 'o-', color=C['UCCSD'], mfc='white', ms=3, lw=0.7,
        label='UCCSD parameters')
ax.plot(q, [r['GVA-k params'] for r in scal], 'o-', color=C['GVA'], ms=3, lw=0.7, label='GVA-k parameters')
ax.plot(q, [r['UCCSD CNOT/layer'] for r in scal], 's', ls=(0, (3.2, 1.6)), color=C['UCCSD'], mfc='white',
        ms=3, lw=0.7, label='UCCSD CNOT')
ax.plot(q, [r['GVA-k CNOT/layer'] for r in scal], 's', ls=(0, (3.2, 1.6)), color=C['GVA'], ms=3, lw=0.7,
        label='GVA-k CNOT')
for r in scal:
    if r['N_k'] > 1:
        ax.annotate(f"{r['Mesh'][0]}$^3$", (r['Qubits'], r['GVA-k params']), textcoords='offset points',
                    xytext=(2, -7), fontsize=5.6, color='#777777')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('qubits (3D mesh, 2+2 bands per cell)')
ax.set_ylabel('count per layer')
ax.legend(loc='upper left', fontsize=5.8)

# minor ticks only at unlabelled decades
for a in axs:
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
    a.xaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
    a.tick_params(which='minor', labelleft=False, labelbottom=False)

fig.canvas.draw()
for a, lab in zip(axs, 'abc'):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'Fig4')
print('molecules', [(r['System / mesh'], r['Group order h'], round(r['Param. ratio'], 3), round(r['CNOT ratio'], 3))
                    for r in series[3][0]])
print('cisplatin', [r['System / mesh'] for r in series[4][0]], 'qubits', q)
