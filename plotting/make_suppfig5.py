"""Supplementary Fig. 5 | Error budget under a global depolarizing model (redrawn from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: unmitigated energy bias of the cisplatin circuits against the two-qubit depolarizing probability p
   (sheet SuppFig5a_Cisplatin_Bias; CNOT counts in the legend from SuppFig5_NoiseBudget).
b: largest tolerable p without mitigation, with Richardson ZNE and with fidelity rescaling, for the cheapest
   GVA circuit of each system with noiseless |dE| < 1.6 mHa (sheet SuppFig5_NoiseBudget).
Usage: python make_suppfig5.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullLocator
import numpy as np
import openpyxl

from gva_style import C, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CHEM = 1.6
P_HW = 1.84e-3      # average two-qubit infidelity, Moses et al. PRX 13, 041052 (2023), as in the original

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


bias = sheet('SuppFig5a_Cisplatin_Bias')
budget = sheet('SuppFig5_NoiseBudget')
CIS = 'Cisplatin CAS(10e,10o)'
SYSTEMS = ['H2', 'H4', 'LiH', 'BeH2', CIS]
TEX = ['H$_2$', 'H$_4$', 'LiH', 'BeH$_2$', 'Cisplatin']


def cnot(system, an, L):
    return next(r['CNOT'] for r in budget if (r['System'], r['Ansatz'], r['Layers']) == (system, an, L))


fig, axs = plt.subplots(1, 2, figsize=(COL2 * 0.82, 62 / 25.4), gridspec_kw={'wspace': 0.32})
fig.subplots_adjust(left=0.09, right=0.985, top=0.88, bottom=0.2)

# ---------------- a ----------------
ax = axs[0]
p = np.array([r['p (two-qubit depolarizing)'] for r in bias], dtype=float)
for an, L, ls in (('GVA', 1, '-'), ('GVA', 2, '-'), ('UCCSD', 1, (0, (3.2, 1.6))), ('UCCSD', 2, (0, (3.2, 1.6)))):
    nc = cnot(CIS, an, L)
    ax.plot(p, [r[f'{an} {L}L bias (mHa)'] for r in bias], ls=ls, color=C[an], lw=0.8 if L == 1 else 1.4,
            label=f'{an} {L}L ({nc / 1e3:.0f}k CNOT)')
ax.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
ax.axvline(P_HW, color='#777777', ls=(0, (3.2, 1.6)), lw=0.7)
ax.text(P_HW * 0.8, 3e4, '1.8×10$^{-3}$', ha='right', fontsize=5.8, color='#777777')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-4, 2e5)
ax.set_xlabel('Two-qubit depolarizing probability p')
ax.set_ylabel('Energy bias, unmitigated (mHa)')
ax.set_title('Cisplatin CAS(10e,10o)', fontsize=6.6, pad=3)
ax.legend(loc='lower left', bbox_to_anchor=(0.26, 0.01), fontsize=5.6)   # empty area below the curves

# ---------------- b ----------------
ax = axs[1]
sel = []
for s in SYSTEMS:
    rr = sorted((r for r in budget if r['System'] == s and r['Ansatz'] == 'GVA'
                 and abs(r['Noiseless dE (mHa)']) < CHEM), key=lambda r: r['CNOT'])
    sel.append(rr[0])
xx = np.arange(len(sel))
w = 0.27
# mitigation schemes are not methods: three greys, light (none) to dark (strongest)
for k, (key, col, lab) in enumerate((('p_max unmitigated', '#C8C8C8', 'Unmitigated, bias ≤ 1.6 mHa'),
                                     ('p_max Richardson (1,3,5)', '#8C8C8C', 'Richardson ZNE (λ = 1, 3, 5)'),
                                     ('p_max for 1/F^2 <= 100', '#3A3A3A', 'Rescaling, overhead ≤ 100'))):
    ax.bar(xx + (k - 1) * w, [r[key] for r in sel], w, color=col, lw=0, label=lab)
ax.set_yscale('log')
ax.set_ylim(1e-10, 10)
ax.set_xticks(xx)
ax.set_xticklabels([f"{t}\n{r['Layers']}L, {r['CNOT']:,}" for t, r in zip(TEX, sel)], fontsize=5.8,
                   linespacing=1.15)
ax.tick_params(axis='x', length=0, pad=3)
ax.set_ylabel('Largest tolerable p')
ax.axhline(P_HW, color='#777777', ls=(0, (3.2, 1.6)), lw=0.7)
ax.legend(loc='upper right', fontsize=5.8, handlelength=1.0, handleheight=0.8)

# minor ticks only at unlabelled decades
for a in axs:
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
    a.tick_params(axis='y', which='minor', labelleft=False)
axs[0].xaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
axs[0].tick_params(axis='x', which='minor', labelbottom=False)

fig.canvas.draw()
for a, lab in zip(axs, 'ab'):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'SuppFig5')
print('a CNOT', [cnot(CIS, an, L) for an, L in (('GVA', 1), ('GVA', 2), ('UCCSD', 1), ('UCCSD', 2))],
      'bias@1e-10', [bias[0][k] for k in list(bias[0])[1:]])
print('b', [(r['System'], r['Layers'], r['CNOT'], '%.2e %.2e %.2e' % (r['p_max unmitigated'],
       r['p_max Richardson (1,3,5)'], r['p_max for 1/F^2 <= 100'])) for r in sel])
