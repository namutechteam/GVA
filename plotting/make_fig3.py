"""Fig. 3 | Cisplatin CAS(10e,10o) (redrawn with matplotlib from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: deviation of active MP2 natural-orbital occupations from 2 (occupied) or 0 (virtual) (sheet Fig3a_ActiveNOs).
b: |dE| vs CASCI during L-BFGS from the HF start, UCCSD/GVA 1-3 layers (sheet Fig3b_Cisplatin_Convergence),
   with CISD/CCSD/CCSD(T) errors (sheet Table2_Cisplatin_References).
c: ADAPT-VQE with the GVA pool, error after each operator (sheet Fig3c_Cisplatin_ADAPT);
   GVA 1L reference from Runs_Cisplatin (HF start).
Usage: python make_fig3.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullLocator
import numpy as np
import openpyxl

from gva_style import C, COL2, REF, REF_LS, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CHEM = 1.6

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


def mha(x):
    return np.abs(np.asarray(x, dtype=float))   # sheet values are already in mHa


nos = sheet('Fig3a_ActiveNOs')
conv = sheet('Fig3b_Cisplatin_Convergence')
adapt = sheet('Fig3c_Cisplatin_ADAPT')
refs = {r['Method (same CAS(10e,10o))']: r['E - E_CASCI (mHa)'] for r in sheet('Table2_Cisplatin_References')}
gva1 = next(r for r in sheet('Runs_Cisplatin')
            if (r['Ansatz'], r['Layers'], r['Seed (0 = HF start)']) == ('GVA', 1, 0))['dE vs CASCI (mHa)']

fig = plt.figure(figsize=(COL2, 62 / 25.4))
gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.25, 0.95], wspace=0.55,
                      left=0.07, right=0.975, top=0.90, bottom=0.17)

# ---------------- a ----------------
ax = fig.add_subplot(gs[0])
occ = np.array([r['MP2 occupation'] for r in nos])
dev = np.where(occ > 1, 2 - occ, occ)
# neutral greys: orange/blue are reserved for UCCSD/GVA in every other panel
OCC, VIR = '#4D4D4D', '#BDBDBD'
ax.bar(np.arange(len(occ)), dev, 0.78, color=[OCC if o > 1 else VIR for o in occ], lw=0)
ax.legend(handles=[Patch(color=OCC, label='Strongly occupied'), Patch(color=VIR, label='Weakly occupied')],
          loc='upper right', handlelength=1.0, handleheight=0.8, fontsize=5.8)
ax.set_yscale('log')
ax.set_ylim(5e-3, 0.2)
ax.set_xticks(np.arange(len(occ)))
ax.set_xticklabels([r['Irrep (C2v)'].replace('1', '$_1$').replace('2', '$_2$') for r in nos], fontsize=5.8)
ax.tick_params(axis='x', length=0, pad=2.5)
ax.set_xlabel('Active MP2 natural orbital')
ax.set_ylabel('2 − n (occupied), n (virtual)')

# ---------------- b ----------------
axb = fig.add_subplot(gs[1])
it = np.array([r['Iteration'] for r in conv], dtype=float)
for an, ls in (('UCCSD', (0, (3.2, 1.6))), ('GVA', '-')):
    for L, lw in ((1, 0.7), (2, 1.1), (3, 1.6)):
        col = f'{an} {L}L dE (mHa)'
        h = np.array([np.nan if r[col] is None else r[col] for r in conv], dtype=float)
        ok = ~np.isnan(h)
        axb.plot(it[ok], mha(h[ok]), ls=ls, color=C[an], lw=lw, label=f'{an} {L}L', zorder=5 - L)
for key in ('CISD', 'CCSD', 'CCSD(T)'):
    v = abs(refs[key])
    axb.axhline(v, color=REF[key], lw=0.8, ls=REF_LS[key], zorder=1)
    axb.text(1.02, v, key, color=REF[key], fontsize=5.8, va='center', transform=axb.get_yaxis_transform())
axb.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7, zorder=1)
axb.set_yscale('log')
axb.set_xscale('symlog', linthresh=10)
axb.set_ylim(5e-3, 60)
axb.set_xlim(left=0)
axb.set_xlabel('Optimizer iteration')
axb.set_ylabel('|ΔE| relative to CASCI (mHa)')
axb.legend(ncol=2, loc='lower left', fontsize=5.8, handlelength=2.2, columnspacing=1.2)

# ---------------- c ----------------
axc = fig.add_subplot(gs[2])
n = [r['Operators'] for r in adapt]
axc.plot(n, mha([r['dE (mHa)'] for r in adapt]), '-', color=C['ADAPT'], lw=1.0)
axc.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
axc.axhline(gva1, color=C['GVA'], ls=(0, (3.2, 1.6)), lw=0.7)
axc.text(n[-1], gva1 * 0.75, 'GVA 1L', color=C['GVA'], fontsize=5.8, ha='right', va='top')
axc.set_yscale('log')
axc.set_ylim(5e-3, 60)
axc.set_xlabel('ADAPT operators added')
axc.set_ylabel('|ΔE| (mHa)')

# minor ticks only at the unlabelled decades (no 2..9 sub-ticks)
for a in (ax, axb, axc):
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=20))
    a.tick_params(axis='y', which='minor', labelleft=False)
axb.xaxis.set_minor_locator(NullLocator())

fig.canvas.draw()
for a, lab in ((ax, 'a'), (axb, 'b'), (axc, 'c')):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'Fig3')
print('refs', refs, 'GVA 1L', gva1, 'occ', occ.round(4).tolist())
