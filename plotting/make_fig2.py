"""Fig. 2 | Small-molecule benchmark (redrawn with matplotlib from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: |dE| vs FCI for UCCSD/GVA 1L/2L (median over 5 starts, bar = min..max) and ADAPT-VQE (HF start),
   values below 1e-6 mHa drawn at that floor (sheet Runs_SmallMolecules).
b, c: parameters and CNOT gates per layer at the reference geometries (sheet Table1_SmallMolecules).
d: |dE| during L-BFGS for BeH2 x2.0, HF start (sheet Fig2d_BeH2x2_Convergence).
Usage: python make_fig2.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullLocator
import numpy as np
import openpyxl

from gva_style import C, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
FLOOR, CHEM = 1e-6, 1.6
MOLS = ['H2', 'H4', 'LiH', 'BeH2']
SCALES = [1.0, 1.5, 2.0]
TEX = {'H2': r'H$_2$', 'H4': r'H$_4$', 'LiH': 'LiH', 'BeH2': r'BeH$_2$'}

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


runs = sheet('Runs_SmallMolecules')
tab1 = sheet('Table1_SmallMolecules')
conv = sheet('Fig2d_BeH2x2_Convergence')


def mha(x):
    return np.maximum(np.abs(np.asarray(x, dtype=float)), FLOOR)   # sheet values are already in mHa


def sel(**kw):
    out = []
    for r in runs:
        if all(r[{'system': 'System', 'scale': 'Bond scale', 'ansatz': 'Ansatz', 'layers': 'Layers'}[k]] == v
               for k, v in kw.items()):
            out.append(r)
    return out


fig = plt.figure(figsize=(COL2, 118 / 25.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1], hspace=0.40, wspace=0.48,
                      left=0.075, right=0.985, top=0.93, bottom=0.085)

# ---------------- a ----------------
ax = fig.add_subplot(gs[0, :])
spec = [('UCCSD', 1, 'o', False), ('GVA', 1, 'o', True), ('UCCSD', 2, 's', False), ('GVA', 2, 's', True)]
offs = [-0.27, -0.09, 0.09, 0.27]
xs, labels, groups = [], [], []
x = 0.0
for m in MOLS:
    g0 = x
    for sc in SCALES:
        for (an, L, mk, filled), dx in zip(spec, offs):
            e = mha([r['dE vs FCI (mHa)'] for r in sel(system=m, scale=sc, ansatz=an, layers=L)])
            col = C[an]
            ax.plot([x + dx] * 2, [e.min(), e.max()], color=col, lw=0.6, solid_capstyle='butt', zorder=2)
            ax.plot(x + dx, np.median(e), mk, ms=3.4, mec=col, mew=0.7, mfc=col if filled else 'white', zorder=3)
        ad = sel(system=m, scale=sc, ansatz='ADAPT-VQE')
        ax.plot(x + 0.42, mha(ad[0]['dE vs FCI (mHa)']), '^', ms=3.8, color=C['ADAPT'], mew=0, zorder=3)
        xs.append(x)
        labels.append(f'{TEX[m]}\n' + f'×{sc:.1f}')
        x += 1
    groups.append((g0 - 0.5, x - 0.5))
    x += 0.5
for (_, a1), (b0, _) in zip(groups[:-1], groups[1:]):      # thin separators between molecules
    ax.axvline((a1 + b0) / 2, color='#D0D0D0', lw=0.5, zorder=0)
ax.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7, zorder=1)
ax.text(xs[-1] + 0.68, CHEM * 1.35, '1.6 mHa', fontsize=5.8, color='#777777', ha='right', va='bottom')
ax.set_yscale('log')
ax.set_ylim(5e-7, 200)
ax.set_yticks([1e-6, 1e-4, 1e-2, 1e0, 1e2])
ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=6.3, linespacing=1.05)
ax.tick_params(axis='x', length=0, pad=3)
ax.set_xlim(-0.6, xs[-1] + 0.75)
ax.set_ylabel('|ΔE| relative to FCI (mHa)')
handles = [Line2D([], [], marker=mk, ls='', ms=3.4, mec=C[an], mfc=C[an] if f else 'white', mew=0.7,
                  label=f'{an} {L}L') for an, L, mk, f in spec]
handles.append(Line2D([], [], marker='^', ls='', ms=3.8, color=C['ADAPT'], mew=0, label='ADAPT-VQE'))
leg = ax.legend(handles=handles, ncol=5, loc='lower left', bbox_to_anchor=(0.0, 1.0), columnspacing=1.6,
          handletextpad=0.3, borderaxespad=0.2)
fig.canvas.draw()
lb = leg.get_window_extent().transformed(fig.transFigure.inverted())
panel_label(fig, 0.005, lb.y1 + 0.004, 'a')

# ---------------- b, c ----------------
ref = {r['System']: r for r in tab1 if r['Bond scale'] == 1}
xx = np.arange(len(MOLS))
w = 0.38
axb = fig.add_subplot(gs[1, 0])
pu = [ref[m]['UCCSD params/layer'] for m in MOLS]
pg = [ref[m]['GVA params/layer'] for m in MOLS]
axb.bar(xx - w / 2, pu, w, color=C['UCCSD'], label='UCCSD', lw=0)
axb.bar(xx + w / 2, pg, w, color=C['GVA'], label='GVA', lw=0)
for i in range(len(MOLS)):
    axb.text(xx[i] - w / 2, pu[i] * 1.18, str(pu[i]), ha='center', va='bottom', fontsize=5.6)
    axb.text(xx[i] + w / 2, pg[i] * 1.18, str(pg[i]), ha='center', va='bottom', fontsize=5.6)
axb.set_yscale('log')
axb.set_ylim(0.6, 800)
axb.set_ylabel('Parameters per layer')
axb.legend(loc='upper left', handlelength=1.0, handleheight=0.8)

axc = fig.add_subplot(gs[1, 1])
cu = [ref[m]['UCCSD CNOT/layer'] for m in MOLS]
cg = [ref[m]['GVA CNOT/layer'] for m in MOLS]
axc.bar(xx - w / 2, cu, w, color=C['UCCSD'], lw=0)
axc.bar(xx + w / 2, cg, w, color=C['GVA'], lw=0)
for i in range(len(MOLS)):
    red = 100 * (1 - cg[i] / cu[i])
    axc.text(xx[i], max(cu[i], cg[i]) * 1.45, f'−{red:.0f}%', ha='center', va='bottom', fontsize=5.6)
axc.set_yscale('log')
axc.set_ylim(20, 2e5)
axc.set_ylabel('CNOT gates per layer')
for a in (axb, axc):
    a.set_xticks(xx)
    a.set_xticklabels([TEX[m] for m in MOLS])
    a.tick_params(axis='x', length=0, pad=3)

# ---------------- d ----------------
axd = fig.add_subplot(gs[1, 2])
it = np.array([r['Iteration'] for r in conv], dtype=float)
for an, L, ls in (('UCCSD', 1, (0, (3.2, 1.6))), ('GVA', 1, '-'), ('UCCSD', 2, (0, (3.2, 1.6))), ('GVA', 2, '-')):
    col = f'{an} {L}L dE (mHa)'
    h = np.array([np.nan if r[col] is None else r[col] for r in conv], dtype=float)
    ok = ~np.isnan(h)
    # as in the original: 1L thin, 2L thick; 1L drawn on top where the curves coincide
    axd.plot(it[ok], mha(h[ok]), ls=ls, color=C[an], lw=0.8 if L == 1 else 1.5,
             label=f'{an} {L}L', zorder=3 if L == 1 else 2)
axd.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
axd.set_yscale('log')
axd.set_xscale('symlog', linthresh=10)
axd.set_ylim(5e-7, 500)
axd.set_yticks([1e-6, 1e-4, 1e-2, 1e0, 1e2])
axd.set_xlim(left=0)
axd.set_xlabel('Optimizer iteration')
axd.set_ylabel('|ΔE| (mHa)')
axd.set_title(r'BeH$_2$, bonds ×2.0', fontsize=6.6, pad=3)
axd.legend(loc='lower left', fontsize=5.8, handlelength=2.2)

for a, lab in ((axb, 'b'), (axc, 'c'), (axd, 'd')):
    fig.canvas.draw()
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.01, a.get_position().y1 + 0.035, lab)

# minor ticks only at the unlabelled decades (no 2..9 sub-ticks)
for a in (ax, axd):
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=20))
    a.tick_params(axis='y', which='minor', labelleft=False)
for a in (axb, axc):
    a.yaxis.set_minor_locator(NullLocator())
axd.xaxis.set_minor_locator(NullLocator())

save(fig, OUT, 'Fig2')
print('b params', list(zip(MOLS, pu, pg)), 'c CNOT', list(zip(MOLS, cu, cg)))
