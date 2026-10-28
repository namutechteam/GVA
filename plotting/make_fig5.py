"""Fig. 5 | VQE on periodic hydrogen chains (redrawn with matplotlib from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: dimerized chain (r = 0.74, d = 1.50 A) vs N_k; b: uniform chains at three spacings, N_k = 5 (solid) and 3 (faint);
c: best error vs CNOT gates for the four N_k = 5 chains. Single-layer and ADAPT-VQE: HF start; two layers:
median and min..max over starts (a, b) or best start (c). CCSD/CCSD(T) as absolute values.
Data: sheet Fig5_SuppFig7_Data; c0^2 from TableS11_PeriodicSystems.
Usage: python make_fig5.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator, LogLocator, NullLocator
import numpy as np
import openpyxl

from gva_style import C, COL2, REF as RC, REF_LS, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CHEM = 1.6

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


data = sheet('Fig5_SuppFig7_Data')
c0sq = {r['System']: r['c0^2'] for r in sheet('TableS11_PeriodicSystems')}


def tag(r, d, nk):
    return f'H2chain_r{r:.2f}_d{d:.2f}_Nk{nk}'


def row(system, method, layers):
    return next(x for x in data if (x['System'], x['Method'], x['Layers']) == (system, method, layers))


def err(system, method, layers):
    """(min, median, max) in mHa: HF start for 1L and ADAPT, over starts for 2L."""
    x = row(system, method, layers)
    if layers in (1, None):
        e = abs(x['Error, HF start (mHa)'])
        return e, e, e
    return tuple(abs(x[k]) for k in ('Error min (mHa)', 'Error median (mHa)', 'Error max (mHa)'))


spec = [('UCCSD', 1, C['UCCSD'], 'o', False, 'UCCSD 1L'),
        ('GVA-k', 1, C['GVA'], 'o', True, 'GVA-k 1L'),
        ('UCCSD', 2, C['UCCSD'], 's', False, 'UCCSD 2L'),
        ('GVA-k', 2, C['GVA'], 's', True, 'GVA-k 2L'),
        ('ADAPT-VQE', None, C['ADAPT'], '^', True, 'ADAPT (GVA-k pool)')]
REF = [(k, RC[k]) for k in ('CCSD', 'CCSD(T)')]   # classical references: shared colours/dash patterns
offs = np.linspace(-0.2, 0.2, len(spec))


def point(ax, x, lo, me, hi, col, mk, filled, alpha=1.0):
    ax.errorbar(x, me, yerr=[[me - lo], [hi - me]], fmt=mk, color=col, mfc=col if filled else 'white',
                mew=0.7, ms=3.2, lw=0.6, capsize=1.2, capthick=0.6, alpha=alpha, zorder=3)


fig = plt.figure(figsize=(COL2, 64 / 25.4))
gs = fig.add_gridspec(1, 3, wspace=0.45, left=0.07, right=0.985, top=0.88, bottom=0.23)

# ---------------- a ----------------
ax = fig.add_subplot(gs[0])
nks = [2, 3, 4, 5]
for (an, L, col, mk, filled, _), dx in zip(spec, offs):
    for nk in nks:
        point(ax, nk + dx, *err(tag(0.74, 1.50, nk), an, L), col, mk, filled)
for key, shade in REF:
    ax.plot(nks, [err(tag(0.74, 1.50, nk), key, None)[1] for nk in nks], ls=REF_LS[key], color=shade, lw=0.9,
            zorder=2)
ax.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
ax.set_yscale('log')
ax.set_xticks(nks)
ax.set_xlabel('$N_k$ (dimerized chain)')
ax.set_ylabel('|E − E$_\mathrm{FCI}$| (mHa)')

# ---------------- b ----------------
axb = fig.add_subplot(gs[1])
Rs = [1.0, 1.5, 2.0]
for (an, L, col, mk, filled, _), dx in zip(spec, offs):
    for nk, solid in ((5, True), (3, False)):
        for i, R in enumerate(Rs):
            point(axb, i + dx + (0.0 if solid else 0.04), *err(tag(R, R, nk), an, L), col, mk,
                  filled and solid, alpha=1.0 if solid else 0.45)
for key, shade in REF:
    axb.plot(range(len(Rs)), [err(tag(R, R, 5), key, None)[1] for R in Rs], ls=REF_LS[key], color=shade, lw=0.9,
             zorder=2)
axb.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
axb.set_yscale('log')
axb.set_xticks(range(len(Rs)))
axb.set_xticklabels([f'{R:.1f}' for R in Rs])
axb.set_xlabel('H–H spacing (Å), uniform chain')
axb.set_xlim(-0.4, 2.4)
axb.set_ylabel('|E − E$_\mathrm{FCI}$| (mHa)')
for i, R in enumerate(Rs):
    axb.text(i, 1.02, f'$c_0^2$ = {c0sq[tag(R, R, 5)]:.2f}', transform=axb.get_xaxis_transform(),
             ha='center', va='bottom', fontsize=5.8, color='#777777')
axb.legend(handles=[Line2D([], [], marker='o', ls='none', color='#555555', mew=0, ms=3.4, label='$N_k$ = 5'),
                    Line2D([], [], marker='o', ls='none', color='#555555', mew=0, alpha=0.45, ms=3.4, label='$N_k$ = 3')],
           loc='upper left', fontsize=5.6, handletextpad=0.2, borderaxespad=0.3, labelspacing=0.3)

# ---------------- c ----------------
axc = fig.add_subplot(gs[2])
CHAINS = (tag(0.74, 1.50, 5), tag(1.0, 1.0, 5), tag(1.5, 1.5, 5), tag(2.0, 2.0, 5))
SHADES = (1.0, 0.62, 0.38, 0.22)       # dark to light: dimerized, 1.0, 1.5, 2.0 A (wider steps than before)
for an, L, col, mk, filled, _ in spec:
    for t, shade in zip(CHAINS, SHADES):
        x = row(t, an, L)
        axc.plot(x['CNOT'], err(t, an, L)[0], mk, color=col, mfc=col if filled else 'white', mew=0.9,
                 alpha=shade, ms=3.5, zorder=3)
axc.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
axc.set_xscale('log')
axc.set_yscale('log')
axc.xaxis.set_major_locator(FixedLocator([1e4, 3e4, 1e5, 3e5]))
axc.xaxis.set_major_formatter(FixedFormatter(['10$^4$', '3×10$^4$', '10$^5$', '3×10$^5$']))
axc.xaxis.set_minor_locator(NullLocator())
axc.set_xlim(1e4, 3.5e5)
axc.legend(handles=[Line2D([], [], marker='s', ls='none', color='#333333', mew=0, alpha=a, ms=3.4, label=lab)
                    for a, lab in zip(SHADES, ('dimerized', '1.0 Å', '1.5 Å', '2.0 Å'))],
           loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=4, fontsize=5.6, handletextpad=0.1,
           columnspacing=0.9, borderaxespad=0.2)
axc.set_xlabel('CNOT gates (20 qubits, $N_k$ = 5)')
axc.set_ylabel('best |E − E$_\mathrm{FCI}$| (mHa)')

# minor ticks only at unlabelled decades
for a in (ax, axb, axc):
    a.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=20))
    a.tick_params(axis='y', which='minor', labelleft=False)

handles = [Line2D([], [], marker=mk, color=col, mfc=col if f else 'white', mew=0.7, ls='none', ms=3.5, label=lab)
           for _, _, col, mk, f, lab in spec]
handles += [Line2D([], [], color=shade, lw=0.9, ls=REF_LS[key], label=key) for key, shade in REF]
handles.append(Line2D([], [], color='#777777', ls=(0, (1.2, 1.6)), lw=0.7, label='1.6 mHa'))
fig.legend(handles=handles, loc='lower center', ncol=8, bbox_to_anchor=(0.5, 0.0), fontsize=5.8,
           handletextpad=0.3, columnspacing=1.1)

fig.canvas.draw()
for a, lab in ((ax, 'a'), (axb, 'b'), (axc, 'c')):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'Fig5')
print('c0^2', [round(c0sq[tag(R, R, 5)], 2) for R in Rs],
      'a GVA-k 2L Nk5', err(tag(0.74, 1.50, 5), 'GVA-k', 2))
