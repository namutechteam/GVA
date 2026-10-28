"""Supplementary Fig. 7 | Tying ablation for the chains (redrawn with matplotlib from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: single-layer energy from the HF start relative to the untied momentum-screened ansatz (k only);
   symlog axis, linear between +-1e-3 mHa.
b: largest <S^2> over the starts made for each single-layer ansatz (floor 1e-16).
Data: sheet TableS14_PeriodicRuns; r, d, N_k from TableS11_PeriodicSystems.
Usage: python make_suppfig7.py <GVA_Source_Data.xlsx> <out_dir>
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


runs = sheet('TableS14_PeriodicRuns')
systems = sheet('TableS11_PeriodicSystems')
tags = [s['System'] for s in systems]
labels = [f"$N_k$={s['N_k']}" for s in systems]
geo = [f"{s['r (A)']:.2f}/{s['d (A)']:.2f}" for s in systems]
groups = []                                  # (first index, last index, r/d) of runs sharing one geometry
for i, g in enumerate(geo):
    if groups and groups[-1][2] == g:
        groups[-1][1] = i
    else:
        groups.append([i, i, g])


def sel(system, ansatz):
    return [r for r in runs if r['System'] == system and r['Ansatz'] == ansatz and r['Layers'] == 1]


def e_hf(system, ansatz):
    return next(r['E (Ha)'] for r in sel(system, ansatz) if r['Seed'] == 0)


# UCCSD orange and final GVA-k blue keep their meaning; intermediate tying variants in neutral colours
var = [('UCCSD', C['UCCSD'], 'o', False),
       ('k only', C['PG'], 'o', True),
       ('k + spin', C['SPIN'], 'D', True),
       ('k + TR', '#6E6E6E', 'v', True),
       ('GVA-k', C['GVA'], 's', True)]
offs = np.linspace(-0.28, 0.28, len(var))
x = np.arange(len(tags))

fig, axs = plt.subplots(2, 1, figsize=(COL2, 95 / 25.4), sharex=True, gridspec_kw={'hspace': 0.15})
fig.subplots_adjust(left=0.085, right=0.985, top=0.91, bottom=0.17)
for (an, col, mk, filled), dx in zip(var, offs):
    e = [(e_hf(t, an) - e_hf(t, 'k only')) * 1e3 for t in tags]
    s2 = [max(r['<S^2>'] for r in sel(t, an)) for t in tags]
    kw = dict(color=col, mfc=col if filled else 'white', mew=0.7, ms=3.2, ls='none', zorder=3)
    axs[0].plot(x + dx, e, mk, label=an, **kw)
    axs[1].plot(x + dx, np.maximum(s2, 1e-16), mk, **kw)

ax = axs[0]
ax.axhline(0, color='#999999', lw=0.5, zorder=1)
ax.set_yscale('symlog', linthresh=1e-3, linscale=0.6)
ax.set_ylim(-3, 30)
ax.set_ylabel('E − E(k only) (mHa), 1L, HF start')
ax.legend(ncol=5, loc='lower left', bbox_to_anchor=(0.0, 1.0), fontsize=5.8, handletextpad=0.3,
          columnspacing=1.6, borderaxespad=0.3)
for a in axs:                                # thin separators between geometries, as in Fig. 2a
    for (_, i1, _), (j0, _, _) in zip(groups[:-1], groups[1:]):
        a.axvline((i1 + j0) / 2, color='#D0D0D0', lw=0.5, zorder=0)

ax = axs[1]
ax.set_yscale('log')
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=30))
ax.tick_params(axis='y', which='minor', labelleft=False)
ax.set_ylabel(r'largest $\langle\hat{S}^2\rangle$ over starts, 1L')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=5.8, linespacing=1.15)
ax.tick_params(axis='x', length=0, pad=3)
for i0, i1, g in groups:                     # geometry written once per group, under the N_k ticks
    ax.text((i0 + i1) / 2, -0.135, g, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=6.2)
ax.set_xlabel('$r$/$d$ (Å)', labelpad=11)
axs[0].tick_params(axis='x', length=0)

fig.canvas.draw()
x0 = min(a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted()).x0 for a in axs)
for a, lab in zip(axs, 'ab'):
    panel_label(fig, x0 - 0.005, a.get_position().y1 + 0.03, lab)

save(fig, OUT, 'SuppFig7')
for an in ('UCCSD', 'k + spin', 'k + TR', 'GVA-k'):
    print(an, 'dE(Nk5, 2.00/2.00) mHa =', round((e_hf(tags[-1], an) - e_hf(tags[-1], 'k only')) * 1e3, 4),
          ' maxS2 =', '%.2e' % max(r['<S^2>'] for r in sel(tags[-1], an)))
print('k only maxS2 (first, last):', ['%.2e' % max(r['<S^2>'] for r in sel(t, 'k only')) for t in (tags[0], tags[-1])])
