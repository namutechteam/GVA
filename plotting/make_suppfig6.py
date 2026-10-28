"""Supplementary Fig. 6 | Pt-Cl bond stretching in cisplatin, CAS(10e,10o) (redrawn from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a: energies relative to the CASCI minimum along the Pt-Cl distance (Cs);
b: errors relative to CASCI for CISD, CCSD, CCSD(T), UCCSD 1L, GVA 1L and GVA 2L (HF start);
c: reference weight c0^2 (left) and active-space correlation energy -E_corr = E_RHF - E_CASCI (right).
Data: sheet SuppFig6_PtCl_Scan.
Usage: python make_suppfig6.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import numpy as np
import openpyxl

from gva_style import C, REF, REF_LS, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CHEM = 1.6
HA2KCAL = 627.509

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
rows = list(wb['SuppFig6_PtCl_Scan'].iter_rows(values_only=True))
pes = [dict(zip(rows[0], r)) for r in rows[1:] if r[0] is not None]

R = np.array([r['R(Pt-Cl) (A)'] for r in pes], dtype=float)
e0 = min(r['E CASCI'] for r in pes)
col = lambda key: np.array([r[key] for r in pes], dtype=float)   # noqa: E731

# classical methods: shared colours/dash patterns (as in Fig. 3, Fig. 5); RHF grey and CASCI black are plain references
G = {'RHF': '#8C8C8C', **REF, 'CASCI': '#000000'}
LS = {'RHF': '-', **REF_LS, 'CASCI': '-'}

fig, axs = plt.subplots(1, 3, figsize=(COL2, 68 / 25.4), gridspec_kw={'wspace': 0.5})
fig.subplots_adjust(left=0.07, right=0.94, top=0.86, bottom=0.16)

# ---------------- a ----------------
ax = axs[0]
for key, lab in (('E RHF (Ha)', 'RHF'), ('E CCSD', 'CCSD'), ('E CCSD(T)', 'CCSD(T)'), ('E CASCI', 'CASCI')):
    ax.plot(R, (col(key) - e0) * HA2KCAL, ls=LS[lab], color=G[lab], marker='o', ms=2, lw=0.9, label=lab)
ax.plot(R, (col('E GVA 2L') - e0) * HA2KCAL, 'x', color=C['GVA'], ms=4, mew=1.0, label='GVA 2L')
ax.set_xlabel('Pt–Cl distance (Å)')
ax.set_ylabel('Energy relative to CASCI minimum (kcal mol$^{-1}$)')
ax.legend(loc='upper left', fontsize=5.8)

# ---------------- b ----------------
ax = axs[1]
for key, lab in (('CISD err (mHa)', 'CISD'), ('CCSD err (mHa)', 'CCSD'), ('CCSD(T) err (mHa)', 'CCSD(T)')):
    ax.plot(R, np.abs(col(key)), ls=LS[lab], marker='o', ms=2, lw=0.9, color=G[lab], label=lab)
for key, c, ls, lab in (('UCCSD 1L err (mHa)', C['UCCSD'], (0, (3.2, 1.6)), 'UCCSD 1L'),
                        ('GVA 1L err (mHa)', C['GVA'], (0, (3.2, 1.6)), 'GVA 1L'),
                        ('GVA 2L err (mHa)', C['GVA'], '-', 'GVA 2L')):
    ax.plot(R, np.abs(col(key)), ls=ls, marker='s', ms=2, color=c, label=lab, lw=0.9)
ax.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
ax.set_yscale('log')
ax.set_ylim(0.015, 40)
ax.set_xlabel('Pt–Cl distance (Å)')
ax.set_ylabel('|ΔE| relative to CASCI (mHa)')
ax.legend(loc='upper center', fontsize=5.8, ncol=2, handlelength=2.0, columnspacing=1.2)
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=20))

# ---------------- c ----------------
ax = axs[2]
ax.plot(R, col('c0^2'), 'o-', ms=2, lw=0.9, color='#000000', label='$c_0^2$ (CASCI)')
ax.set_xlabel('Pt–Cl distance (Å)')
ax.set_ylabel('Reference weight $c_0^2$')
sec = ax.twinx()
SEC = '#8C8C8C'          # grey, not GVA blue: this curve is not a GVA result
sec.spines['right'].set_visible(True)
sec.spines['right'].set_color(SEC)
sec.plot(R, (col('E RHF (Ha)') - col('E CASCI')) * 1e3, 's', ls=(0, (3.2, 1.6)), ms=2, lw=0.9, color=SEC,
         label='−$E_\\mathrm{corr}$ (active space)')
sec.set_ylabel('−$E_\\mathrm{corr}$, active space (mHa)', color=SEC)
sec.tick_params(axis='y', colors=SEC)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = sec.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=5.8,
          columnspacing=1.2, borderaxespad=0.3)

fig.canvas.draw()
for a, lab in zip(axs, 'abc'):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'SuppFig6')
print('R', R.tolist())
print('a RHF kcal', ((col('E RHF (Ha)') - e0) * HA2KCAL).round(1).tolist())
print('b GVA2L', np.abs(col('GVA 2L err (mHa)')).round(4).tolist())
print('c c0^2', col('c0^2').round(4).tolist(), '-Ecorr', ((col('E RHF (Ha)') - col('E CASCI')) * 1e3).round(2).tolist())
