"""Supplementary Fig. 3 | Accuracy against circuit resources for cisplatin (redrawn from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
a, b: |dE| vs CASCI against CNOT count and parameter count for UCCSD and GVA with 1-3 layers
      (median, bar = all starts; symbol size grows with layers), screening-only and tying-only 1L (diamonds),
      ADAPT-VQE (triangle) (sheet Runs_Cisplatin); CCSD and CCSD(T) errors (sheet Table2_Cisplatin_References).
Usage: python make_suppfig3.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import openpyxl

from gva_style import C, REF as RC, REF_LS, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CHEM, FLOOR = 1.6, 1e-6

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


def mha(x):
    return np.maximum(np.abs(np.asarray(x, dtype=float)), FLOOR)   # sheet values are already in mHa


runs = sheet('Runs_Cisplatin')
refs = {r['Method (same CAS(10e,10o))']: r['E - E_CASCI (mHa)'] for r in sheet('Table2_Cisplatin_References')}
COLS = {'UCCSD': C['UCCSD'], 'GVA': C['GVA'], 'PG only': C['PG'], 'spin only': C['SPIN']}
REF = [(k, RC[k]) for k in ('CCSD', 'CCSD(T)')]   # classical references: shared colours/dash patterns
DASHDOT = (0, (4, 1.5, 1, 1.5))

fig, axs = plt.subplots(1, 2, figsize=(COL2 * 0.82, 62 / 25.4), gridspec_kw={'wspace': 0.32})
fig.subplots_adjust(left=0.09, right=0.98, top=0.90, bottom=0.25)
for ax, key, xl in ((axs[0], 'CNOT', 'CNOT gates (ladder compilation)'),
                    (axs[1], 'Parameters', 'Variational parameters')):
    for an, mk in (('UCCSD', 'o'), ('GVA', 'o'), ('PG only', 'D'), ('spin only', 'D')):
        for L in (1, 2, 3):
            rr = [r for r in runs if r['Ansatz'] == an and r['Layers'] == L]
            if not rr:
                continue
            e = mha([r['dE vs CASCI (mHa)'] for r in rr])
            col = COLS[an]
            abl = an in ('PG only', 'spin only')
            me = np.median(e)
            spread = e.max() - e.min() > 1e-4
            ax.errorbar(rr[0][key], me, yerr=[[me - e.min()], [e.max() - me]] if spread else None,
                        fmt=mk, ms=9 if abl else 3.5 + L, color=col,
                        mfc='none' if (abl or an == 'UCCSD') else col,
                        mew=0.9, capsize=1.5, capthick=0.6, lw=0.6, zorder=2 if abl else 3)
            if an in ('UCCSD', 'GVA'):
                ax.annotate(f'{L}L', (rr[0][key], me), textcoords='offset points',
                            xytext=(5, -9) if (an, L) == ('GVA', 2) else (5, 3),   # clear of CCSD(T) line
                            fontsize=5.8, color=col)
    ad = next(r for r in runs if r['Ansatz'] == 'ADAPT-VQE')
    ax.plot(ad[key], mha(ad['dE vs CASCI (mHa)']), '^', color=C['ADAPT'], ms=5, mew=0)
    ax.axhline(CHEM, color='#777777', ls=(0, (1.2, 1.6)), lw=0.7)
    for k2, shade in REF:
        ax.axhline(abs(refs[k2]), color=shade, lw=0.8, ls=REF_LS[k2], zorder=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xl)
    ax.set_ylabel('|ΔE| relative to CASCI (mHa)')
    ax.minorticks_off()
    if key == 'CNOT':
        ax.set_xticks([2e4, 5e4, 1e5, 2e5, 4e5])
        ax.set_xticklabels(['2×10$^4$', '5×10$^4$', '10$^5$', '2×10$^5$', '4×10$^5$'])
    else:
        ax.set_xticks([100, 300, 1000, 3000])
        ax.set_xticklabels(['100', '300', '1,000', '3,000'])
    ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1])
    ax.set_yticklabels(['0.01', '0.03', '0.1', '0.3', '1'])

handles = [Line2D([], [], marker='o', ls='', mec=C['UCCSD'], mfc='white', mew=0.9, label='UCCSD'),
           Line2D([], [], marker='o', ls='', color=C['GVA'], label='GVA'),
           Line2D([], [], marker='D', ls='', mec=C['PG'], mfc='none', mew=0.9, ms=6, label='Screening only 1L'),
           Line2D([], [], marker='D', ls='', mec=C['SPIN'], mfc='none', mew=0.9, ms=6, label='Tying only 1L'),
           Line2D([], [], marker='^', ls='', color=C['ADAPT'], mew=0, label='ADAPT-VQE')]
handles += [Line2D([], [], ls=REF_LS[k2], lw=0.8, color=shade, label=k2) for k2, shade in REF]
fig.legend(handles=handles, loc='lower center', ncol=7, fontsize=5.8, bbox_to_anchor=(0.5, 0.0),
           handletextpad=0.3, columnspacing=1.0)

fig.canvas.draw()
for a, lab in zip(axs, 'ab'):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'SuppFig3')
for an in ('UCCSD', 'GVA'):
    for L in (1, 2, 3):
        e = [r['dE vs CASCI (mHa)'] for r in runs if r['Ansatz'] == an and r['Layers'] == L]
        print(an, L, 'n', len(e), 'median', round(float(np.median(e)), 4), 'min', round(min(e), 4), 'max', round(max(e), 4))
print('ADAPT', ad['CNOT'], ad['Parameters'], round(ad['dE vs CASCI (mHa)'], 4),
      'CCSD', round(refs['CCSD'], 4), 'CCSD(T)', round(refs['CCSD(T)'], 4))
