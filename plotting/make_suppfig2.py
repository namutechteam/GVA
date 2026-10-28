"""Supplementary Fig. 2 | Resource scaling for cisplatin active spaces (redrawn from GVA_Source_Data.xlsx).

Content identical to the original figure; only the design differs.
Parameters (a) and CNOT gates (b) per layer for UCCSD and GVA against the number of active orbitals
(sheet SuppFig2_Scaling); slopes fitted in log-log to the six largest spaces; grey squares (right axis):
GVA/UCCSD ratio.
Usage: python make_suppfig2.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullLocator
import numpy as np
import openpyxl

from gva_style import C, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])
GREY = '#808080'

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)


def sheet(name):
    rows = list(wb[name].iter_rows(values_only=True))
    return [dict(zip(rows[0], r)) for r in rows[1:] if any(v is not None for v in r)]


scal = sheet('SuppFig2_Scaling')
n = np.array([r['Active orbitals'] for r in scal], dtype=float)


def fit_slope(x, y):
    return np.polyfit(np.log(x), np.log(y), 1)[0]


fig, axs = plt.subplots(1, 2, figsize=(COL2 * 0.75, 60 / 25.4))
fig.subplots_adjust(left=0.10, right=0.90, top=0.90, bottom=0.18, wspace=0.62)
slopes = {}
for ax, ku, kg, yl in ((axs[0], 'UCCSD params', 'GVA params', 'Parameters per layer'),
                       (axs[1], 'UCCSD CNOT/layer', 'GVA CNOT/layer', 'CNOT gates per layer')):
    yu = np.array([r[ku] for r in scal], dtype=float)
    yg = np.array([r[kg] for r in scal], dtype=float)
    su, sg = fit_slope(n[3:], yu[3:]), fit_slope(n[3:], yg[3:])
    slopes[yl] = (su, sg)
    sec = ax.twinx()
    sec.plot(n, yg / yu, 's', ls=(0, (1.2, 1.6)), color=GREY, ms=2.5, lw=0.7, zorder=1)
    sec.set_ylim(0, 0.6)
    sec.set_ylabel('GVA / UCCSD', color=GREY)
    sec.spines['right'].set_visible(True)
    sec.spines['right'].set_color(GREY)
    sec.tick_params(axis='y', colors=GREY)
    ax.set_zorder(sec.get_zorder() + 1)      # method curves above the ratio
    ax.patch.set_visible(False)
    ax.plot(n, yu, 'o-', color=C['UCCSD'], mfc='white', mew=0.8, ms=3.5, lw=1.0, label=f'UCCSD (slope {su:.2f})')
    ax.plot(n, yg, 'o-', color=C['GVA'], ms=3.5, lw=1.0, label=f'GVA (slope {sg:.2f})')
    ax.set_xscale('log')
    ax.set_yscale('log')
    tk = [4, 6, 8, 10, 12, 16, 20]
    ax.set_xticks(tk)
    ax.set_xticklabels([str(v) for v in tk])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(1.0,), numticks=20))
    ax.tick_params(axis='y', which='minor', labelleft=False)
    ax.set_xlabel('Active orbitals (qubits = 2×)')
    ax.set_ylabel(yl)
    ax.legend(loc='lower right', handlelength=1.8)

fig.canvas.draw()
for a, lab in zip(axs, 'ab'):
    bb = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    panel_label(fig, bb.x0 - 0.005, 0.99, lab)

save(fig, OUT, 'SuppFig2')
print('slopes', {k: tuple(round(v, 2) for v in s) for k, s in slopes.items()})
