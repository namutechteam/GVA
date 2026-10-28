"""Fig. 1 | Construction of GVA (redrawn with matplotlib from GVA_Source_Data.xlsx).

a: workflow schematic (numbered step boxes with tinted title bars and equations).
b: parameters and CNOT gates per layer for cisplatin CAS(10e,10o) in C2v after each step
   (values read from sheet Runs_Cisplatin, single-layer Hartree-Fock-start rows).
Usage: python make_fig1.py <GVA_Source_Data.xlsx> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import openpyxl

from gva_style import C, COL2, panel_label, save

XLSX, OUT = Path(sys.argv[1]), Path(sys.argv[2])

# ---------------- data (panel b) ----------------
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
rows = list(wb['Runs_Cisplatin'].iter_rows(values_only=True))
col = {h: i for i, h in enumerate(rows[0])}


def pick(ansatz, layers=1, seed=0):
    for r in rows[1:]:
        if (r[col['Ansatz']], r[col['Layers']], r[col['Seed (0 = HF start)']]) == (ansatz, layers, seed):
            return int(r[col['Parameters']]), int(r[col['CNOT']])
    raise KeyError(ansatz)


steps = [('UCCSD', 'UCCSD', C['UCCSD']),
         ('+ point-group screening', 'PG only', C['PG']),
         ('+ spin-flip tying (GVA)', 'GVA', C['GVA'])]
vals = [(lab, *pick(key), colr) for lab, key, colr in steps]

# ---------------- figure ----------------
fig = plt.figure(figsize=(COL2, 62 / 25.4))
axa = fig.add_axes([0.0, 0.0, 0.60, 1.0])
axa.set_xlim(0, 100)
axa.set_ylim(0, 38)
axa.axis('off')

BW, BH, HB = 29.5, 12.6, 3.3            # box width, height, title-bar height (data units)
xs = [1.5, 35.25, 69.0]
ytop, ybot = 21.0, 3.2


def step(x, y, n, title, body, color, lw=0.6):
    """Plain box: thin coloured outline, bold title (number in the method colour), body text."""
    axa.add_patch(FancyBboxPatch((x, y), BW, BH, boxstyle='round,pad=0,rounding_size=0.7',
                                 fc='white', ec=color, lw=lw))
    axa.text(x + 1.3, y + BH - 2.1, f'{n}', color=color, fontsize=6.4, fontweight='bold', va='center', ha='left')
    axa.text(x + 3.2, y + BH - 2.1, title, color='#111111', fontsize=6.4, fontweight='bold', va='center', ha='left')
    axa.text(x + BW / 2, y + (BH - 3.6) / 2, body, fontsize=6.2, color='#111111',
             va='center', ha='center', linespacing=1.7)


def arrow(p, q):
    axa.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>,head_length=3.0,head_width=1.7', lw=0.7,
                                  color='#555555', shrinkA=0, shrinkB=0))


# colour = meaning, identical to panel b: grey = classical preparation / generic circuit,
# vermillion = UCCSD pool, sky blue = point-group screening, blue = spin-flip tying (GVA)
step(xs[0], ytop, 1, 'Molecular geometry', 'RHF and MP2 natural\norbitals', C['PREP'])
step(xs[1], ytop, 2, 'Active-space Hamiltonian',
     r'$h_{pq}$,  $(pq|rs)$' + '\n' + r'orbital irreps $\Gamma_p$', C['PREP'])
step(xs[2], ytop, 3, 'UCCSD pool',
     r'$\hat{G} = \hat{T} - \hat{T}^{\dagger}$' + '\nsingles + doubles', C['UCCSD'])
step(xs[2], ybot, 4, 'Point-group screening',
     r'keep $\hat{G}$ if' + '\n' + r'$\Gamma_a \otimes \Gamma_b \otimes \Gamma_i \otimes \Gamma_j \ni A_1$',
     C['PG'], lw=1.2)
step(xs[1], ybot, 5, 'Spin-flip tying',
     r'$\theta(\alpha\alpha) = \theta(\beta\beta)$' + '\n' + r'$\theta(a_\alpha b_\beta) = \theta(b_\alpha a_\beta)$',
     C['GVA'], lw=1.2)
step(xs[0], ybot, 6, 'Layered product, L layers',
     r'$U(\theta) = \Pi_{l,k}\ \exp(\theta_{lk}\hat{G}_k)$' + '\nenergy and adjoint gradient', C['PREP'])

# emphasise the two GVA steps: dashed bracket around boxes 4-5 with a label
pad = 1.1
axa.add_patch(FancyBboxPatch((xs[1] - pad, ybot - pad), xs[2] + BW - xs[1] + 2 * pad, BH + 2 * pad,
                             boxstyle='round,pad=0,rounding_size=1.2', fc='none', ec=C['GVA'],
                             lw=0.6, ls=(0, (3, 2))))
axa.text(xs[1] - pad + 1.0, ybot + BH + pad + 0.25, 'GVA', color=C['GVA'],
         fontsize=6.2, fontweight='bold', va='bottom', ha='left')

gap = 0.8
ym_t, ym_b = ytop + (BH - 3.6) / 2, ybot + (BH - 3.6) / 2
arrow((xs[0] + BW + gap, ym_t), (xs[1] - gap, ym_t))
arrow((xs[1] + BW + gap, ym_t), (xs[2] - gap, ym_t))
arrow((xs[2] + BW / 2, ytop - gap), (xs[2] + BW / 2, ybot + BH + pad + 0.3))
arrow((xs[2] - gap, ym_b), (xs[1] + BW + gap, ym_b))
arrow((xs[1] - pad - 0.2, ym_b), (xs[0] + BW + gap, ym_b))
panel_label(fig, 0.0, 0.99, 'a')

# ---------------- panel b: parameters and CNOT gates per layer, side by side ----------------
labels = ['UCCSD', '+ point group', '+ spin tying (GVA)']
colors = [v[3] for v in vals]
y = list(range(len(vals)))[::-1]
specs = [([v[1] for v in vals], 'Parameters\nper layer', 1.0),
         ([v[2] / 1000 for v in vals], 'CNOT gates\nper layer (' + r'$\times10^3$' + ')', 1000.0)]
B0, LBL, W, G = 0.655, 0.10, 0.088, 0.052
x0 = B0 + LBL
ybot_b, ytop_b = 0.22, (21.0 + 12.6) / 38
axes = [fig.add_axes([x0, ybot_b, W, ytop_b - ybot_b]), fig.add_axes([x0 + W + G, ybot_b, W, ytop_b - ybot_b])]
for i, (ax, (data, xl, scale)) in enumerate(zip(axes, specs)):
    ax.barh(y, data, color=colors, height=0.56, edgecolor='none')
    for yy, v in zip(y, data):
        ax.text(v + max(data) * 0.04, yy, f'{v * scale:,.0f}', va='center', ha='left', fontsize=6.2,
                color='#222222')
    ax.set_xlim(0, max(data) * 1.75)
    ax.set_yticks(y)
    ax.set_yticklabels(labels if i == 0 else [], fontsize=6.5)
    ax.tick_params(axis='y', length=0, pad=3)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel(xl, linespacing=1.3)
    ax.xaxis.grid(True, lw=0.35, color=C['LIGHT'])
    ax.set_axisbelow(True)
    ax.locator_params(axis='x', nbins=3)
fig.canvas.draw()
bb = axes[0].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
panel_label(fig, bb.x0 - 0.012, 0.99, 'b')
xc = (axes[0].get_position().x0 + axes[1].get_position().x1) / 2
fig.text(xc, 0.975, r'Cisplatin CAS(10e,10o), C$_{2v}$', fontsize=6.8, ha='center', va='top')
save(fig, OUT, 'Fig1')
print('values used (label, params, CNOT):', [(v[0], v[1], v[2]) for v in vals])
