"""Shared figure style for the GVA npj QI figures (Nature-family conventions).

- Matplotlib-bundled DejaVu Sans (portable release default), 6-7 pt text, 0.5 pt axes
- bold lower-case panel labels (a, b, ...) as used in the manuscript text ("Fig. 1a")
- Okabe-Ito colour-blind-safe palette with one fixed colour per method across all figures
- single column 88 mm, double column 180 mm; vector PDF/SVG with editable text
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = Path(__file__).resolve().parent
for f in (HERE / 'fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(f))

MM = 1 / 25.4
COL1, COL2 = 88 * MM, 180 * MM

# one colour per method, used identically in every figure
C = {
    'UCCSD': '#D55E00',      # vermillion
    'GVA': '#0072B2',        # blue
    'PG': '#56B4E9',         # sky blue (screening only)
    'SPIN': '#CC79A7',       # reddish purple (tying only)
    'ADAPT': '#009E73',      # bluish green
    'CLASSICAL': '#4D4D4D',  # dark grey (CISD/CCSD/CCSD(T))
    'PREP': '#8C8C8C',       # grey (classical preparation steps)
    'LIGHT': '#E6E6E6',
}

# classical references: the original figures' colours (none collides with UCCSD/GVA/ADAPT),
# all drawn dash-dotted (as the captions state) with a different dash pattern each
REF = {'CISD': '#CC79A7', 'CCSD': '#999933', 'CCSD(T)': '#882255'}
REF_LS = {'CISD': (0, (6, 1.5, 1, 1.5)), 'CCSD': (0, (3, 1.5, 1, 1.5)), 'CCSD(T)': (0, (3, 1.2, 1, 1.2, 1, 1.2))}

# optional muted variant (GVA_PALETTE=muted): same hues, lower saturation, sits better beside greys
if os.environ.get('GVA_PALETTE') == 'muted':
    C.update({'UCCSD': '#B85C38', 'GVA': '#2E5F8A', 'PG': '#7EA6C8', 'SPIN': '#A7799A', 'ADAPT': '#4C8C6B'})

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'DejaVu Sans', 'mathtext.it': 'DejaVu Sans:italic',
    'mathtext.bf': 'DejaVu Sans:bold', 'mathtext.sf': 'DejaVu Sans',
    'axes.linewidth': 0.5, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.edgecolor': '#222222',
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.labelpad': 2.5,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'legend.fontsize': 6.5, 'legend.frameon': False, 'legend.handlelength': 1.6,
    'lines.linewidth': 1.0, 'lines.markersize': 3.5,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.dpi': 600, 'figure.dpi': 150,
})


def panel_label(fig_or_ax, x, y, s, **kw):
    """Bold lower-case panel label in figure coordinates (fig) or axes coordinates (ax)."""
    tr = fig_or_ax.transFigure if hasattr(fig_or_ax, 'transFigure') else fig_or_ax.transAxes
    fig_or_ax.text(x, y, s, transform=tr, fontsize=8, fontweight='bold', va='top', ha='left', **kw)


def save(fig, out_dir, name):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(out / f'{name}.{ext}', bbox_inches='tight', pad_inches=0.02)
