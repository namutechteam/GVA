import matplotlib as mpl

# Okabe-Ito colour-blind-safe palette
C = {"GVA": "#0072B2", "UCCSD": "#D55E00", "ADAPT": "#009E73",
     "PG": "#56B4E9", "spin": "#E69F00", "grey": "#7F7F7F", "black": "#000000",
     "CISD": "#CC79A7", "CCSD": "#999933", "CCSDT": "#882255", "HF": "#7F7F7F"}

MM = 1 / 25.4
SINGLE = 89 * MM
DOUBLE = 183 * MM


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["FreeSans", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "lines.linewidth": 1.0, "lines.markersize": 3.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.dpi": 600, "figure.dpi": 150,
        "mathtext.default": "regular", "pdf.fonttype": 42,
    })


def panel(ax, letter, x=-0.18, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")
