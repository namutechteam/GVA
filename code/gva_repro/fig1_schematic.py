"""Figure 1: GVA construction and parameter funnel for cisplatin CAS(10e,10o)."""
import json

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from gva import style

style.apply()
cis = json.load(open("results/cisplatin_eq.json"))
cnt = cis["counts"]

fig = plt.figure(figsize=(style.DOUBLE, 62 * style.MM))
axa = fig.add_axes([0.0, 0.0, 0.62, 1.0])
axb = fig.add_axes([0.775, 0.2, 0.215, 0.66])
axa.set_xlim(0, 100)
axa.set_ylim(0, 60)
axa.axis("off")

boxes = [
    (2, 36, 28, 18, "Molecular geometry\nRHF and MP2 natural\norbitals", "#F2F2F2"),
    (35, 36, 28, 18, "Active-space Hamiltonian\n$h_{pq}$, $(pq|rs)$\norbital irreps $\\Gamma_p$", "#F2F2F2"),
    (68, 36, 30, 18, "UCCSD pool\n$\\hat G=\\hat T-\\hat T^{\\dagger}$\nsingles + doubles", "#FBE3D6"),
    (68, 6, 30, 18, "Point-group screening\nkeep $\\hat G$ if\n$\\Gamma_a\\otimes\\Gamma_b\\otimes\\Gamma_i\\otimes\\Gamma_j\\ni A_1$", "#DCEBF7"),
    (35, 6, 28, 18, "Spin-flip tying\n$\\theta(\\alpha\\alpha)=\\theta(\\beta\\beta)$\n$\\theta(a_\\alpha b_\\beta)=\\theta(b_\\alpha a_\\beta)$", "#DCEBF7"),
    (2, 6, 28, 18, "Layered product, $L$ layers\n$U(\\theta)=\\Pi_{l,k}\\,\\exp(\\theta_{lk}\\hat G_k)$\nenergy and adjoint gradient", "#DCEBF7"),
]
for x, y, w, h, txt, fc in boxes:
    axa.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.5",
                                 fc=fc, ec="#404040", lw=0.6))
    axa.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=6.3,
             linespacing=1.35)
arrow = dict(arrowstyle="-|>", lw=0.7, color="#404040", mutation_scale=7)
axa.annotate("", xy=(34.5, 45), xytext=(30.5, 45), arrowprops=arrow)
axa.annotate("", xy=(67.5, 45), xytext=(63.5, 45), arrowprops=arrow)
axa.annotate("", xy=(83, 24.8), xytext=(83, 35.2), arrowprops=arrow)
axa.annotate("", xy=(63.5, 15), xytext=(67.5, 15), arrowprops=arrow)
axa.annotate("", xy=(30.5, 15), xytext=(34.5, 15), arrowprops=arrow)
axa.text(0, 58, "a", fontsize=8, fontweight="bold", va="top")

labels = ["UCCSD", "+ point group", "+ spin tying (GVA)"]
keys = ["UCCSD", "PG only", "GVA"]
vals = [cnt[k]["n_params"] for k in keys]
cn = [cnt[k]["cnot_1L"] for k in keys]
cols = [style.C["UCCSD"], style.C["PG"], style.C["GVA"]]
y = [2, 1, 0]
axb.barh(y, vals, color=cols, height=0.62)
for yi, v, c in zip(y, vals, cn):
    axb.text(v + 15, yi, f"{v} param.\n{c:,} CNOT", va="center", fontsize=6)
axb.set_yticks(y)
axb.set_yticklabels(labels)
axb.set_xlim(0, 1500)
axb.set_xticks([0, 500, 1000])
axb.set_xlabel("Variational parameters per layer")
axb.set_title("Cisplatin CAS(10e,10o), $C_{2v}$", fontsize=6.5, pad=4)
style.panel(axb, "b", x=-0.78, y=1.12)
fig.savefig("figures/Fig1.png", dpi=600, bbox_inches="tight")
fig.savefig("figures/Fig1.pdf", bbox_inches="tight")
