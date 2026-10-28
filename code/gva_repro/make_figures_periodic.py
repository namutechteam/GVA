"""Figures 9 and 10 and Supplementary Fig. 2 (periodic systems)."""
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from gva import style

style.apply()
C = style.C
CHEM = 1.6
FLOOR = 1e-6

import os
per = json.load(open("results/periodic.json")) if os.path.exists("results/periodic.json") \
    else {"systems": [], "runs": []}
sp = json.load(open("results/scaling_periodic.json"))
sm = json.load(open("results/small_molecules.json"))
scal = json.load(open("results/scaling.json"))
cis = json.load(open("results/cisplatin_eq.json"))
SYS = {s["system"]: s for s in per["systems"]}
H_ORDER = {"D2h": 8, "C2v": 4, "Cs": 2, "C1": 1}


def tag(r, d, nk):
    return f"H2chain_r{r:.2f}_d{d:.2f}_Nk{nk}"


def runs(system, ansatz, layers=None, seed=None):
    out = [x for x in per["runs"] if x["system"] == system and x["ansatz"] == ansatz
           and (layers is None or x["layers"] == layers) and (seed is None or x["seed"] == seed)]
    return out


def mha(x):
    return np.maximum(np.abs(np.asarray(x, float)) * 1e3, FLOOR)


def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def molecular_ratio_points():
    """(h, parameter ratio, CNOT ratio, label) for the molecular data."""
    pts = []
    for s in sm["systems"]:
        if s["scale"] != 1.0:
            continue
        c = s["counts"]
        irr = s["irreps"]
        pts.append((s["system"], c["GVA"]["n_params"] / c["UCCSD"]["n_params"],
                    c["GVA"]["cnot_1L"] / c["UCCSD"]["cnot_1L"]))
    return pts


# ------------------------------------------------------------------ Fig 9
def fig9():
    fig, axs = plt.subplots(1, 3, figsize=(style.DOUBLE, 60 * style.MM),
                            gridspec_kw={"wspace": 0.42})
    h = np.logspace(0, np.log10(160), 50)
    # crystals
    one = [r for r in sp["one_d"] if r["o"] == 1 and r["n_k"] >= 2]
    thr = [r for r in sp["three_d"] if r["o"] == 2 and r["n_k"] >= 8]
    mc = [r for r in sp["molecular_crystal"] if r["n_k"] >= 8]
    # molecules: point-group order from the small-molecule irreps and cisplatin (C2v)
    mol_h = {"H2": 8, "H4": 8, "LiH": 4, "BeH2": 8}
    mol = [(mol_h[s["system"]], s["counts"]["GVA"]["n_params"] / s["counts"]["UCCSD"]["n_params"],
            s["counts"]["GVA"]["cnot_1L"] / s["counts"]["UCCSD"]["cnot_1L"])
           for s in sm["systems"] if s["scale"] == 1.0]
    cisp = [(4, r["GVA_params"] / r["UCCSD_params"], r["GVA_cnot_1L"] / r["UCCSD_cnot_1L"])
            for r in scal if r["n_active_orbitals"] >= 10]
    for k, (ax, key_g, key_u, idx, law, lab) in enumerate((
            (axs[0], "GVA_k_params", "UCCSD_params", 1, 0.5 / h, "parameters"),
            (axs[1], "k_cnot_1L", "UCCSD_cnot_1L", 2, 1 / h, "CNOT gates"))):
        ax.plot(h, law, color=C["grey"], lw=0.8, ls="--", zorder=1)
        ax.plot([r["n_k"] for r in one], [r[key_g] / r[key_u] for r in one], "o-",
                color=C["GVA"], ms=2.5, lw=0.7, label="1D chain, 1+1 bands")
        ax.plot([r["n_k"] for r in thr], [r[key_g] / r[key_u] for r in thr], "s-",
                color=C["ADAPT"], ms=3, lw=0.7, label="3D mesh, 2+2 bands")
        ax.plot([r["n_k"] for r in mc], [r[key_g] / r[key_u] for r in mc], "D",
                color=C["CCSDT"], ms=3, mfc="none", mew=0.8, label="3D mesh, 8+8 bands")
        ax.plot([m[0] for m in mol], [m[idx] for m in mol], "^", color=C["UCCSD"], ms=3.5,
                label="small molecules")
        ax.plot([m[0] for m in cisp], [m[idx] for m in cisp], "v", color=C["spin"], ms=3.5,
                mfc="none", mew=0.8, label="cisplatin CAS")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("symmetry group order $h$ ($h = N_k$ for crystals)")
        ax.set_ylabel(f"GVA / UCCSD, {lab}")
        ax.set_xlim(1.5, 160)
        ax.text(0.97, 0.93, "$1/(2h)$" if k == 0 else "$1/h$", transform=ax.transAxes,
                ha="right", va="top", color=C["grey"])
        style.panel(ax, "ab"[k], x=-0.27)
    axs[0].legend(loc="lower left", fontsize=5.5, handlelength=1.6)
    ax = axs[2]
    rr = [r for r in sp["three_d"] if r["o"] == 2]
    q = [r["n_qubits"] for r in rr]
    ax.plot(q, [r["UCCSD_params"] for r in rr], "o-", color=C["UCCSD"], mfc="none", ms=3,
            lw=0.7, label="UCCSD parameters")
    ax.plot(q, [r["GVA_k_params"] for r in rr], "o-", color=C["GVA"], ms=3, lw=0.7,
            label="GVA-k parameters")
    ax.plot(q, [r["UCCSD_cnot_1L"] for r in rr], "s--", color=C["UCCSD"], mfc="none", ms=3,
            lw=0.7, label="UCCSD CNOT")
    ax.plot(q, [r["k_cnot_1L"] for r in rr], "s--", color=C["GVA"], ms=3, lw=0.7,
            label="GVA-k CNOT")
    for r in rr:
        if r["n_k"] > 1:
            ax.annotate(f"{r['mesh'][0]}$^3$", (r["n_qubits"], r["GVA_k_params"]),
                        textcoords="offset points", xytext=(2, -7), fontsize=5.5,
                        color=C["grey"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("qubits (3D mesh, 2+2 bands per cell)")
    ax.set_ylabel("count per layer")
    ax.legend(loc="upper left", fontsize=5.5)
    style.panel(ax, "c", x=-0.27)
    save(fig, "Fig9")


# ------------------------------------------------------------------ Fig 10
def best_med(system, ansatz, layers):
    # single-layer ansaetze and ADAPT-VQE: Hartree-Fock start (seed 0), as in Table 3
    rr = runs(system, ansatz, layers, 0 if layers in (1, None) else None)
    e = np.array([x["dE"] for x in rr])
    return np.min(e), np.median(e), np.max(e), rr


def fig10():
    fig = plt.figure(figsize=(style.DOUBLE, 62 * style.MM))
    gs = fig.add_gridspec(1, 3, wspace=0.42)
    spec = [("UCCSD", 1, C["UCCSD"], "o", "none", "UCCSD 1L"),
            ("GVA-k", 1, C["GVA"], "o", "full", "GVA-k 1L"),
            ("UCCSD", 2, C["UCCSD"], "s", "none", "UCCSD 2L"),
            ("GVA-k", 2, C["GVA"], "s", "full", "GVA-k 2L"),
            ("ADAPT-VQE", None, C["ADAPT"], "^", "full", "ADAPT (GVA-k pool)")]
    # a: dimerized chain vs N_k
    ax = fig.add_subplot(gs[0, 0])
    nks = [2, 3, 4, 5]
    offs = np.linspace(-0.2, 0.2, len(spec))
    for (an, L, col, mk, fill, lab), dx in zip(spec, offs):
        lo, me, hi = [], [], []
        for nk in nks:
            a, b, c, _ = best_med(tag(0.74, 1.50, nk), an, L)
            lo.append(a); me.append(b); hi.append(c)
        x = np.array(nks) + dx
        ax.errorbar(x, mha(me), yerr=[mha(me) - mha(lo), mha(hi) - mha(me)], fmt=mk,
                    color=col, mfc=col if fill == "full" else "white", ms=3.2, lw=0.6,
                    capsize=1.2, label=lab)
    for key, col, lab in (("E_CCSD", C["CCSD"], "CCSD"), ("E_CCSD(T)", C["CCSDT"], "CCSD(T)")):
        ax.plot(nks, mha([SYS[tag(0.74, 1.50, nk)][key] - SYS[tag(0.74, 1.50, nk)]["E_FCI"]
                          for nk in nks]), "-", color=col, lw=0.8, label=lab)
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.7)
    ax.set_yscale("log")
    ax.set_xticks(nks)
    ax.set_xlabel("$N_k$ (dimerized chain)")
    ax.set_ylabel("|E − E$_{FCI}$| (mHa)")
    style.panel(ax, "a", x=-0.27)
    # b: uniform chains vs R, N_k = 5 (filled) and 3 (open, shifted)
    ax = fig.add_subplot(gs[0, 1])
    Rs = [1.0, 1.5, 2.0]
    for (an, L, col, mk, fill, lab), dx in zip(spec, offs):
        for nk, mfc_on in ((5, True), (3, False)):
            lo, me, hi = [], [], []
            for R in Rs:
                a, b, c, _ = best_med(tag(R, R, nk), an, L)
                lo.append(a); me.append(b); hi.append(c)
            x = np.arange(len(Rs)) + dx + (0.0 if nk == 5 else 0.04)
            ax.errorbar(x, mha(me), yerr=[mha(me) - mha(lo), mha(hi) - mha(me)], fmt=mk,
                        color=col, mfc=(col if fill == "full" else "white") if mfc_on else "white",
                        alpha=1.0 if nk == 5 else 0.45, ms=3.2, lw=0.6, capsize=1.2)
    for key, col in (("E_CCSD", C["CCSD"]), ("E_CCSD(T)", C["CCSDT"])):
        ax.plot(range(len(Rs)), mha([SYS[tag(R, R, 5)][key] - SYS[tag(R, R, 5)]["E_FCI"]
                                     for R in Rs]), "-", color=col, lw=0.8)
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.7)
    ax.set_yscale("log")
    ax.set_xticks(range(len(Rs)))
    ax.set_xticklabels([f"{R:.1f}" for R in Rs])
    ax.set_xlabel("H–H spacing (Å), uniform chain")
    ax.set_xlim(-0.4, 2.4)
    ax.set_ylabel("|E − E$_{FCI}$| (mHa)")
    for i, R in enumerate(Rs):
        ax.text(i, 1.02, f"$c_0^2$ = {SYS[tag(R, R, 5)]['c0_sq']:.2f}",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.5,
                color=C["grey"])
    style.panel(ax, "b", x=-0.27)
    # c: error vs CNOT for all N_k = 5 systems
    ax = fig.add_subplot(gs[0, 2])
    for (an, L, col, mk, fill, lab) in spec:
        for t, shade in ((tag(0.74, 1.50, 5), 1.0), (tag(1.0, 1.0, 5), 0.75),
                         (tag(1.5, 1.5, 5), 0.5), (tag(2.0, 2.0, 5), 0.3)):
            rr = runs(t, an, L, 0 if L in (1, None) else None)
            if not rr:
                continue
            e = min(x["dE"] for x in rr)
            cn = rr[0]["cnot"]
            ax.plot(cn, mha(e), mk, color=col, mfc=col if fill == "full" else "white",
                    alpha=shade, ms=3.5)
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    ax.xaxis.set_major_locator(FixedLocator([1e4, 3e4, 1e5, 3e5]))
    ax.xaxis.set_major_formatter(FixedFormatter(["10$^4$", "3×10$^4$", "10$^5$", "3×10$^5$"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(1e4, 3.5e5)
    ax.set_xlabel("CNOT gates (20 qubits, $N_k$ = 5)")
    ax.set_ylabel("best |E − E$_{FCI}$| (mHa)")
    style.panel(ax, "c", x=-0.27)
    handles = [Line2D([], [], marker=mk, color=col, mfc=col if fill == "full" else "white",
                      ls="none", ms=3.5, label=lab) for (an, L, col, mk, fill, lab) in spec]
    handles += [Line2D([], [], color=C["CCSD"], label="CCSD"),
                Line2D([], [], color=C["CCSDT"], label="CCSD(T)"),
                Line2D([], [], color=C["grey"], ls=":", label="1.6 mHa")]
    fig.legend(handles=handles, loc="lower center", ncol=8, bbox_to_anchor=(0.5, -0.13),
               fontsize=5.8, handletextpad=0.3, columnspacing=0.9)
    save(fig, "Fig10")


# ------------------------------------------------------------------ Fig S2
def figS2():
    """Tying ablation (1L, HF start) and random starts for the untied k-only ansatz."""
    systems = [s["system"] for s in per["systems"]]
    lab = [f"{SYS[t]['r']:.2f}/{SYS[t]['d']:.2f}\n$N_k$={SYS[t]['n_k']}" for t in systems]
    fig, axs = plt.subplots(2, 1, figsize=(style.DOUBLE, 95 * style.MM), sharex=True,
                            gridspec_kw={"hspace": 0.15})
    x = np.arange(len(systems))
    var = [("UCCSD", C["UCCSD"], "o", "none"), ("k only", C["PG"], "o", "full"),
           ("k + spin", C["spin"], "D", "full"), ("k + TR", C["CISD"], "v", "full"),
           ("GVA-k", C["GVA"], "s", "full")]
    offs = np.linspace(-0.28, 0.28, len(var))
    for (an, col, mk, fill), dx in zip(var, offs):
        e = [(runs(t, an, 1, 0)[0]["E"] - runs(t, "k only", 1, 0)[0]["E"]) * 1e3 for t in systems]
        axs[0].plot(x + dx, e, mk, color=col, mfc=col if fill == "full" else "white",
                    ms=3.2, ls="none", label=an)
        s2 = []
        for t in systems:
            rr = runs(t, an, 1)
            s2.append(max(r["S2"] for r in rr))
        axs[1].plot(x + dx, np.maximum(s2, 1e-16), mk, color=col,
                    mfc=col if fill == "full" else "white", ms=3.2, ls="none")
    axs[0].axhline(0, color=C["grey"], ls="-", lw=0.5)
    axs[0].set_yscale("symlog", linthresh=1e-3, linscale=0.6)
    axs[0].set_ylim(-3, 30)
    axs[0].set_ylabel("E − E(k only) (mHa), 1L, HF start")
    axs[0].legend(ncol=5, loc="lower left", fontsize=5.8)
    axs[1].set_yscale("log")
    axs[1].set_ylabel("largest ⟨$\\hat S^2$⟩ over starts, 1L")
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(lab, fontsize=5.5)
    axs[1].set_xlabel("$r$/$d$ (Å)")
    style.panel(axs[0], "a", x=-0.07)
    style.panel(axs[1], "b", x=-0.07)
    save(fig, "FigS2")


if __name__ == "__main__":
    fig9()
    fig10()
    figS2()
    print("periodic figures written")
