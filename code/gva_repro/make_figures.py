"""Figures 2-8 from the computed results (results/*.json)."""
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from gva import style
from gva.noise import bias_unmitigated, budget

style.apply()
C = style.C
CHEM = 1.6          # mHa
FLOOR = 1e-6        # mHa, plotting floor for converged runs
P_HW = 1.84e-3      # average two-qubit infidelity, Moses et al. PRX 13, 041052 (2023)

sm = json.load(open("results/small_molecules.json"))
cis = json.load(open("results/cisplatin_eq.json"))
meas = json.load(open("results/measurement.json"))
scal = json.load(open("results/scaling.json"))
pes = json.load(open("results/pes.json"))

MOLS = ["H2", "H4", "LiH", "BeH2"]
MOLTEX = {"H2": "H$_2$", "H4": "H$_4$", "LiH": "LiH", "BeH2": "BeH$_2$"}
SCALES = [1.0, 1.5, 2.0]


def runs(**kw):
    out = []
    for r in sm["runs"]:
        if all(r.get(k) == v for k, v in kw.items()):
            out.append(r)
    return out


def mha(x):
    return np.maximum(np.abs(np.asarray(x)) * 1e3, FLOOR)


def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ Fig 2
def fig2():
    fig = plt.figure(figsize=(style.DOUBLE, 112 * style.MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.55, wspace=0.42)
    ax = fig.add_subplot(gs[0, :])
    xs, labels = [], []
    x = 0
    spec = [("UCCSD", 1, "o", "none"), ("GVA", 1, "o", "full"),
            ("UCCSD", 2, "s", "none"), ("GVA", 2, "s", "full")]
    offs = [-0.27, -0.09, 0.09, 0.27]
    for m in MOLS:
        for sc in SCALES:
            for (an, L, mk, fill), dx in zip(spec, offs):
                rr = runs(system=m, scale=sc, ansatz=an, layers=L)
                e = mha([r["dE"] for r in rr])
                med = np.median(e)
                col = C[an]
                ax.plot([x + dx] * 2, [e.min(), e.max()], color=col, lw=0.6)
                ax.plot(x + dx, med, mk, ms=3.6, mec=col, mew=0.8,
                        mfc=col if fill == "full" else "white")
            ad = runs(system=m, scale=sc, ansatz="ADAPT-VQE")[0]
            ax.plot(x + 0.42, mha(ad["dE"]), "^", ms=3.4, color=C["ADAPT"])
            xs.append(x)
            labels.append(f"{MOLTEX[m]}\n×{sc:.1f}")
            x += 1
        x += 0.5
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
    ax.text(xs[-1] + 0.65, CHEM * 1.3, "1.6 mHa", fontsize=6, color=C["grey"],
            ha="right")
    ax.set_yscale("log")
    ax.set_ylim(5e-7, 200)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.6, xs[-1] + 0.7)
    ax.set_ylabel("|ΔE| relative to FCI (mHa)")
    handles = [Line2D([], [], marker=mk, ls="", mec=C[an], mfc=C[an] if f == "full" else "white",
                      mew=0.8, label=f"{an} {L}L") for an, L, mk, f in spec]
    handles.append(Line2D([], [], marker="^", ls="", color=C["ADAPT"], label="ADAPT-VQE"))
    ax.legend(handles=handles, ncol=5, loc="upper left", bbox_to_anchor=(0, 1.13))
    style.panel(ax, "a", x=-0.07, y=1.06)

    # b: parameters and CNOTs per layer
    axb = fig.add_subplot(gs[1, 0])
    sysd = {s["system"]: s for s in sm["systems"] if s["scale"] == 1.0}
    w = 0.38
    xx = np.arange(len(MOLS))
    pu = [sysd[m]["counts"]["UCCSD"]["n_params"] for m in MOLS]
    pg = [sysd[m]["counts"]["GVA"]["n_params"] for m in MOLS]
    axb.bar(xx - w / 2, pu, w, color=C["UCCSD"], label="UCCSD")
    axb.bar(xx + w / 2, pg, w, color=C["GVA"], label="GVA")
    for i in range(len(MOLS)):
        axb.text(xx[i] + w / 2, pg[i] * 1.15, str(pg[i]), ha="center", fontsize=5.5)
        axb.text(xx[i] - w / 2, pu[i] * 1.15, str(pu[i]), ha="center", fontsize=5.5)
    axb.set_yscale("log")
    axb.set_ylim(0.6, 800)
    axb.set_xticks(xx)
    axb.set_xticklabels([MOLTEX[m] for m in MOLS])
    axb.set_ylabel("Parameters per layer")
    axb.legend(loc="upper left")
    style.panel(axb, "b", x=-0.32)

    axc = fig.add_subplot(gs[1, 1])
    cu = [sysd[m]["counts"]["UCCSD"]["cnot_1L"] for m in MOLS]
    cg = [sysd[m]["counts"]["GVA"]["cnot_1L"] for m in MOLS]
    axc.bar(xx - w / 2, cu, w, color=C["UCCSD"])
    axc.bar(xx + w / 2, cg, w, color=C["GVA"])
    for i in range(len(MOLS)):
        red = 100 * (1 - cg[i] / cu[i])
        axc.text(xx[i], max(cu[i], cg[i]) * 1.5, f"−{red:.0f}%", ha="center",
                 fontsize=5.5)
    axc.set_yscale("log")
    axc.set_ylim(20, 2e5)
    axc.set_xticks(xx)
    axc.set_xticklabels([MOLTEX[m] for m in MOLS])
    axc.set_ylabel("CNOT gates per layer")
    style.panel(axc, "c", x=-0.32)

    # d: convergence BeH2 x2.0
    axd = fig.add_subplot(gs[1, 2])
    ef = [s for s in sm["systems"] if s["system"] == "BeH2" and s["scale"] == 2.0][0]["E_FCI"]
    for an, L, ls in (("UCCSD", 1, "--"), ("GVA", 1, "-"), ("UCCSD", 2, "--"), ("GVA", 2, "-")):
        r = runs(system="BeH2", scale=2.0, ansatz=an, layers=L, seed=0)[0]
        h = mha(np.array(r["history"]) - ef)
        axd.plot(np.arange(len(h)), h, ls, color=C[an], lw=0.9 if L == 1 else 1.4,
                 label=f"{an} {L}L")
    axd.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
    axd.set_yscale("log")
    axd.set_xscale("symlog", linthresh=10)
    axd.set_ylim(5e-7, 500)
    axd.set_xlim(left=0)
    axd.set_xlabel("Optimizer iteration")
    axd.set_ylabel("|ΔE| (mHa)")
    axd.set_title("BeH$_2$, bonds ×2.0", fontsize=6.5)
    axd.legend(loc="lower left", fontsize=5.5)
    style.panel(axd, "d", x=-0.32)
    save(fig, "Fig2")


# ------------------------------------------------------------------ Fig 3
def best(rec, kind, key):
    vals = rec[kind]
    m = min(vals, key=lambda k: vals[k][key])
    return vals[m][key], m


def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(style.DOUBLE * 0.82, 60 * style.MM))
    labels = ["H$_2$", "H$_4$", "LiH", "BeH$_2$", "Cisplatin\nCAS(10e,10o)"]
    xx = np.arange(len(meas))
    w = 0.38
    q = [best(r, "QWC", "n_groups")[0] for r in meas]
    g = [best(r, "GC", "n_groups")[0] for r in meas]
    ax = axs[0]
    ax.bar(xx - w / 2, q, w, color=C["grey"], label="Qubit-wise commuting")
    ax.bar(xx + w / 2, g, w, color=C["GVA"], label="General commuting")
    for i in range(len(xx)):
        ax.text(xx[i] - w / 2, q[i] * 1.12, str(q[i]), ha="center", fontsize=5.5)
        ax.text(xx[i] + w / 2, g[i] * 1.12, str(g[i]), ha="center", fontsize=5.5)
    ax.set_yscale("log")
    ax.set_ylim(1, 1e4)
    ax.set_xticks(xx)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Measurement groups")
    ax.legend(loc="upper left")
    style.panel(ax, "a", x=-0.2)
    ax = axs[1]
    u = [r["shots_ungrouped"] for r in meas]
    qs = [best(r, "QWC", "shots")[0] for r in meas]
    gs_ = [best(r, "GC", "shots")[0] for r in meas]
    w = 0.27
    ax.bar(xx - w, u, w, color="#BFBFBF", label="No grouping")
    ax.bar(xx, qs, w, color=C["grey"], label="QWC")
    ax.bar(xx + w, gs_, w, color=C["GVA"], label="GC")
    for i in range(len(xx)):
        ax.text(xx[i] + w, gs_[i] * 1.25, f"{gs_[i] / qs[i]:.2f}×", ha="center",
                fontsize=5.5)
    ax.set_yscale("log")
    ax.set_xticks(xx)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Shots per energy evaluation\n(σ = 1.6 mHa)")
    ax.legend(loc="upper left")
    style.panel(ax, "b", x=-0.22)
    fig.tight_layout(w_pad=2.5)
    save(fig, "Fig3")


# ------------------------------------------------------------------ Fig 4
def cis_run(an, L, seed=0):
    for r in cis["runs"]:
        if r["ansatz"] == an and r["layers"] == L and r["seed"] == seed:
            return r
    return None


def fig4():
    fig = plt.figure(figsize=(style.DOUBLE, 62 * style.MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.25, 0.95], wspace=0.62)
    ax = fig.add_subplot(gs[0])
    occ = np.array([o["occupation"] for o in
                    json.load(open("results/cisplatin_orbitals.json"))["orbitals"]])
    dev = np.where(occ > 1, 2 - occ, occ)
    cols = [C["UCCSD"] if o > 1 else C["GVA"] for o in occ]
    ax.bar(np.arange(10), dev, color=cols)
    ax.set_yscale("log")
    ax.set_ylim(5e-3, 0.2)
    ax.set_xticks(np.arange(10))
    orb = json.load(open("results/cisplatin_orbitals.json"))["orbitals"]
    ax.set_xticklabels([o["irrep"].replace("1", "$_1$").replace("2", "$_2$") for o in orb],
                       fontsize=5.5)
    ax.set_xlabel("Active MP2 natural orbital")
    ax.set_ylabel("2 − n (occupied), n (virtual)")
    style.panel(ax, "a", x=-0.36)

    ax = fig.add_subplot(gs[1])
    ef = cis["E_CASCI"]
    for an, ls in (("UCCSD", "--"), ("GVA", "-")):
        for L, lw in ((1, 0.8), (2, 1.2), (3, 1.6)):
            r = cis_run(an, L)
            if r is None:
                continue
            h = mha(np.array(r["history"]) - ef)
            ax.plot(np.arange(len(h)), h, ls, color=C[an], lw=lw, label=f"{an} {L}L")
    for key, col, lab in (("E_CISD", C["CISD"], "CISD"), ("E_CCSD", C["CCSD"], "CCSD"),
                          ("E_CCSD(T)", C["CCSDT"], "CCSD(T)")):
        v = (cis[key] - ef) * 1e3
        ax.axhline(abs(v), color=col, lw=0.7, ls="-.")
        ax.text(1.02, abs(v), lab, color=col, fontsize=5.5, va="center",
                transform=ax.get_yaxis_transform())
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_xscale("symlog", linthresh=10)
    ax.set_ylim(5e-3, 60)
    ax.set_xlim(left=0)
    ax.set_xlabel("Optimizer iteration")
    ax.set_ylabel("|ΔE| relative to CASCI (mHa)")
    ax.legend(ncol=2, loc="lower left", fontsize=5.5)
    style.panel(ax, "b", x=-0.22)

    ax = fig.add_subplot(gs[2])
    ad = [r for r in cis["runs"] if r["ansatz"] == "ADAPT-VQE"][0]
    n = [t["n_ops"] for t in ad["trace"]]
    e = mha([t["E"] - ef for t in ad["trace"]])
    ax.plot(n, e, "-", color=C["ADAPT"], lw=1.0)
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
    ax.axhline(cis_run("GVA", 1)["dE"] * 1e3, color=C["GVA"], ls="--", lw=0.7)
    ax.text(n[-1], cis_run("GVA", 1)["dE"] * 1e3 * 0.75, "GVA 1L", color=C["GVA"],
            fontsize=5.5, ha="right", va="top")
    ax.set_yscale("log")
    ax.set_ylim(5e-3, 60)
    ax.set_xlabel("ADAPT operators added")
    ax.set_ylabel("|ΔE| (mHa)")
    style.panel(ax, "c", x=-0.3)
    save(fig, "Fig4")


# ------------------------------------------------------------------ Fig 5
def fig5():
    fig, axs = plt.subplots(1, 2, figsize=(style.DOUBLE * 0.82, 62 * style.MM))
    for ax, key, xl in ((axs[0], "cnot", "CNOT gates (ladder compilation)"),
                        (axs[1], "n_params", "Variational parameters")):
        for an, mk in (("UCCSD", "o"), ("GVA", "o"), ("PG only", "D"), ("spin only", "D")):
            for L in (1, 2, 3):
                rr = [r for r in cis["runs"] if r["ansatz"] == an and r["layers"] == L]
                if not rr:
                    continue
                e = mha([r["dE"] for r in rr])
                col = C.get(an, C["PG"] if an == "PG only" else C["spin"])
                abl = an in ("PG only", "spin only")
                spread = e.max() - e.min() > 1e-4
                ax.errorbar(rr[0][key], np.median(e),
                            yerr=[[np.median(e) - e.min()], [e.max() - np.median(e)]] if spread else None,
                            fmt=mk, ms=9 if abl else 3.5 + L, color=col,
                            mfc="none" if (abl or an == "UCCSD") else col,
                            mew=0.9, capsize=1.5, lw=0.6, zorder=2 if abl else 3)
                if an in ("UCCSD", "GVA"):
                    ax.annotate(f"{L}L", (rr[0][key], np.median(e)), textcoords="offset points",
                                xytext=(5, 3), fontsize=5.5, color=col)
        ad = [r for r in cis["runs"] if r["ansatz"] == "ADAPT-VQE"][0]
        ax.plot(ad[key], mha(ad["dE"]), "^", color=C["ADAPT"], ms=5)
        ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
        for k2, col in (("E_CCSD", C["CCSD"]), ("E_CCSD(T)", C["CCSDT"])):
            ax.axhline(abs(cis[k2] - cis["E_CASCI"]) * 1e3, color=col, lw=0.6, ls="-.")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel("|ΔE| relative to CASCI (mHa)")
        ax.minorticks_off()
        if key == "cnot":
            ax.set_xticks([2e4, 5e4, 1e5, 2e5, 4e5])
            ax.set_xticklabels(["2×10$^4$", "5×10$^4$", "10$^5$", "2×10$^5$", "4×10$^5$"])
        else:
            ax.set_xticks([100, 300, 1000, 3000])
            ax.set_xticklabels(["100", "300", "1,000", "3,000"])
        ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1])
        ax.set_yticklabels(["0.01", "0.03", "0.1", "0.3", "1"])
    handles = [Line2D([], [], marker="o", ls="", mec=C["UCCSD"], mfc="white", label="UCCSD"),
               Line2D([], [], marker="o", ls="", color=C["GVA"], label="GVA"),
               Line2D([], [], marker="D", ls="", mec=C["PG"], mfc="none", ms=6, label="Screening only 1L"),
               Line2D([], [], marker="D", ls="", mec=C["spin"], mfc="none", ms=6, label="Tying only 1L"),
               Line2D([], [], marker="^", ls="", color=C["ADAPT"], label="ADAPT-VQE"),
               Line2D([], [], ls="-.", color=C["CCSD"], label="CCSD"),
               Line2D([], [], ls="-.", color=C["CCSDT"], label="CCSD(T)")]
    style.panel(axs[0], "a", x=-0.22)
    style.panel(axs[1], "b", x=-0.22)
    fig.tight_layout(w_pad=2.5, rect=(0, 0.1, 1, 1))
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=5.5,
               bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=1.0)
    save(fig, "Fig5")


# ------------------------------------------------------------------ Fig 6
def fit_slope(x, y):
    p = np.polyfit(np.log(x), np.log(y), 1)
    return p[0]


def fig6():
    fig, axs = plt.subplots(1, 2, figsize=(style.DOUBLE * 0.75, 60 * style.MM))
    n = np.array([r["n_active_orbitals"] for r in scal])
    for ax, ku, kg, yl in ((axs[0], "UCCSD_params", "GVA_params", "Parameters per layer"),
                           (axs[1], "UCCSD_cnot_1L", "GVA_cnot_1L", "CNOT gates per layer")):
        yu = np.array([r[ku] for r in scal], float)
        yg = np.array([r[kg] for r in scal], float)
        su, sg = fit_slope(n[3:], yu[3:]), fit_slope(n[3:], yg[3:])
        ax.plot(n, yu, "o-", color=C["UCCSD"], mfc="white", label=f"UCCSD (slope {su:.2f})")
        ax.plot(n, yg, "o-", color=C["GVA"], label=f"GVA (slope {sg:.2f})")
        ax.set_xscale("log")
        ax.set_yscale("log")
        tk = [4, 6, 8, 10, 12, 16, 20]
        ax.set_xticks(tk)
        ax.set_xticklabels([str(v) for v in tk])
        ax.minorticks_off()
        ax.set_xlabel("Active orbitals (qubits = 2×)")
        ax.set_ylabel(yl)
        ax.legend(loc="lower right")
        sec = ax.twinx()
        sec.plot(n, yg / yu, "s:", color=C["grey"], ms=2.5, lw=0.7)
        sec.set_ylim(0, 0.6)
        sec.set_ylabel("GVA / UCCSD", color=C["grey"])
        sec.spines["right"].set_visible(True)
        sec.tick_params(axis="y", colors=C["grey"])
    style.panel(axs[0], "a", x=-0.25)
    style.panel(axs[1], "b", x=-0.25)
    fig.tight_layout(w_pad=3)
    save(fig, "Fig6")


# ------------------------------------------------------------------ Fig 7
def noise_table():
    gap_cis = [m for m in meas if m["system"].startswith("Cisplatin")][0]
    rows = []
    for m in meas:
        name = m["system"]
        gap = m["identity_coeff"] - m["E_exact"]
        if name.startswith("Cisplatin"):
            circuits = [(an, L, cis_run(an, L)["cnot"], cis_run(an, L)["dE"])
                        for an in ("GVA", "UCCSD") for L in (1, 2, 3)
                        if cis_run(an, L) is not None]
        else:
            circuits = []
            for an in ("GVA", "UCCSD"):
                for L in (1, 2):
                    r = runs(system=name, scale=1.0, ansatz=an, layers=L, seed=0)[0]
                    circuits.append((an, L, r["cnot"], r["dE"]))
        for an, L, nc, de in circuits:
            b = budget(nc, gap)
            rows.append({"system": name, "ansatz": an, "layers": L, "cnot": nc,
                         "dE_noiseless_mHa": de * 1e3, "gap_Ha": gap, **b})
    return rows


def fig7(rows):
    fig, axs = plt.subplots(1, 2, figsize=(style.DOUBLE * 0.82, 62 * style.MM))
    ax = axs[0]
    cm = [m for m in meas if m["system"].startswith("Cisplatin")][0]
    gap = cm["identity_coeff"] - cm["E_exact"]
    p = np.logspace(-10, -2, 400)
    for an, L, ls in (("GVA", 1, "-"), ("GVA", 2, "-"), ("UCCSD", 1, "--"), ("UCCSD", 2, "--")):
        nc = cis_run(an, L)["cnot"]
        ax.plot(p, bias_unmitigated(p, nc, gap) * 1e3, ls, color=C[an],
                lw=0.8 if L == 1 else 1.4, label=f"{an} {L}L ({nc/1e3:.0f}k CNOT)")
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
    ax.axvline(P_HW, color=C["grey"], ls="--", lw=0.7)
    ax.text(P_HW * 0.8, 2e-3, "1.8×10$^{-3}$", ha="right", fontsize=5.5, color=C["grey"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 2e5)
    ax.set_xlabel("Two-qubit depolarizing probability p")
    ax.set_ylabel("Energy bias, unmitigated (mHa)")
    ax.set_title("Cisplatin CAS(10e,10o)", fontsize=6.5)
    ax.legend(loc="upper left", fontsize=5.5)
    style.panel(ax, "a", x=-0.22)

    ax = axs[1]
    sel = []
    for name in [m["system"] for m in meas]:
        rr = [r for r in rows if r["system"] == name and r["ansatz"] == "GVA"
              and abs(r["dE_noiseless_mHa"]) < CHEM]
        rr.sort(key=lambda r: r["cnot"])
        sel.append(rr[0])
    xx = np.arange(len(sel))
    w = 0.27
    for k, (key, col, lab) in enumerate((("p_max_unmitigated", C["UCCSD"], "Unmitigated, bias ≤ 1.6 mHa"),
                                         ("p_max_richardson", C["spin"], "Richardson ZNE (λ = 1, 3, 5)"),
                                         ("p_max_rescaling_overhead", C["GVA"], "Rescaling, overhead ≤ 100"))):
        ax.bar(xx + (k - 1) * w, [r[key] for r in sel], w, color=col, label=lab)
    ax.set_yscale("log")
    ax.set_ylim(1e-10, 10)
    ax.set_xticks(xx)
    ax.set_xticklabels([f"{t}\n{r['layers']}L, {r['cnot']:,}" for t, r in
                        zip(["H$_2$", "H$_4$", "LiH", "BeH$_2$", "Cisplatin"], sel)],
                       fontsize=5.5)
    ax.set_ylabel("Largest tolerable p")
    ax.axhline(P_HW, color=C["grey"], ls="--", lw=0.7)
    ax.legend(loc="upper right", fontsize=5.5)
    style.panel(ax, "b", x=-0.22)
    fig.tight_layout(w_pad=2.5)
    save(fig, "Fig7")


# ------------------------------------------------------------------ Fig 8
def fig8():
    R = np.array([r["R_PtCl"] for r in pes])
    e0 = min(r["E_CASCI"] for r in pes)
    fig, axs = plt.subplots(1, 3, figsize=(style.DOUBLE, 60 * style.MM))
    ax = axs[0]
    for key, col, lab, ls in (("E_ref", C["HF"], "RHF", "-"), ("E_CCSD", C["CCSD"], "CCSD", "-"),
                              ("E_CCSD(T)", C["CCSDT"], "CCSD(T)", "-"),
                              ("E_CASCI", C["black"], "CASCI", "-")):
        ax.plot(R, [(r[key] - e0) * 627.509 for r in pes], ls, color=col, label=lab,
                marker="o", ms=2)
    ax.plot(R, [(r["GVA_2L"]["E"] - e0) * 627.509 for r in pes], "x", color=C["GVA"],
            ms=4, label="GVA 2L")
    ax.set_xlabel("Pt–Cl distance (Å)")
    ax.set_ylabel("Energy relative to CASCI minimum (kcal mol$^{-1}$)")
    ax.legend(loc="upper left")
    style.panel(ax, "a", x=-0.28)
    ax = axs[1]
    ef = np.array([r["E_CASCI"] for r in pes])
    for key, col, lab in (("E_CISD", C["CISD"], "CISD"), ("E_CCSD", C["CCSD"], "CCSD"),
                          ("E_CCSD(T)", C["CCSDT"], "CCSD(T)")):
        ax.plot(R, mha([r[key] for r in pes] - ef), "o-", ms=2, color=col, label=lab)
    for key, col, ls, lab in (("UCCSD_1L", C["UCCSD"], "--", "UCCSD 1L"),
                              ("GVA_1L", C["GVA"], "--", "GVA 1L"),
                              ("GVA_2L", C["GVA"], "-", "GVA 2L")):
        ax.plot(R, mha([r[key]["dE"] for r in pes]), ls, marker="s", ms=2, color=col,
                label=lab, lw=0.9)
    ax.axhline(CHEM, color=C["grey"], ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(0.015, 40)
    ax.set_xlabel("Pt–Cl distance (Å)")
    ax.set_ylabel("|ΔE| relative to CASCI (mHa)")
    ax.legend(loc="upper center", fontsize=5.5, ncol=2, frameon=True, framealpha=1,
              edgecolor="none")
    style.panel(ax, "b", x=-0.28)
    ax = axs[2]
    ax.plot(R, [r["c0_sq"] for r in pes], "o-", ms=2, color=C["black"], label="$c_0^2$ (CASCI)")
    ax.set_xlabel("Pt–Cl distance (Å)")
    ax.set_ylabel("Reference weight $c_0^2$")
    sec = ax.twinx()
    sec.spines["right"].set_visible(True)
    sec.plot(R, [(r["E_ref"] - r["E_CASCI"]) * 1e3 for r in pes], "s--", ms=2,
             color=C["GVA"], label="−$E_{corr}$ (active space)")
    sec.set_ylabel("−$E_{corr}$, active space (mHa)", color=C["GVA"])
    sec.tick_params(axis="y", colors=C["GVA"])
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = sec.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=5.5)
    style.panel(ax, "c", x=-0.28)
    fig.tight_layout(w_pad=1.8)
    save(fig, "Fig8")


def figS1():
    fig, ax = plt.subplots(figsize=(style.SINGLE, 62 * style.MM))
    for r in sm["runs"]:
        if r["ansatz"] not in ("UCCSD", "GVA"):
            continue
        col = C[r["ansatz"]]
        ax.plot(max(abs(r["dE"]) * 1e3, FLOOR), max(r["S2"], 1e-12),
                "o" if r["layers"] == 1 else "s", ms=3, mec=col,
                mfc="white" if r["ansatz"] == "UCCSD" else col, mew=0.7, alpha=0.85)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(CHEM, color=C["grey"], ls=":", lw=0.8)
    ax.set_xlabel("|ΔE| relative to FCI (mHa)")
    ax.set_ylabel("⟨S$^2$⟩ of final state")
    handles = [Line2D([], [], marker="o", ls="", mec=C["UCCSD"], mfc="white", label="UCCSD 1L"),
               Line2D([], [], marker="s", ls="", mec=C["UCCSD"], mfc="white", label="UCCSD 2L"),
               Line2D([], [], marker="o", ls="", color=C["GVA"], label="GVA 1L"),
               Line2D([], [], marker="s", ls="", color=C["GVA"], label="GVA 2L")]
    ax.legend(handles=handles, loc="upper left")
    save(fig, "FigS1")


if __name__ == "__main__":
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    rows = noise_table()
    json.dump(rows, open("results/noise_budget.json", "w"), indent=1, default=float)
    fig7(rows)
    fig8()
    figS1()
