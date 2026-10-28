"""Every number quoted in the manuscript and SI is computed here from results/."""
import json
from collections import defaultdict

import numpy as np

HA2KCAL = 627.509


def load():
    sm = json.load(open("results/small_molecules.json"))
    cis = json.load(open("results/cisplatin_eq.json"))
    meas = json.load(open("results/measurement.json"))
    scal = json.load(open("results/scaling.json"))
    pes = json.load(open("results/pes.json"))
    noise = json.load(open("results/noise_budget.json"))
    geo = json.load(open("results/geometries.json"))
    return sm, cis, meas, scal, pes, noise, geo


def orbital_sensitivity():
    return json.load(open("results/orbital_sensitivity.json"))


def fmt_sig(x, n=2):
    """Round to n significant figures, plain notation."""
    if x == 0:
        return "0"
    from math import floor, log10
    d = n - 1 - int(floor(log10(abs(x))))
    v = round(x, d)
    return f"{v:.{max(d, 0)}f}"


def compute():
    sm, cis, meas, scal, pes, noise, geo = load()
    F = {}
    # ------------------------------------------------ small molecules
    g = defaultdict(list)
    for r in sm["runs"]:
        g[(r["system"], r["scale"], r["ansatz"], r["layers"])].append(r)
    F["small_groups"] = g
    F["small_systems"] = sm["systems"]
    sysd = {(s["system"], s["scale"]): s for s in sm["systems"]}
    F["sysd"] = sysd
    # parameter / CNOT reductions at scale 1.0
    red = {}
    for m in ("H2", "H4", "LiH", "BeH2"):
        c = sysd[(m, 1.0)]["counts"]
        red[m] = (c["UCCSD"]["n_params"], c["GVA"]["n_params"],
                  c["UCCSD"]["cnot_1L"], c["GVA"]["cnot_1L"])
    F["red"] = red
    # max |GVA1L - UCCSD1L| for HF start (seed 0)
    diffs = []
    for m in ("H2", "H4", "LiH", "BeH2"):
        for sc in (1.0, 1.5, 2.0):
            u = [r for r in g[(m, sc, "UCCSD", 1)] if r["seed"] == 0][0]["dE"]
            v = [r for r in g[(m, sc, "GVA", 1)] if r["seed"] == 0][0]["dE"]
            diffs.append((v - u, m, sc))
    F["gva_minus_uccsd_1L_seed0"] = diffs
    F["max_gva_penalty_1L_mHa"] = max(d for d, _, _ in diffs) * 1e3
    # spin contamination events
    cont = []
    for (m, sc, an, L), rr in g.items():
        if an in ("UCCSD", "GVA"):
            for r in rr:
                if r["S2"] > 1e-2:
                    cont.append((m, sc, an, L, r["seed"], r["S2"], r["dE"]))
    F["spin_contaminated"] = cont
    F["gva_max_S2"] = max(r["S2"] for (m, sc, an, L), rr in g.items() if an == "GVA" for r in rr)
    # 2L medians
    two = {}
    for m in ("H2", "H4", "LiH", "BeH2"):
        for sc in (1.0, 1.5, 2.0):
            for an in ("UCCSD", "GVA"):
                e = np.array([r["dE"] for r in g[(m, sc, an, 2)]]) * 1e3
                it = np.array([r["nit"] for r in g[(m, sc, an, 2)]])
                two[(m, sc, an)] = (np.median(e), e.min(), e.max(), np.median(it))
    F["two"] = two
    F["max_2L_err_mHa"] = max(v[2] for v in two.values())
    one_fail = sorted({(m, sc) for (m, sc, an, L), rr in g.items()
                       if an == "GVA" and L == 1 and min(abs(r["dE"]) for r in rr) * 1e3 > 1.6})
    F["gva1L_fail"] = one_fail
    # lowest triplet gap of stretched H4 (FCI, Sz = 0 sector)
    from pyscf import fci
    from pyscf.fci import spin_op
    from gva.systems import small_molecule_hamiltonian
    Hh = small_molecule_hamiltonian("H4", 2.0)
    ee, cc = fci.direct_spin1.kernel(Hh.h1, Hh.eri, Hh.norb, Hh.nelec, nroots=4)
    s2 = [spin_op.spin_square(c, Hh.norb, Hh.nelec)[0] for c in cc]
    trip = min(e for e, s in zip(ee, s2) if abs(s - 2) < 1e-6)
    F["h4_triplet_gap_mHa"] = (trip - ee[0]) * 1e3
    # ------------------------------------------------ cisplatin
    F["cis"] = cis
    ef = cis["E_CASCI"]
    F["cis_corr_mHa"] = (ef - cis["E_ref"]) * 1e3
    for k in ("E_CISD", "E_CCSD", "E_CCSD(T)"):
        F[f"cis_{k}_err_mHa"] = (cis[k] - ef) * 1e3
    cr = defaultdict(list)
    for r in cis["runs"]:
        cr[(r["ansatz"], r["layers"])].append(r)
    F["cis_runs"] = cr
    F["cis_med"] = {k: (np.median([r["dE"] for r in v]) * 1e3,
                        min(r["dE"] for r in v) * 1e3, max(r["dE"] for r in v) * 1e3,
                        v[0]["n_params"], v[0]["cnot"], len(v),
                        np.median([r["nit"] for r in v]) if "nit" in v[0] else None,
                        np.median([r["wall_s"] for r in v]))
                    for k, v in cr.items()}
    # ------------------------------------------------ measurement
    F["meas"] = meas
    # ------------------------------------------------ scaling
    F["scal"] = scal
    n = np.array([r["n_active_orbitals"] for r in scal])
    for key in ("UCCSD_params", "GVA_params", "UCCSD_cnot_1L", "GVA_cnot_1L"):
        y = np.array([r[key] for r in scal], float)
        F[f"slope_{key}"] = np.polyfit(np.log(n[3:]), np.log(y[3:]), 1)[0]
    F["ratio_params"] = [r["GVA_params"] / r["UCCSD_params"] for r in scal]
    F["ratio_cnot"] = [r["GVA_cnot_1L"] / r["UCCSD_cnot_1L"] for r in scal]
    # ------------------------------------------------ PES
    F["pes"] = pes
    # ------------------------------------------------ noise
    F["noise"] = noise
    F["geo"] = geo
    return F


if __name__ == "__main__":
    F = compute()
    for k, v in F.items():
        if k in ("small_groups", "cis", "cis_runs", "meas", "scal", "pes", "noise", "geo",
                 "small_systems", "sysd"):
            continue
        print(k, v)


def S_names(S):
    return [s["system"] for s in S]


def periodic():
    """Numbers for the periodic section (results/periodic.json and friends)."""
    import os
    per = json.load(open("results/periodic.json"))
    P = {"per": per, "sys": {s["system"]: s for s in per["systems"]}}
    P["scal"] = json.load(open("results/scaling_periodic.json"))
    for name in ("measurement_periodic", "tying_check_periodic", "verify_periodic"):
        f = f"results/{name}.json"
        P[name] = json.load(open(f)) if os.path.exists(f) else None

    def tag(r, d, nk):
        return f"H2chain_r{r:.2f}_d{d:.2f}_Nk{nk}"
    P["tag"] = tag
    g = defaultdict(list)
    for r in per["runs"]:
        g[(r["system"], r["ansatz"], r["layers"])].append(r)
    P["g"] = g

    def stat(system, ansatz, layers):
        rr = g[(system, ansatz, layers)]
        e = np.array([x["dE"] for x in rr]) * 1e3
        hf = [x for x in rr if x["seed"] == 0][0]
        return {"hf": hf["dE"] * 1e3, "min": e.min(), "med": float(np.median(e)), "max": e.max(),
                "n": len(rr), "S2max": max(x["S2"] for x in rr), "params": hf["n_params"],
                "cnot": hf["cnot"], "ov_hf": hf["overlap_sq_FCI"], "runs": rr}
    P["stat"] = stat
    # validation maxima
    S = per["systems"]
    P["max_fci_bloch_diff"] = max(abs(s["E_FCI_bloch"] - s["E_FCI"]) for s in S)
    P["max_fci_supercell_diff"] = max(abs(s["E_FCI_supercell_gamma"] - s["E_FCI"]) for s in S)
    P["max_hf_diff"] = max(abs(s["E_HF_from_integrals"] - s["E_KRHF_supercell"]) for s in S)
    P["max_mom_viol"] = max(s["momentum_violation_max"] for s in S)
    P["max_tr_viol"] = max(s["time_reversal_max"] for s in S)
    P["max_outside_K"] = max(1 - s["FCI_momentum_sector_weights"][s["K_HF"]] for s in S)
    P["max_forbidden_grad"] = max(r.get("max_grad_forbidden", 0.0) for r in per["runs"])
    P["n_runs"] = len(per["runs"])
    # 1L tying penalty (GVA-k minus k only, HF start) and GVA-k minus UCCSD (real) at 1L
    pen, d1 = [], []
    for s in S:
        t = s["system"]
        e_gk = stat(t, "GVA-k", 1)["hf"]
        pen.append(e_gk - stat(t, "k only", 1)["hf"])
        d1.append(e_gk - stat(t, "UCCSD", 1)["hf"])
    P["tying_penalty_1L_mHa"] = (min(pen), max(pen))
    P["gva_minus_uccsd_1L_mHa"] = (min(d1), max(d1))
    # correlation regimes: weak = reference weight above 0.85 (dimerized, 1.0 A)
    weak = [s["system"] for s in S if s["c0_sq"] > 0.85]
    strong = [s["system"] for s in S if s["c0_sq"] <= 0.85]
    P["weak"], P["strong"] = weak, strong
    P["min_c0_weak"] = min(P["sys"][t]["c0_sq"] for t in weak)
    dw = [stat(t, "GVA-k", 1)["hf"] - stat(t, "UCCSD", 1)["hf"] for t in weak]
    P["d1_weak_uHa"] = max(abs(x) for x in dw) * 1e3
    P["d1_strong_mHa"] = {t: stat(t, "GVA-k", 1)["hf"] - stat(t, "UCCSD", 1)["hf"] for t in strong}
    P["tr_cost_mHa"] = {t: stat(t, "k + TR", 1)["hf"] - stat(t, "k only", 1)["hf"] for t in S_names(S)}
    P["spin_cost_mHa"] = {t: stat(t, "k + spin", 1)["hf"] - stat(t, "k only", 1)["hf"] for t in S_names(S)}
    gross = defaultdict(int)
    s2max = defaultdict(float)
    for r in per["runs"]:
        gross[r["ansatz"]] += r["S2"] > 0.5
        s2max[r["ansatz"]] = max(s2max[r["ansatz"]], r["S2"])
    P["ratio_2L_strong"] = {t: stat(t, "GVA-k", 2)["min"] / stat(t, "UCCSD", 2)["min"] for t in strong}
    P["S2_gross"] = dict(gross)
    P["S2_max"] = dict(s2max)
    # spin contamination counts (<S^2> > 1e-2) by ansatz
    cont = defaultdict(int)
    tot = defaultdict(int)
    for r in per["runs"]:
        key = r["ansatz"]
        tot[key] += 1
        cont[key] += r["S2"] > 1e-2
    P["contaminated"] = dict(cont)
    P["n_by_ansatz"] = dict(tot)
    return P
