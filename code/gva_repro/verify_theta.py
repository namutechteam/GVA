"""Rebuild every final VQE / ADAPT-VQE state from the stored parameters.

For each run the ansatz is reconstructed from its definition (variant, layers,
or, for ADAPT-VQE, the stored operator labels), the state U(theta)|Phi_0> is
prepared with the stored theta, and its energy and <S^2> are compared with the
values written by the run scripts. Tolerances: 1e-9 Ha for the energy, 1e-8
for <S^2>. Writes results/logs/verify_theta.log and
results/verify_theta.json; exits non-zero on any mismatch.
"""
import json
import sys

import numpy as np

from gva.systems import (small_molecule_hamiltonian, cisplatin_cas_hamiltonian,
                         ActiveHamiltonian)
from gva.ansatz import uccsd_pool, gva_groups, excitation_label
from gva.engine import DeterminantSpace, Ansatz

TOL_E, TOL_S2 = 1e-9, 1e-8
VAR = {"UCCSD": (False, False), "GVA": (True, True), "PG only": (True, False),
       "spin only": (False, True)}
lines, rows, bad = [], [], 0


def adapt_groups(pool, labels):
    lut = {excitation_label(e): e for e in pool}
    return [tuple(lut[l] for l in g) for g in labels]


def check(tag, space, groups, layers, theta, E_ref, S2_ref):
    global bad
    A = Ansatz(space, groups, layers=layers)
    th = np.asarray(theta)
    assert len(th) == A.n_params, (tag, len(th), A.n_params)
    psi = A.state(th)
    E = space.energy(psi)
    s2 = space.s2(psi)
    dE = abs(E - E_ref)
    dS = abs(s2 - S2_ref) if S2_ref is not None else 0.0
    ok = dE < TOL_E and dS < TOL_S2
    bad += 0 if ok else 1
    lines.append(f"{tag:45s} n_par={A.n_params:5d}  |E - E_stored|={dE:.1e} Ha  "
                 f"|S2 - S2_stored|={dS:.1e}  {'OK' if ok else 'MISMATCH'}")
    rows.append({"run": tag, "n_params": A.n_params, "E_rebuilt": E, "E_stored": E_ref,
                 "abs_dE": dE, "S2_rebuilt": s2, "S2_stored": S2_ref, "ok": ok})


# ------------------------------------------------------------ small molecules
sm = json.load(open("results/small_molecules.json"))
for s in sm["systems"]:
    H = small_molecule_hamiltonian(s["system"], s["scale"])
    space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    pool = uccsd_pool(H.nocc, H.norb)
    for r in sm["runs"]:
        if r["system"] != s["system"] or r["scale"] != s["scale"]:
            continue
        tag = f"{r['system']} x{r['scale']:.1f} {r['ansatz']} {r['layers']}L s{r['seed']}"
        if r["ansatz"] == "ADAPT-VQE":
            groups = adapt_groups(pool, r["operators"])
            check(tag, space, groups, 1, r["theta"], r["E"], None)
        else:
            check(tag, space, gva_groups(pool, H.orbsym, *VAR[r["ansatz"]]), r["layers"],
                  r["theta"], r["E"], r["S2"])
    print(s["system"], s["scale"], "checked", flush=True)

# ------------------------------------------------------------ cisplatin C2v
d = np.load("results/cisplatin_cas10_10.npz")
H = ActiveHamiltonian(label="cis", norb=10, nelec=tuple(int(x) for x in d["nelec"]),
                      ecore=float(d["ecore"]), h1=d["h1"], eri=d["eri"],
                      orbsym=d["orbsym"], irrep_names=[], e_rhf_total=np.nan)
space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
pool = uccsd_pool(H.nocc, H.norb)
cis = json.load(open("results/cisplatin_eq.json"))
for r in cis["runs"]:
    tag = f"cisplatin C2v {r['ansatz']} {r['layers']}L s{r['seed']}"
    if r["ansatz"] == "ADAPT-VQE":
        check(tag, space, adapt_groups(pool, r["operators"]), 1, r["theta"], r["E"], None)
    else:
        check(tag, space, gva_groups(pool, H.orbsym, *VAR[r["ansatz"]]), r["layers"],
              r["theta"], r["E"], r["S2"])
print("cisplatin checked", flush=True)

# ------------------------------------------------------------ Pt-Cl scan
pes = {round(r["R_PtCl"], 3): r for r in json.load(open("results/pes.json"))}
prev, states = None, {}
for r in [2.32, 2.5, 2.7, 2.9, 3.2, 3.5, 3.8, "eq", 2.2, 2.0]:
    if r == "eq":
        prev = states[2.32]
        continue
    H = cisplatin_cas_hamiltonian(r1=r, r2=2.32, prev=prev, symmetry="Cs")
    states[r] = prev = {"mol": H.extra["mol"], "mf": H.extra["mf"],
                        "C_act_ao": H.extra["C_act_ao"]}
    space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    pool = uccsd_pool(H.nocc, H.norb)
    rec = pes[round(r, 3)]
    for key, var, L in (("GVA_1L", "GVA", 1), ("GVA_2L", "GVA", 2), ("UCCSD_1L", "UCCSD", 1)):
        check(f"Pt-Cl R={r:.2f} {key}", space, gva_groups(pool, H.orbsym, *VAR[var]), L,
              rec[key]["theta"], rec[key]["E"], rec[key]["S2"])
    print("R", r, "checked", flush=True)

with open("results/logs/verify_theta.log", "w") as f:
    f.write(f"Rebuilt {len(rows)} final states from stored parameters; "
            f"{len(rows) - bad} OK, {bad} mismatches "
            f"(tolerances {TOL_E:g} Ha, {TOL_S2:g} for <S^2>)\n\n")
    f.write("\n".join(lines) + "\n")
json.dump(rows, open("results/verify_theta.json", "w"), indent=1, default=float)
print(f"{len(rows)} states rebuilt, {bad} mismatches")
sys.exit(1 if bad else 0)
