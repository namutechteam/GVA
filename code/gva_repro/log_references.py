"""Re-run every electronic-structure reference step with full PySCF output.

The production scripts ran PySCF with verbose=0, so no SCF/MP2/CC output was
kept. This script repeats exactly the same code path (the gva.systems functions)
with PySCF's verbose level raised to 4, writes one log per system to
results/logs/reference/, and checks every energy against the values stored in
results/*.json. Any deviation larger than 1e-7 Ha raises an error.

VQE optimizations are not repeated here (see PROVENANCE.md).
"""
import json
import os
import time

import numpy as np
from pyscf import gto

import gva.systems as S

LOGDIR = "results/logs/reference"
os.makedirs(LOGDIR, exist_ok=True)
TOL = 1e-7
_orig_M = gto.M
_state = {"stem": None, "k": 0}


def _logged_M(*args, **kw):
    _state["k"] += 1
    kw["verbose"] = 4
    kw["output"] = os.path.join(LOGDIR, f"{_state['stem']}__part{_state['k']}.log")
    return _orig_M(*args, **kw)


gto.M = _logged_M          # S.gto is pyscf.gto, so the production code is unchanged


def start(stem):
    _state["stem"], _state["k"] = stem, 0


def merge(stem, summary_lines):
    parts = sorted(f for f in os.listdir(LOGDIR) if f.startswith(stem + "__part"))
    with open(os.path.join(LOGDIR, stem + ".log"), "w") as out:
        out.write(f"# {stem}: reference calculations re-run {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n")
        for line in summary_lines:
            out.write("# " + line + "\n")
        for p in parts:
            out.write(f"\n######## {p} ########\n")
            out.write(open(os.path.join(LOGDIR, p)).read())
            os.remove(os.path.join(LOGDIR, p))


def check(label, new, old):
    d = abs(new - old)
    status = "OK" if d < TOL else "MISMATCH"
    line = f"{label:14s} new={new:.10f}  stored={old:.10f}  |diff|={d:.1e}  {status}"
    if d >= TOL:
        raise AssertionError(line)
    return line


report = []
# ------------------------------------------------------------ small molecules
sm = json.load(open("results/small_molecules.json"))
for s in sm["systems"]:
    stem = f"{s['system']}_x{s['scale']:.1f}"
    start(stem)
    H = S.small_molecule_hamiltonian(s["system"], s["scale"])
    ref = S.classical_references(H)
    lines = [check("E_RHF", ref["E_ref"], s["E_RHF"]),
             check("E_CISD", ref["E_CISD"], s["E_CISD"]),
             check("E_CCSD", ref["E_CCSD"], s["E_CCSD"]),
             check("E_CCSD(T)", ref["E_CCSD(T)"], s["E_CCSD(T)"]),
             check("E_FCI", ref["E_FCI"], s["E_FCI"])]
    merge(stem, lines)
    report += [stem] + ["  " + l for l in lines]
    print(stem, "ok", flush=True)

# ------------------------------------------------------------ cisplatin, C2v
cis = json.load(open("results/cisplatin_eq.json"))
start("cisplatin_C2v_eq")
H = S.cisplatin_cas_hamiltonian()
ref = S.classical_references(H)
lines = [check("E_RHF", H.e_rhf_total, cis["E_RHF"]),
         check("E_MP2_corr", H.extra["e_mp2_corr"], cis["E_MP2_corr"]),
         check("E_ref", ref["E_ref"], cis["E_ref"]),
         check("E_CISD", ref["E_CISD"], cis["E_CISD"]),
         check("E_CCSD", ref["E_CCSD"], cis["E_CCSD"]),
         check("E_CCSD(T)", ref["E_CCSD(T)"], cis["E_CCSD(T)"]),
         check("E_CASCI", ref["E_FCI"], cis["E_CASCI"]),
         f"point group {H.extra['mol'].groupname}, n_core {H.extra['n_core']}, "
         f"active irreps {H.irrep_names}"]
d = np.load("results/cisplatin_cas10_10.npz")
lines.append(f"max |h1 - stored h1| = {np.abs(H.h1 - d['h1']).max():.1e}, "
             f"max |eri - stored eri| = {np.abs(H.eri - d['eri']).max():.1e}")
assert np.abs(H.h1 - d["h1"]).max() < 1e-8 and np.abs(H.eri - d["eri"]).max() < 1e-8
merge("cisplatin_C2v_eq", lines)
report += ["cisplatin_C2v_eq"] + ["  " + l for l in lines]
print("cisplatin eq ok", flush=True)

# ------------------------------------------------------------ Pt-Cl scan, Cs
pes = {round(r["R_PtCl"], 3): r for r in json.load(open("results/pes.json"))}
R_EQ = 2.32
prev = None
order = [2.32, 2.5, 2.7, 2.9, 3.2, 3.5, 3.8, "eq", 2.2, 2.0]
states = {}
for r in order:
    if r == "eq":
        prev = states[R_EQ]
        continue
    stem = f"cisplatin_Cs_R{r:.2f}"
    start(stem)
    H = S.cisplatin_cas_hamiltonian(r1=r, r2=R_EQ, prev=prev, label=stem, symmetry="Cs")
    ref = S.classical_references(H)
    old = pes[round(r, 3)]
    lines = [check("E_RHF", H.e_rhf_total, old["E_RHF"]),
             check("E_ref", ref["E_ref"], old["E_ref"]),
             check("E_CISD", ref["E_CISD"], old["E_CISD"]),
             check("E_CCSD", ref["E_CCSD"], old["E_CCSD"]),
             check("E_CCSD(T)", ref["E_CCSD(T)"], old["E_CCSD(T)"]),
             check("E_CASCI", ref["E_FCI"], old["E_CASCI"]),
             f"tracked from R = {('none (occupation selection)' if prev is None else 'previous point')}"]
    merge(stem, lines)
    report += [stem] + ["  " + l for l in lines]
    states[r] = {"mol": H.extra["mol"], "mf": H.extra["mf"], "C_act_ao": H.extra["C_act_ao"]}
    prev = states[r]
    print(stem, "ok", flush=True)

with open(os.path.join(LOGDIR, "SUMMARY.txt"), "w") as f:
    f.write("Reference energies recomputed with full PySCF output and compared with "
            f"stored values (tolerance {TOL:g} Ha)\n\n")
    f.write("\n".join(report) + "\n")
print("all reference energies reproduced")
