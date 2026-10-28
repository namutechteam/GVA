"""Cisplatin CAS(10e,10o) benchmark at the C2v reference geometry.

Writes results/cisplatin_eq.json (incrementally) and the active-space
Hamiltonian to results/cisplatin_cas10_10.npz.
"""
import json
import time

import numpy as np

from datetime import datetime, timezone

from gva.systems import (cisplatin_cas_hamiltonian, classical_references, cisplatin_atoms,
                         set_pyscf_log)
from gva.ansatz import (uccsd_pool, gva_groups, counts_summary, cnot_count_ladder,
                        parameter_map, excitation_label)


def now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")

from gva.engine import DeterminantSpace, Ansatz, run_vqe, adapt_vqe

OUT = "results/cisplatin_eq.json"
t_all = time.time()
set_pyscf_log("results/logs/pyscf", "cisplatin_C2v_eq")
H = cisplatin_cas_hamiltonian()
ref = classical_references(H)
ef = ref["E_FCI"]
np.savez("results/cisplatin_cas10_10.npz", h1=H.h1, eri=H.eri, ecore=H.ecore,
         orbsym=H.orbsym, nelec=np.array(H.nelec), civec=ref["civec"],
         C_act_ao=H.extra["C_act_ao"], no_occ=H.extra["act_no_occ"])
space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
pool = uccsd_pool(H.nocc, H.norb)
cs = counts_summary(H.nocc, H.norb, H.orbsym)
mol = H.extra["mol"]
summary = {
    "geometry_xyz_angstrom": [(a, list(map(float, c))) for a, c in cisplatin_atoms()],
    "point_group": mol.groupname, "n_electrons_total": int(mol.nelectron),
    "n_ao": int(mol.nao), "E_RHF": H.e_rhf_total,
    "E_MP2_corr": H.extra["e_mp2_corr"], "n_core_orbitals": H.extra["n_core"],
    "active_irreps": H.irrep_names, "active_NO_occupations": H.extra["act_no_occ"].tolist(),
    "active_orbital_energies": H.extra["active_eps"].tolist(),
    "E_ref": ref["E_ref"], "E_CISD": ref["E_CISD"], "E_CCSD": ref["E_CCSD"],
    "E_CCSD(T)": ref["E_CCSD(T)"], "E_CASCI": ef, "c0_sq": ref["c0_ref_weight"],
    "T1_diag": ref["t1_diag"], "counts": cs, "dim": space.dim,
    "parameter_maps_1L": {v: parameter_map(gva_groups(pool, H.orbsym, *f)) for v, f in
                          {"UCCSD": (False, False), "GVA": (True, True),
                           "PG only": (True, False), "spin only": (False, True)}.items()},
    "runs": []}


def dump():
    json.dump(summary, open(OUT, "w"), indent=1, default=float)


import os
done = set()
if os.path.exists(OUT):
    old = json.load(open(OUT))
    summary["runs"] = [r for r in old.get("runs", []) if r["ansatz"] != "ADAPT-VQE"]
    done = {(r["ansatz"], r["layers"], r["seed"]) for r in summary["runs"]}
dump()
plan = []
SEEDS = {("GVA", 1): (0, 1, 2), ("GVA", 2): (0, 1, 2), ("GVA", 3): (0,),
         ("UCCSD", 1): (0, 1, 2), ("UCCSD", 2): (0, 1), ("UCCSD", 3): (0,)}
for var in ("GVA", "UCCSD"):
    for L in (1, 2, 3):
        plan += [(var, L, s) for s in SEEDS[(var, L)]]
plan += [("PG only", 1, 0), ("spin only", 1, 0)]
VAR = {"UCCSD": (False, False), "GVA": (True, True), "PG only": (True, False),
       "spin only": (False, True)}
for var, L, s in plan:
    if (var, L, s) in done:
        continue
    groups = gva_groups(pool, H.orbsym, *VAR[var])
    A = Ansatz(space, groups, layers=L)
    th0 = np.zeros(A.n_params) if s == 0 else \
        np.random.default_rng(2000 + s).normal(0, 0.05, A.n_params)
    t_start = now()
    r = run_vqe(A, th0, e_exact=ef, maxiter=2500)
    summary["runs"].append({
        "ansatz": var, "layers": L, "seed": s, "n_params": A.n_params,
        "cnot": L * cs[var]["cnot_1L"], "E": r["E"], "dE": r["dE"],
        "nit": r["nit"], "nfev": r["nfev"], "iter_to_chem_acc": r["iter_to_chem_acc"],
        "S2": r["S2"], "wall_s": r["wall_s"], "grad_norm": r["grad_norm"],
        "message": r["message"], "started": t_start, "finished": now(),
        "theta": r["theta"].tolist(), "history": r["history"].tolist()})
    dump()
    print(var, L, s, f"dE={r['dE']*1e3:.4f} mHa nit={r['nit']} t={r['wall_s']:.0f}s",
          flush=True)

t_start = now()
ad = adapt_vqe(space, gva_groups(pool, H.orbsym, True, True), e_exact=ef,
               grad_tol=1e-3, max_ops=120, verbose=True)
summary["runs"].append({"ansatz": "ADAPT-VQE", "layers": None, "seed": 0,
                        "n_params": ad["n_ops"],
                        "cnot": sum(cnot_count_ladder(e) for g in ad["chosen"] for e in g),
                        "E": ad["E"], "dE": ad["dE"], "trace": ad["trace"], "stop_reason": ad["stop_reason"],
                        "wall_s": ad["wall_s"], "started": t_start, "finished": now(),
                        "operators": [[excitation_label(e) for e in g] for g in ad["chosen"]],
                        "theta": np.asarray(ad["theta"]).tolist()})
dump()
print("total", time.time() - t_all)
