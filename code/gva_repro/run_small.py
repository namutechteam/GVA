"""Small-molecule benchmark: H2, H4, LiH, BeH2 (STO-3G), three geometries each.

Ansaetze: UCCSD and GVA, 1 and 2 layers, 5 initialisations each
(seed 0 = Hartree-Fock start, theta = 0; seeds 1-4 = N(0, 0.05^2) start).
Ablation (1 layer, HF start): PG screening only, spin tying only.
ADAPT-VQE with the spin-complemented, symmetry-screened pool.
"""
import json
import pickle
import time

import numpy as np

from datetime import datetime, timezone

from gva.systems import small_molecule_hamiltonian, classical_references, set_pyscf_log
from gva.ansatz import (uccsd_pool, gva_groups, counts_summary, cnot_count_ladder,
                        parameter_map, excitation_label)


def now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")

from gva.engine import DeterminantSpace, Ansatz, run_vqe, adapt_vqe

MOLS = ["H2", "H4", "LiH", "BeH2"]
SCALES = [1.0, 1.5, 2.0]
SEEDS = range(5)
VARIANTS = {"UCCSD": (False, False), "GVA": (True, True),
            "PG only": (True, False), "spin only": (False, True)}

records, systems = [], []
t_all = time.time()
for name in MOLS:
    for sc in SCALES:
        set_pyscf_log("results/logs/pyscf", f"{name}_x{sc:.1f}")
        H = small_molecule_hamiltonian(name, sc)
        ref = classical_references(H)
        space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
        pool = uccsd_pool(H.nocc, H.norb)
        cs = counts_summary(H.nocc, H.norb, H.orbsym)
        ef = ref["E_FCI"]
        systems.append({
            "system": name, "scale": sc, "n_qubits": 2 * H.norb,
            "n_electrons": sum(H.nelec), "irreps": H.irrep_names,
            "E_RHF": ref["E_ref"], "E_CISD": ref["E_CISD"],
            "E_CCSD": ref["E_CCSD"], "E_CCSD(T)": ref["E_CCSD(T)"],
            "E_FCI": ef, "E_corr": ef - ref["E_ref"],
            "c0_sq": ref["c0_ref_weight"], "counts": cs,
            "parameter_maps_1L": {v: parameter_map(gva_groups(pool, H.orbsym, *f))
                                  for v, f in VARIANTS.items()}})
        for var, (pg, spn) in VARIANTS.items():
            groups = gva_groups(pool, H.orbsym, pg, spn)
            layer_list = [1, 2] if var in ("UCCSD", "GVA") else [1]
            for L in layer_list:
                seeds = SEEDS if var in ("UCCSD", "GVA") else [0]
                for s in seeds:
                    A = Ansatz(space, groups, layers=L)
                    th0 = np.zeros(A.n_params) if s == 0 else \
                        np.random.default_rng(1000 + s).normal(0, 0.05, A.n_params)
                    t_start = now()
                    r = run_vqe(A, th0, e_exact=ef)
                    records.append({
                        "system": name, "scale": sc, "ansatz": var, "layers": L,
                        "seed": s, "n_params": A.n_params,
                        "cnot": L * cs[var]["cnot_1L"], "E": r["E"],
                        "dE": r["dE"], "nit": r["nit"], "nfev": r["nfev"],
                        "iter_to_chem_acc": r["iter_to_chem_acc"],
                        "S2": r["S2"], "wall_s": r["wall_s"],
                        "grad_norm": r["grad_norm"], "message": r["message"],
                        "started": t_start, "finished": now(),
                        "theta": r["theta"].tolist(),
                        "history": r["history"].tolist()})
        t_start = now()
        ad = adapt_vqe(space, gva_groups(pool, H.orbsym, True, True), e_exact=ef,
                       grad_tol=1e-3, max_ops=120)
        records.append({"system": name, "scale": sc, "ansatz": "ADAPT-VQE",
                        "layers": None, "seed": 0, "n_params": ad["n_ops"],
                        "cnot": sum(cnot_count_ladder(e) for g in ad["chosen"] for e in g),
                        "E": ad["E"], "dE": ad["dE"], "trace": ad["trace"], "stop_reason": ad["stop_reason"],
                        "wall_s": ad["wall_s"], "started": t_start, "finished": now(),
                        "operators": [[excitation_label(e) for e in g] for g in ad["chosen"]],
                        "theta": np.asarray(ad["theta"]).tolist()})
        print(f"{name} {sc}: done ({time.time()-t_all:.0f}s)", flush=True)

json.dump({"systems": systems, "runs": records},
          open("results/small_molecules.json", "w"), indent=1, default=float)
print("total", time.time() - t_all)
