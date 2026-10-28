"""Pt-Cl bond-stretching scan of cisplatin, CAS(10e,10o).

One Pt-Cl distance is varied (the other ligands fixed), lowering the point
group to Cs. The whole scan, including 2.32 A, is run in Cs so that the ansatz
is defined identically at every geometry. Active orbitals are chosen by
occupation at 2.32 A and then tracked by maximum overlap with the previous
geometry, outward and inward. Phase 1 builds all active-space Hamiltonians
(with tracking); phase 2 runs the classical references and VQE, resuming from
results/pes.json if present.
"""
import json
import os
import time

import numpy as np
from pyscf import gto

from datetime import datetime, timezone

from gva.systems import cisplatin_cas_hamiltonian, classical_references, set_pyscf_log
from gva.ansatz import uccsd_pool, gva_groups, counts_summary


def now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")

from gva.engine import DeterminantSpace, Ansatz, run_vqe

R_EQ = 2.32
OUTWARD = [2.32, 2.5, 2.7, 2.9, 3.2, 3.5, 3.8]
INWARD = [2.2, 2.0]
OUT = "results/pes.json"

rows = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {round(r["R_PtCl"], 3) for r in rows}


def tracking_overlap(prev, H):
    """Smallest singular value of the overlap between successive active spaces."""
    S = gto.intor_cross("int1e_ovlp", prev["mol"], H.extra["mol"])
    O = prev["C_act_ao"].T @ S @ H.extra["C_act_ao"]
    return float(np.linalg.svd(O, compute_uv=False).min())


# ---------------------------------------------------------------- phase 1
hams = {}
prev = None
for path in (OUTWARD, INWARD):
    prev = hams[R_EQ]["state"] if path is INWARD else None
    for r in path:
        if r in hams:
            continue
        set_pyscf_log("results/logs/pyscf", f"cisplatin_Cs_R{r:.2f}_build")
        H = cisplatin_cas_hamiltonian(r1=r, r2=R_EQ, prev=prev,
                                      label=f"cisplatin_R{r:.2f}", symmetry="Cs")
        ov = tracking_overlap(prev, H) if prev is not None else 1.0
        state = {"mol": H.extra["mol"], "mf": H.extra["mf"], "C_act_ao": H.extra["C_act_ao"]}
        hams[r] = {"H": H, "state": state, "overlap": ov}
        prev = state
        print(f"built R={r:.2f} overlap={ov:.4f}", flush=True)

# ---------------------------------------------------------------- phase 2
for r in sorted(hams):
    H = hams[r]["H"]
    set_pyscf_log("results/logs/pyscf", f"cisplatin_Cs_R{r:.2f}_refs")
    ref = classical_references(H)
    ef = ref["E_FCI"]
    if round(r, 3) in done:
        for row in rows:          # refresh tracking metadata of resumed points
            if round(row["R_PtCl"], 3) == round(r, 3):
                row["tracking_min_singular_value"] = hams[r]["overlap"]
                assert abs(row["E_CASCI"] - ef) < 1e-8, (r, row["E_CASCI"], ef)
        continue
    t = time.time()
    t_start = now()
    space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    pool = uccsd_pool(H.nocc, H.norb)
    cs = counts_summary(H.nocc, H.norb, H.orbsym)
    rec = {"R_PtCl": r, "point_group": H.extra["point_group"],
           "E_RHF": H.e_rhf_total, "E_ref": ref["E_ref"], "E_CISD": ref["E_CISD"],
           "E_CCSD": ref["E_CCSD"], "CCSD_converged": ref["CCSD_converged"],
           "E_CCSD(T)": ref["E_CCSD(T)"], "E_CASCI": ef,
           "c0_sq": ref["c0_ref_weight"], "T1_diag": ref["t1_diag"],
           "active_NO_occupations": H.extra["act_no_occ"].tolist(),
           "active_irreps": H.irrep_names, "counts": cs,
           "tracking_min_singular_value": hams[r]["overlap"]}
    for L in (1, 2):
        A = Ansatz(space, gva_groups(pool, H.orbsym, True, True), layers=L)
        rv = run_vqe(A, e_exact=ef, maxiter=2500, record=False)
        rec[f"GVA_{L}L"] = {"E": rv["E"], "dE": rv["dE"], "n_params": A.n_params,
                            "nit": rv["nit"], "S2": rv["S2"], "wall_s": rv["wall_s"],
                            "message": rv["message"], "theta": rv["theta"].tolist()}
    A = Ansatz(space, gva_groups(pool, H.orbsym, False, False), layers=1)
    rv = run_vqe(A, e_exact=ef, maxiter=2500, record=False)
    rec["UCCSD_1L"] = {"E": rv["E"], "dE": rv["dE"], "n_params": A.n_params,
                       "nit": rv["nit"], "S2": rv["S2"], "wall_s": rv["wall_s"],
                       "message": rv["message"], "theta": rv["theta"].tolist()}
    rec["wall_s"] = time.time() - t
    rec["started"], rec["finished"] = t_start, now()
    rows.append(rec)
    rows.sort(key=lambda x: x["R_PtCl"])
    json.dump(rows, open(OUT, "w"), indent=1, default=float)
    print(f"R={r:.2f} dCCSD={1e3*(ref['E_CCSD']-ef):.3f} dGVA2={1e3*rec['GVA_2L']['dE']:.4f} "
          f"c0^2={ref['c0_ref_weight']:.3f} t={rec['wall_s']:.0f}s", flush=True)
json.dump(rows, open(OUT, "w"), indent=1, default=float)
print("done")
