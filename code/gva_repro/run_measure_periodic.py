"""Measurement cost of periodic Hamiltonians in the Bloch and real supercell bases.

For the 20-qubit chains (N_k = 5) the Jordan-Wigner Hamiltonian is built in
both bases with Pauli terms below 1e-8 removed (the threshold OpenFermion
applies by default, as in the molecular counts). Terms are grouped by sorted
insertion under qubit-wise (QWC) and general (GC) commutativity, and the shot
count for a 1.6 mHa standard error is the optimal-allocation estimate with the
exact ground state as the trial state. OpenFermion's default pruning is kept while
the operator is built (as for the molecular counts). Writes
results/measurement_periodic.json.
"""
import json
import time

import numpy as np

from gva.periodic import (periodic_hamiltonian, CDeterminantSpace, qubit_hamiltonian_c,
                          civec_to_statevector_c)
from gva.qubit import (expectation, conflict_matrix, sorted_insertion, group_sigmas,
                       shots_for_precision, civec_to_statevector)
from pyscf import fci

EPS = 1.6e-3
out = []
for (r, d, nk) in [(0.74, 1.50, 5), (1.0, 1.0, 5)]:
    H = periodic_hamiltonian(r, d, nk)
    e_r, c_r = fci.direct_spin1.kernel(H.h1_real, H.eri_real, H.norb, H.nelec,
                                       ecore=H.ecore, conv_tol=1e-12)
    sp = CDeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    w, V = sp.ground_state(k=1)
    for basis in ("Bloch", "real"):
        t = time.time()
        if basis == "Bloch":
            qh = qubit_hamiltonian_c(H.h1, H.eri, H.ecore, prune=1e-8, exact=False)
            psi = civec_to_statevector_c(V[:, 0], H.norb, H.nelec)
            e_ex = float(w[0])
        else:
            qh = qubit_hamiltonian_c(H.h1_real.astype(complex), H.eri_real.astype(complex),
                                     H.ecore, prune=1e-8, exact=False)
            psi = civec_to_statevector(c_r, H.norb, H.nelec)
            e_ex = float(e_r)
        e, vals = expectation(qh, psi)
        sig_terms = np.abs(qh["coef"]) * np.sqrt(np.clip(1 - vals ** 2, 0, None))
        rec = {"system": f"H2chain_r{r:.2f}_d{d:.2f}_Nk{nk}", "basis": basis,
               "n_qubits": qh["n_qubits"], "n_pauli_terms": int(len(qh["coef"])) + 1,
               "one_norm": float(np.abs(qh["coef"]).sum()), "E_exact": e_ex,
               "E_pauli": float(e), "abs_diff_pruned": float(abs(e - e_ex)),
               "shots_ungrouped": shots_for_precision(sig_terms, EPS)}
        for kind in ("QWC", "GC"):
            C = conflict_matrix(qh, kind)
            col = sorted_insertion(qh, C)
            sig = group_sigmas(qh, psi, col)
            rec[kind] = {"n_groups": int(col.max() + 1), "shots": shots_for_precision(sig, EPS)}
        rec["wall_s"] = time.time() - t
        out.append(rec)
        json.dump(out, open("results/measurement_periodic.json", "w"), indent=1)
        print(rec["system"], basis, rec["n_pauli_terms"], f"1-norm {rec['one_norm']:.3f}",
              f"|E_pauli-E| {rec['abs_diff_pruned']:.1e}",
              {k: (rec[k]["n_groups"], f"{rec[k]['shots']:.3g}") for k in ("QWC", "GC")},
              f"{rec['wall_s']:.0f}s", flush=True)
