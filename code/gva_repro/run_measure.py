"""Pauli grouping (QWC vs general commutativity) and shot estimates.

Shots are the optimal-allocation estimate N = (sum_g sigma_g)^2 / eps^2 for a
target standard error eps = 1.6 mHa on one energy evaluation, using the exact
ground state of each Hamiltonian (FCI or CASCI) as a proxy for the converged
variational state.
"""
import json
import time

import numpy as np

from gva.systems import (small_molecule_hamiltonian, classical_references,
                         ActiveHamiltonian)
from gva.qubit import (qubit_hamiltonian, civec_to_statevector, expectation,
                       group_terms, group_sigmas, shots_for_precision)

EPS = 1.6e-3
out = []


def analyse(label, H, civec, e_exact):
    t = time.time()
    qh = qubit_hamiltonian(H)
    psi = civec_to_statevector(civec, H.norb, H.nelec)
    e, vals = expectation(qh, psi)
    assert abs(e - e_exact) < 1e-8, (label, e, e_exact)
    sig_terms = np.abs(qh["coef"]) * np.sqrt(np.clip(1 - vals ** 2, 0, None))
    rec = {"system": label, "n_qubits": qh["n_qubits"],
           "n_pauli_terms": int(len(qh["coef"])) + 1,      # + identity
           "identity_coeff": qh["const"], "E_exact": e_exact,
           "shots_ungrouped": shots_for_precision(sig_terms, EPS), "QWC": {}, "GC": {}}
    for kind in ("QWC", "GC"):
        for meth, col in group_terms(qh, kind).items():
            sig = group_sigmas(qh, psi, col)
            rec[kind][meth] = {"n_groups": int(col.max() + 1),
                               "shots": shots_for_precision(sig, EPS)}
    rec["wall_s"] = time.time() - t
    out.append(rec)
    json.dump(out, open("results/measurement.json", "w"), indent=1, default=float)
    print(label, rec["n_pauli_terms"], {k: {m: (v["n_groups"], f"{v['shots']:.3g}")
                                           for m, v in rec[k].items()} for k in ("QWC", "GC")},
          f"{rec['wall_s']:.0f}s", flush=True)


for name in ("H2", "H4", "LiH", "BeH2"):
    H = small_molecule_hamiltonian(name, 1.0)
    ref = classical_references(H)
    analyse(name, H, ref["civec"], ref["E_FCI"])

d = np.load("results/cisplatin_cas10_10.npz")
H = ActiveHamiltonian(label="cisplatin", norb=10, nelec=tuple(d["nelec"]),
                      ecore=float(d["ecore"]), h1=d["h1"], eri=d["eri"],
                      orbsym=d["orbsym"], irrep_names=[], e_rhf_total=np.nan)
from pyscf import fci
e_cas = fci.direct_spin1.energy(H.h1, H.eri, d["civec"], 10, H.nelec) + H.ecore
analyse("Cisplatin CAS(10e,10o)", H, d["civec"], e_cas)
