"""Jordan-Wigner energy of the cisplatin CASCI vector with and without
OpenFermion's internal 1e-8 term pruning (explains the 1e-8 Ha tolerance of
the original check in run_measure.py). Writes results/logs/check_jw_cisplatin.log."""
import numpy as np
from pyscf import fci

from gva.periodic import qubit_hamiltonian_c
from gva.qubit import qubit_hamiltonian, civec_to_statevector, expectation
from gva.systems import ActiveHamiltonian

d = np.load("results/cisplatin_cas10_10.npz")
H = ActiveHamiltonian(label="cisplatin", norb=10, nelec=tuple(int(x) for x in d["nelec"]),
                      ecore=float(d["ecore"]), h1=d["h1"], eri=d["eri"], orbsym=d["orbsym"],
                      irrep_names=[], e_rhf_total=np.nan)
e_cas = fci.direct_spin1.energy(H.h1, H.eri, d["civec"], 10, H.nelec) + H.ecore
psi = civec_to_statevector(d["civec"], 10, H.nelec)
lines = [f"E_CASCI (determinant space) = {e_cas:.12f} Ha"]
qh = qubit_hamiltonian(H)
e1, _ = expectation(qh, psi)
lines.append(f"OpenFermion default (pruning 1e-8): {len(qh['coef']) + 1} Pauli terms, "
             f"|E_JW - E_CASCI| = {abs(e1 - e_cas):.2e} Ha")
qh2 = qubit_hamiltonian_c(H.h1.astype(complex), H.eri.astype(complex), H.ecore)
e2, _ = expectation(qh2, psi)
lines.append(f"pruning disabled: {len(qh2['coef']) + 1} Pauli terms, "
             f"|E_JW - E_CASCI| = {abs(e2 - e_cas):.2e} Ha")
print("\n".join(lines))
open("results/logs/check_jw_cisplatin.log", "w").write("\n".join(lines) + "\n")
