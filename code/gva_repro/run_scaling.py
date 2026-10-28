"""Parameter and CNOT counts versus active-space size for cisplatin.

Active spaces CAS(k e, k o), k = 4..20, are taken from the same symmetry-
adapted MP2 natural orbitals (k/2 strongly and k/2 weakly occupied NOs).
"""
import json

from gva.systems import (build_cisplatin_mol, run_rhf, mp2_natural_orbitals,
                         select_active_by_occupation, cas_hamiltonian_from_orbitals)
from gva.ansatz import counts_summary

mol = build_cisplatin_mol()
mf = run_rhf(mol)
nos, occs, syms, _ = mp2_natural_orbitals(mol, mf)
rows = []
for k in range(4, 22, 2):
    act = select_active_by_occupation(occs, k // 2, k // 2)
    H = cas_hamiltonian_from_orbitals(mol, mf, nos, occs, syms, act, f"CAS({k},{k})")
    cs = counts_summary(H.nocc, H.norb, H.orbsym)
    rows.append({"n_active_orbitals": k, "n_qubits": 2 * k,
                 "UCCSD_params": cs["UCCSD"]["n_params"],
                 "UCCSD_formula": cs["formula_UCCSD"],
                 "UCCSD_cnot_1L": cs["UCCSD"]["cnot_1L"],
                 "PGonly_params": cs["PG only"]["n_params"],
                 "spinonly_params": cs["spin only"]["n_params"],
                 "GVA_params": cs["GVA"]["n_params"],
                 "GVA_cnot_1L": cs["GVA"]["cnot_1L"],
                 "irreps": H.irrep_names})
    print(rows[-1], flush=True)
json.dump(rows, open("results/scaling.json", "w"), indent=1)
