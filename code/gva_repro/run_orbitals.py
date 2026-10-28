"""Irreps, occupations and composition of the cisplatin active MP2 natural orbitals."""
import json
import numpy as np
from gva.systems import (build_cisplatin_mol, run_rhf, mp2_natural_orbitals,
                         select_active_by_occupation)

mol = build_cisplatin_mol()
mf = run_rhf(mol)
nos, occs, syms, _ = mp2_natural_orbitals(mol, mf)
act = select_active_by_occupation(occs, 5, 5)
names = dict(zip(mol.irrep_id, mol.irrep_name))
S = mol.intor("int1e_ovlp")
labels = mol.ao_labels(fmt=False)   # (atom index, element, shell label, m)
out = []
for k in act:
    c = nos[:, k]
    pop = c * (S @ c)
    # aggregate by atom + shell type (e.g. 'Pt 5d')
    # aggregate gross populations by element and angular momentum, summing
    # over symmetry-equivalent atoms (shell numbers of contracted functions
    # carry no physical meaning, so only s/p/d/f is kept)
    merged = {}
    for i, (ia, el, nl, m) in enumerate(labels):
        key = f"{el} {nl[-1]}"
        merged[key] = merged.get(key, 0.0) + pop[i]
    top = sorted(merged.items(), key=lambda x: -abs(x[1]))[:4]
    out.append({"irrep": names[syms[k]], "occupation": float(occs[k]),
                "composition": ", ".join(f"{k_} {v:.2f}" for k_, v in top)})
json.dump({"orbitals": out}, open("results/cisplatin_orbitals.json", "w"), indent=1)
for o in out:
    print(o)
