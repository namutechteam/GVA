import numpy as np, time
from gva.systems import small_molecule_hamiltonian, classical_references
from gva.qubit import *
for name in ["H2","LiH","H4","BeH2"]:
    H = small_molecule_hamiltonian(name,1.0); ref = classical_references(H)
    qh = qubit_hamiltonian(H)
    psi = civec_to_statevector(ref["civec"], H.norb, H.nelec)
    e,_ = expectation(qh, psi)
    print(name, "terms", len(qh["coef"]), "E_qubit-E_FCI = %.2e"%(e-ref["E_FCI"]), "const", qh["const"])
    for kind in ["QWC","GC"]:
        cols = group_terms(qh, kind)
        res = {k:(c.max()+1, shots_for_precision(group_sigmas(qh,psi,c),1.6e-3)) for k,c in cols.items()}
        print("  ",kind, {k:(v[0], "%.2e"%v[1]) for k,v in res.items()})
