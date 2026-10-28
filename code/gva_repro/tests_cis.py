import numpy as np, time, pickle
from gva.systems import *
t=time.time()
H = cisplatin_cas_hamiltonian()
mol = H.extra["mol"]
print("PG", mol.groupname, "topgroup", mol.topgroup, "nelec", mol.nelectron, "nao", mol.nao, "E_RHF", H.e_rhf_total, "t", time.time()-t)
print("MP2 corr", H.extra["e_mp2_corr"])
print("active NO occ", np.round(H.extra["act_no_occ"],4))
print("active irreps", H.irrep_names, "eps", np.round(H.extra["active_eps"],3))
print("ncore", H.extra["n_core"], "E_ref", H.extra["e_ref_total"])
ref = classical_references(H)
for k in ["E_ref","E_CISD","E_CCSD","E_CCSD(T)","E_FCI"]: print(k, ref[k], (ref[k]-ref["E_FCI"])*1e3)
print("T1", ref["t1_diag"], "c0^2", ref["c0_ref_weight"], "t", time.time()-t)
# Mulliken-like composition of active orbitals
from pyscf import lo
C = H.extra["C_act_ao"]; S = mol.intor("int1e_ovlp")
labels = mol.ao_labels()
for k in range(H.norb):
    pop = C[:,k]*(S@C[:,k])
    top = np.argsort(-np.abs(pop))[:4]
    print(k, H.irrep_names[k], [(labels[i].strip(), round(pop[i],2)) for i in top])
pickle.dump({"H":{k:getattr(H,k) for k in ["label","norb","nelec","ecore","h1","eri","orbsym","irrep_names","e_rhf_total"]},"ref":{k:v for k,v in ref.items()}, "C":C, "occ":H.extra["act_no_occ"]}, open("results/cis_eq_tmp.pkl","wb"))
