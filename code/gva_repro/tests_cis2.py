import numpy as np, time, pickle
from gva.systems import ActiveHamiltonian
from gva.ansatz import *
from gva.engine import *
d = pickle.load(open("results/cis_eq_tmp.pkl","rb"))
H = ActiveHamiltonian(**d["H"], extra={})
sp_ = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
print("dim", sp_.dim, "Eref", sp_.energy(sp_.reference()))
print(counts_summary(H.nocc, H.norb, H.orbsym))
pool = uccsd_pool(H.nocc,H.norb)
for lab,(pg,spn) in {"GVA":(True,True),"UCCSD":(False,False)}.items():
    t=time.time(); A = Ansatz(sp_, gva_groups(pool,H.orbsym,pg,spn), layers=1); t1=time.time()-t
    th = np.zeros(A.n_params); t=time.time(); 
    for _ in range(3): A.energy_and_grad(th+0.01)
    print(lab, A.n_params, len(A.ops), "build %.1fs"%t1, "E+grad %.2fs"%((time.time()-t)/3))
