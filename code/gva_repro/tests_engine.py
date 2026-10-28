import numpy as np, time
from gva.systems import small_molecule_hamiltonian, classical_references
from gva.ansatz import uccsd_pool, gva_groups, counts_summary
from gva.engine import DeterminantSpace, Ansatz, run_vqe

for name, sc in [("H2",1.0),("LiH",1.0),("H4",1.0),("BeH2",2.0)]:
    H = small_molecule_hamiltonian(name, sc)
    ref = classical_references(H)
    sp_ = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    assert abs(sp_.energy(sp_.reference()) - ref["E_ref"]) < 1e-9, (sp_.energy(sp_.reference()), ref["E_ref"])
    assert abs(ref["E_ref"]-H.e_rhf_total) < 1e-8
    pool = uccsd_pool(H.nocc, H.norb)
    cs = counts_summary(H.nocc, H.norb, H.orbsym)
    print(name, sc, H.irrep_names, {k:v for k,v in cs.items()})
    # gradient check
    for groups,lab in [(gva_groups(pool,H.orbsym,False,False),'UCCSD'),(gva_groups(pool,H.orbsym,True,True),'GVA')]:
        A = Ansatz(sp_, groups, layers=2)
        th = np.random.default_rng(1).normal(0,0.1,A.n_params)
        e,g = A.energy_and_grad(th)
        k = np.argmax(np.abs(g)); h=1e-6; tp=th.copy(); tp[k]+=h; tm=th.copy(); tm[k]-=h
        fd = (A.energy_and_grad(tp)[0]-A.energy_and_grad(tm)[0])/(2*h)
        r = run_vqe(Ansatz(sp_, groups, layers=1), e_exact=ref["E_FCI"])
        r2 = run_vqe(Ansatz(sp_, groups, layers=2), e_exact=ref["E_FCI"])
        print(f"  {lab}: grad {g[k]:.8f} fd {fd:.8f} | 1L dE={r['dE']*1e3:.5f} mHa nit={r['nit']} S2={r['S2']:.2e} | 2L dE={r2['dE']*1e3:.6f}")
    print("  refs dE(mHa): CISD %.4f CCSD %.4f CCSD(T) %.4f  Ecorr %.3f"%tuple([(ref[k]-ref['E_FCI'])*1e3 for k in ['E_CISD','E_CCSD','E_CCSD(T)','E_ref']]))
