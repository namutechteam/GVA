"""Validation of the periodic code (gva/periodic.py). Exits non-zero on failure.

1. Integrals: momentum conservation, Hermiticity, time-reversal relations;
   Hartree-Fock energy from the integrals equals the KRHF energy x N_k.
2. Three routes to the FCI energy agree: real supercell basis (PySCF
   direct_spin1), Bloch basis (complex direct_nosym, Lanczos) and an explicit
   Gamma-point supercell built with pyscf.pbc.tools.super_cell.
3. Adjoint gradients of the complex ansatz agree with central finite
   differences for every tying variant, one and two layers.
4. Momentum-violating excitations have zero energy gradient on states built
   from momentum-conserving ones.
5. Jordan-Wigner energy of a complex Bloch-basis state equals the
   determinant-space energy (OpenFermion's internal 1e-8 term pruning is
   disabled in qubit_hamiltonian_c; with it, differences of ~1e-9 Ha appear
   for systems with many integrals below 1e-8).
6. Tying relations: in the untied k-only optimum (HF start) the amplitudes of
   time-reversal and spin-flip partners satisfy the imposed relations up to
   Trotter-ordering effects; the largest deviation is reported.
"""
import sys

import numpy as np
from pyscf import fci

from gva.ansatz import uccsd_pool, excitation_label
from gva.periodic import (periodic_hamiltonian, supercell_gamma_check, CDeterminantSpace,
                          CAnsatz, run_vqe_c, variant_orbits, VARIANTS_K, exc_momentum,
                          pool_gradients, qubit_hamiltonian_c, civec_to_statevector_c,
                          tr_partner, spin_partner)
from gva.qubit import expectation

fail = 0


def check(name, ok, detail=""):
    global fail
    fail += 0 if ok else 1
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)


for (r, d, nk) in [(0.74, 1.50, 3), (1.0, 1.0, 3), (0.74, 1.50, 4)]:
    H = periodic_hamiltonian(r, d, nk)
    tag = f"r={r} d={d} Nk={nk}"
    check(f"[{tag}] momentum conservation", H.extra["momentum_violation_max"] < 1e-12,
          f"max |violating integral| = {H.extra['momentum_violation_max']:.1e}")
    check(f"[{tag}] Hermiticity", H.extra["hermiticity_max"] < 1e-12,
          f"{H.extra['hermiticity_max']:.1e}")
    check(f"[{tag}] time reversal", H.extra["time_reversal_max"] < 1e-12,
          f"{H.extra['time_reversal_max']:.1e}")
    dhf = abs(H.extra["E_HF_from_integrals"] - H.extra["E_KRHF_supercell"])
    check(f"[{tag}] E_HF(integrals) = N_k E_KRHF", dhf < 1e-10, f"{dhf:.1e} Ha")
    e_r, _ = fci.direct_spin1.kernel(H.h1_real, H.eri_real, H.norb, H.nelec,
                                     ecore=H.ecore, conv_tol=1e-13)
    sp = CDeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    w, V = sp.ground_state(k=1)
    sc = supercell_gamma_check(r, d, nk)
    check(f"[{tag}] FCI real = Bloch", abs(w[0] - e_r) < 1e-10, f"{abs(w[0]-e_r):.1e} Ha")
    check(f"[{tag}] FCI real = explicit supercell", abs(sc["E_FCI_supercell"] - e_r) < 1e-9,
          f"{abs(sc['E_FCI_supercell']-e_r):.1e} Ha")
    pool = uccsd_pool(H.nocc, H.norb)
    rng = np.random.default_rng(7)
    for v in VARIANTS_K:
        orbs, _ = variant_orbits(H, v)
        for L in (1, 2):
            A = CAnsatz(sp, orbs, L)
            th = rng.normal(0, 0.2, A.n_params)
            _, g = A.energy_and_grad(th)
            h = 1e-5
            idx = rng.choice(A.n_params, min(8, A.n_params), replace=False)
            fd = np.array([(A.energy_and_grad(th + h * np.eye(A.n_params)[i])[0] -
                            A.energy_and_grad(th - h * np.eye(A.n_params)[i])[0]) / (2 * h)
                           for i in idx])
            err = np.abs(fd - g[idx]).max()
            check(f"[{tag}] gradient {v} {L}L", err < 1e-7, f"max |FD - adjoint| = {err:.1e}")
    orbs, _ = variant_orbits(H, "GVA-k")
    A = CAnsatz(sp, orbs, 1)
    psi = A.state(rng.normal(0, 0.3, A.n_params))
    forb = [e for e in pool if exc_momentum(e, H.kidx, nk) != 0]
    gf = pool_gradients(sp, psi, forb)
    check(f"[{tag}] zero gradient of momentum-violating excitations", gf.max() < 1e-12,
          f"{len(forb)} excitations, max {gf.max():.1e}")
    if H.norb <= 8:
        qh = qubit_hamiltonian_c(H.h1, H.eri, H.ecore)      # exact (no pruning)
        ejw, _ = expectation(qh, civec_to_statevector_c(psi, H.norb, H.nelec))
        check(f"[{tag}] Jordan-Wigner energy", abs(ejw - sp.energy(psi)) < 1e-10,
              f"{abs(ejw - sp.energy(psi)):.1e} Ha")
    # tying relations in the untied optimum
    orbs_k, _ = variant_orbits(H, "k only")
    Ak = CAnsatz(sp, orbs_k, 1)
    res = run_vqe_c(Ak)
    amp = {o.members[0][0]: res["theta"][2 * i] + 1j * res["theta"][2 * i + 1]
           for i, o in enumerate(orbs_k)}
    dev_tr = dev_sp = 0.0
    tmax = max(abs(t) for t in amp.values())
    for e, t in amp.items():
        f, s = tr_partner(e, H.partner)
        dev_tr = max(dev_tr, abs(amp[f] - s * np.conj(t)))
        f, s = spin_partner(e)
        dev_sp = max(dev_sp, abs(amp[f] - s * t))
    print(f"INFO  [{tag}] untied optimum: max |t| = {tmax:.3e}; max deviation from "
          f"t(TR partner) = s t* : {dev_tr:.1e}; from t(spin partner) = t : {dev_sp:.1e}")
    check(f"[{tag}] tying relations hold to < 5% of max |t|",
          dev_tr < 0.05 * tmax and dev_sp < 0.05 * tmax)

print("FAILURES:", fail)
sys.exit(1 if fail else 0)
