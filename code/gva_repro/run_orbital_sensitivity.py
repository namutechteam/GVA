"""Dependence of the cisplatin CAS(10e,10o) results on the choice of active orbitals.

Three choices at the same C2v geometry, basis and ECP:
  A  symmetry-adapted MP2 natural orbitals (as in the main text): the five
     strongly occupied NOs with the lowest occupations and the five weakly
     occupied NOs with the highest occupations;
  B  canonical RHF orbitals: HOMO-4 ... HOMO and LUMO ... LUMO+4;
  C  CASSCF(10e,10o) orbitals optimized from guess A (state-specific ground
     state), expressed as active natural orbitals.
For each: correlation energy in the active space, reference weight, classical
references, parameter/CNOT counts and VQE errors (GVA 1L, 2L; UCCSD 1L; HF start).
Writes results/orbital_sensitivity.json and PySCF logs.
"""
import json
import time
from datetime import datetime, timezone

import numpy as np
from pyscf import mcscf, symm

from gva.systems import (build_cisplatin_mol, run_rhf, mp2_natural_orbitals,
                         select_active_by_occupation, cas_hamiltonian_from_orbitals,
                         classical_references, set_pyscf_log)
from gva.ansatz import uccsd_pool, gva_groups, counts_summary
from gva.engine import DeterminantSpace, Ansatz, run_vqe

OUT = "results/orbital_sensitivity.json"


def now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def composition(mol, C):
    S = mol.intor("int1e_ovlp")
    lab = mol.ao_labels(fmt=False)
    out = []
    for k in range(C.shape[1]):
        c = C[:, k]
        pop = c * (S @ c)
        agg = {}
        for i, (ia, el, nl, m) in enumerate(lab):
            key = f"{el} {nl[-1]}"
            agg[key] = agg.get(key, 0.0) + pop[i]
        top = sorted(agg.items(), key=lambda x: -abs(x[1]))[:3]
        out.append(", ".join(f"{a} {v:.2f}" for a, v in top))
    return out


def evaluate(label, H, mol, extra):
    t0 = time.time()
    ref = classical_references(H)
    ef = ref["E_FCI"]
    space = DeterminantSpace(H.norb, H.nelec, H.h1, H.eri, H.ecore)
    pool = uccsd_pool(H.nocc, H.norb)
    cs = counts_summary(H.nocc, H.norb, H.orbsym)
    rec = {"choice": label, "E_ref": ref["E_ref"], "E_CASCI": ef,
           "E_corr_active_mHa": (ef - ref["E_ref"]) * 1e3, "c0_sq": ref["c0_ref_weight"],
           "CISD_err_mHa": (ref["E_CISD"] - ef) * 1e3, "CCSD_err_mHa": (ref["E_CCSD"] - ef) * 1e3,
           "CCSD(T)_err_mHa": (ref["E_CCSD(T)"] - ef) * 1e3,
           "active_irreps": H.irrep_names,
           "active_composition": composition(mol, H.extra["C_act_ao"]),
           "counts": cs, "started": now()}
    rec.update(extra)
    for key, var, L in (("GVA_1L", (True, True), 1), ("UCCSD_1L", (False, False), 1),
                        ("GVA_2L", (True, True), 2)):
        A = Ansatz(space, gva_groups(pool, H.orbsym, *var), layers=L)
        rv = run_vqe(A, e_exact=ef, maxiter=2500, record=False)
        rec[key] = {"E": rv["E"], "dE_mHa": rv["dE"] * 1e3, "n_params": A.n_params,
                    "nit": rv["nit"], "S2": rv["S2"], "wall_s": rv["wall_s"],
                    "theta": rv["theta"].tolist()}
        print(label, key, f"{rv['dE'] * 1e3:.4f} mHa", flush=True)
    rec["finished"] = now()
    rec["wall_s"] = time.time() - t0
    return rec


set_pyscf_log("results/logs/pyscf", "orbital_sensitivity")
mol = build_cisplatin_mol()
mf = run_rhf(mol)
irr = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff))
nocc = mol.nelectron // 2
rows = []

# ---- A: MP2 natural orbitals
nos, occs, syms, e_mp2 = mp2_natural_orbitals(mol, mf)
actA = select_active_by_occupation(occs, 5, 5)
HA = cas_hamiltonian_from_orbitals(mol, mf, nos, occs, syms, actA, "A")
rows.append(evaluate("A: MP2 natural orbitals", HA, mol,
                     {"selected": f"NOs ranked {actA.min() + 1}-{actA.max() + 1} by occupation "
                                  f"(of {len(occs)}); occupations "
                                  + ", ".join(f"{o:.4f}" for o in occs[actA])}))
json.dump(rows, open(OUT, "w"), indent=1, default=float)

# ---- B: canonical RHF frontier orbitals
actB = np.arange(nocc - 5, nocc + 5)
occB = mf.mo_occ.copy()
HB = cas_hamiltonian_from_orbitals(mol, mf, mf.mo_coeff, occB, irr, actB, "B")
rows.append(evaluate("B: canonical RHF (HOMO-4 to LUMO+4)", HB, mol,
                     {"selected": f"canonical MOs {nocc - 4}-{nocc + 5} (HOMO = {nocc}); "
                                  "orbital energies (Ha) "
                                  + ", ".join(f"{e:.3f}" for e in mf.mo_energy[actB])}))
json.dump(rows, open(OUT, "w"), indent=1, default=float)

# ---- C: CASSCF(10e,10o) from guess A
ncore = HA.extra["n_core"]
core = np.array([k for k in np.where(occs > 1.0)[0] if k not in set(actA)])
virt = np.array([k for k in range(len(occs)) if k not in set(core) | set(actA)])
order = np.concatenate([core, actA, virt])
mc = mcscf.CASSCF(mf, 10, 10)
mc.natorb = True
mc.conv_tol = 1e-9
mc.max_cycle_macro = 100
t = time.time()
mc.kernel(nos[:, order])
mo = mc.mo_coeff
irrC = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mo))
occC = np.zeros(mo.shape[1])
occC[:ncore] = 2.0
nat = np.asarray(mc.mo_occ)[ncore:ncore + 10] if mc.mo_occ is not None else None
occC[ncore:ncore + 10] = nat if nat is not None else np.r_[np.full(5, 1.99), np.full(5, 0.01)]
actC = np.arange(ncore, ncore + 10)
HC = cas_hamiltonian_from_orbitals(mol, mf, mo, occC, irrC, actC, "C")
rows.append(evaluate("C: CASSCF(10e,10o) natural orbitals", HC, mol,
                     {"selected": f"CASSCF converged={bool(mc.converged)}, "
                                  f"E_CASSCF = {mc.e_tot:.8f} Ha, {time.time() - t:.0f} s; "
                                  "natural occupations "
                                  + ", ".join(f"{o:.4f}" for o in occC[actC]),
                      "E_CASSCF": float(mc.e_tot)}))
json.dump(rows, open(OUT, "w"), indent=1, default=float)
print("done", flush=True)
