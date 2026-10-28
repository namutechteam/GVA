"""Rebuild every periodic VQE / ADAPT-VQE state from the stored parameters.

For each run in results/periodic.json the ansatz is reconstructed from its
definition (variant and layers, or for ADAPT-VQE the stored orbit labels), the
integrals are read from results/periodic_integrals_<system>.npz, and energy
and <S^2> are compared with the stored values (tolerances 1e-9 Ha, 1e-8).
For the Hartree-Fock-start single-layer UCCSD and GVA-k states the energy is
also evaluated with the exact Jordan-Wigner qubit Hamiltonian (OpenFermion
pruning disabled) for systems with up to 16 qubits. Writes results/verify_periodic.json and
results/logs/verify_periodic.log; exits non-zero on any mismatch.
"""
import json
import sys

import numpy as np

from gva.ansatz import uccsd_pool, excitation_label
from gva.engine import DeterminantSpace, Ansatz
from gva.periodic import (CDeterminantSpace, CAnsatz, Orbit, build_orbits, VARIANTS_K,
                          qubit_hamiltonian_c, civec_to_statevector_c)
from gva.qubit import expectation, civec_to_statevector

TOL_E, TOL_S2 = 1e-9, 1e-8
data = json.load(open("results/periodic.json"))
sysd = {s["system"]: s for s in data["systems"]}
rows, lines, bad = [], [], 0


def orbits_from_labels(pool, labels):
    lut = {excitation_label(e): e for e in pool}
    out = []
    for o in labels:
        out.append(Orbit([(lut[l], int(s), bool(c)) for l, s, c in o["members"]], o["mode"]))
    return out


cache = {}
for r in data["runs"]:
    tag = r["system"]
    if tag not in cache:
        z = np.load(f"results/periodic_integrals_{tag}.npz")
        s = sysd[tag]
        nelec = tuple(int(x) for x in z["nelec"])
        norb = z["h1"].shape[0]
        csp = CDeterminantSpace(norb, nelec, z["h1"], z["eri"], float(z["ecore"]))
        rsp = DeterminantSpace(norb, nelec, z["h1_real"], z["eri_real"], float(z["ecore"]))
        pool = uccsd_pool(nelec[0], norb)
        orbs = {v: build_orbits(pool, z["kidx"], s["n_k"], z["partner"], *f)[0]
                for v, f in VARIANTS_K.items()}
        cache[tag] = (z, csp, rsp, pool, orbs, norb, nelec)
        print(tag, "loaded", flush=True)
    z, csp, rsp, pool, orbs, norb, nelec = cache[tag]
    th = np.asarray(r["theta"])
    if r["ansatz"] == "UCCSD":
        A = Ansatz(rsp, [(e,) for e in pool], layers=r["layers"])
        sp = rsp
    elif r["ansatz"] == "ADAPT-VQE":
        A = CAnsatz(csp, orbits_from_labels(pool, r["parameter_labels"]), 1)
        sp = csp
    else:
        A = CAnsatz(csp, orbs[r["ansatz"]], layers=r["layers"])
        sp = csp
    assert A.n_params == len(th), (tag, r["ansatz"], A.n_params, len(th))
    psi = A.state(th)
    E = sp.energy(psi)
    s2 = sp.s2(psi)
    dE, dS = abs(E - r["E"]), abs(s2 - r["S2"])
    rec = {"system": tag, "ansatz": r["ansatz"], "layers": r["layers"], "seed": r["seed"],
           "n_params": A.n_params, "E_rebuilt": E, "E_stored": r["E"], "abs_dE": dE,
           "S2_rebuilt": s2, "S2_stored": r["S2"], "abs_dS2": dS}
    ok = dE < TOL_E and dS < TOL_S2
    # exact JW check up to 16 qubits (building the unpruned 20-qubit operator with
    # OpenFermion takes >10 min each); 20-qubit states carry the run-time JW check
    if r["ansatz"] in ("UCCSD", "GVA-k") and r["layers"] == 1 and r["seed"] == 0 and norb <= 8:
        if r["ansatz"] == "GVA-k":
            qh = qubit_hamiltonian_c(z["h1"], z["eri"], float(z["ecore"]))
            sv = civec_to_statevector_c(psi, norb, nelec)
        else:
            qh = qubit_hamiltonian_c(z["h1_real"].astype(complex), z["eri_real"].astype(complex),
                                     float(z["ecore"]))
            sv = civec_to_statevector(np.asarray(psi).reshape(rsp.sa, rsp.sb), norb, nelec)
        ejw, _ = expectation(qh, sv)
        rec["E_jordan_wigner_exact"] = float(ejw)
        rec["jw_abs_diff"] = float(abs(ejw - E))
        rec["n_pauli_terms_exact"] = int(len(qh["coef"])) + 1
        ok = ok and rec["jw_abs_diff"] < 1e-9
    rec["ok"] = bool(ok)
    bad += 0 if ok else 1
    rows.append(rec)
    line = (f"{tag:28s} {r['ansatz']:9s} {str(r['layers']):4s} s{r['seed']} n_par={A.n_params:4d} "
            f"|dE|={dE:.1e} |dS2|={dS:.1e}"
            + (f" |E_JW-E|={rec['jw_abs_diff']:.1e} ({rec['n_pauli_terms_exact']} terms)"
               if "jw_abs_diff" in rec else "") + ("  OK" if ok else "  MISMATCH"))
    lines.append(line)
    print(line, flush=True)

summary = (f"{len(rows)} states rebuilt; max |dE| = {max(x['abs_dE'] for x in rows):.1e} Ha; "
           f"max |dS2| = {max(x['abs_dS2'] for x in rows):.1e}; "
           f"max |E_JW - E| = {max(x.get('jw_abs_diff', 0) for x in rows):.1e} Ha; mismatches: {bad}")
print(summary)
json.dump({"summary": summary, "runs": rows}, open("results/verify_periodic.json", "w"), indent=1)
open("results/logs/verify_periodic.log", "w").write("\n".join(lines + [summary]) + "\n")
sys.exit(1 if bad else 0)
