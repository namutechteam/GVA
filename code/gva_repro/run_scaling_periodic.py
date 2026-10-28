"""Parameter and CNOT counts for crystals versus k-point mesh.

1. Closed-form/grouped counts (gva.periodic.periodic_counts) for
   1D meshes N_k = 1..32 and 3D meshes n x n x n (n = 1..5) with
   o = v = 1, 2, 4 bands per cell, and for a molecular-crystal example
   (four molecules per cell, two occupied and two virtual orbitals each).
2. Check of the closed forms against explicit enumeration of the pool and of
   the spin-flip/time-reversal orbits (build_orbits) on small 1D, 2D and 3D
   meshes, including meshes with even sizes (time-reversal-invariant k-points
   other than Gamma).
Writes results/scaling_periodic.json.
"""
import json
import time
from datetime import datetime, timezone

from gva.ansatz import uccsd_pool, cnot_count_ladder
from gva.periodic import (periodic_counts, synthetic_mesh_structure, build_orbits,
                          orbit_counts)


def now():
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


t0 = time.time()
out = {"one_d": [], "three_d": [], "molecular_crystal": [], "enumeration_check": []}
for ov in (1, 2):
    for nk in range(1, 33):
        out["one_d"].append(periodic_counts((nk,), ov, ov))
for ov in (1, 2, 4):
    for n in range(1, 6):
        out["three_d"].append(periodic_counts((n, n, n), ov, ov))
for n in (1, 2, 3, 4):
    out["molecular_crystal"].append(periodic_counts((n, n, n), 8, 8))
print(now(), f"closed-form tables done ({time.time() - t0:.0f}s)", flush=True)

checks = [((2,), 1, 1), ((3,), 1, 1), ((4,), 1, 1), ((5,), 1, 1), ((6,), 1, 1),
          ((3,), 1, 2), ((4,), 2, 1), ((2, 2, 1), 1, 1), ((3, 3, 1), 1, 1),
          ((4, 2, 1), 1, 1), ((2, 2, 2), 1, 1), ((4, 4, 1), 1, 1)]
for mesh, o, v in checks:
    t = time.time()
    pc = periodic_counts(mesh, o, v)
    kidx, partner, allowed, nocc, norb = synthetic_mesh_structure(mesh, o, v)
    pool = uccsd_pool(nocc, norb)
    rec = {"mesh": list(mesh), "o": o, "v": v, "pool": len(pool),
           "UCCSD_params_enum": len(pool),
           "UCCSD_cnot_enum": sum(cnot_count_ladder(e) for e in pool)}
    for name, flags in (("k only", (True, False, False)), ("k + spin", (True, True, False)),
                        ("k + TR", (True, False, True)), ("GVA-k", (True, True, True))):
        orbs, nzero = build_orbits(pool, kidx, None, partner, *flags, allowed=allowed)
        c = orbit_counts(orbs)
        rec[name] = {"n_params": c["n_params"], "cnot_1L": c["cnot_1L"],
                     "n_generators": c["n_generators"], "n_zero_orbits": nzero,
                     "n_real_orbits": c["n_real_orbits"]}
    rec["closed_form"] = {k: pc[k] for k in ("UCCSD_params", "UCCSD_cnot_1L", "k_only_params",
                                             "k_spin_params", "k_TR_params", "GVA_k_params",
                                             "k_cnot_1L", "k_generators")}
    rec["agree"] = bool(
        rec["UCCSD_params_enum"] == pc["UCCSD_params"] and
        rec["UCCSD_cnot_enum"] == pc["UCCSD_cnot_1L"] and
        rec["k only"]["n_params"] == pc["k_only_params"] and
        rec["k + spin"]["n_params"] == pc["k_spin_params"] and
        rec["k + TR"]["n_params"] == pc["k_TR_params"] and
        rec["GVA-k"]["n_params"] == pc["GVA_k_params"] and
        rec["GVA-k"]["cnot_1L"] == pc["k_cnot_1L"] and
        rec["GVA-k"]["n_generators"] == pc["k_generators"] and
        rec["GVA-k"]["n_zero_orbits"] == 0)
    rec["wall_s"] = time.time() - t
    out["enumeration_check"].append(rec)
    print(now(), mesh, o, v, "pool", len(pool), "GVA-k", rec["GVA-k"]["n_params"],
          "closed", pc["GVA_k_params"], "agree", rec["agree"], f"{rec['wall_s']:.0f}s",
          flush=True)
    json.dump(out, open("results/scaling_periodic.json", "w"), indent=1)
assert all(r["agree"] for r in out["enumeration_check"])
print(now(), "all closed forms agree with enumeration", flush=True)
