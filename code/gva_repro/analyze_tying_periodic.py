"""How closely the untied optimum obeys the tying relations (periodic systems).

For every system, the single-layer 'k only' run from the Hartree-Fock start
(each momentum-conserving excitation with its own complex amplitude) is read
from results/periodic.json. For each excitation the amplitude of its
time-reversal partner is compared with s t^* and that of its spin-flip partner
with t. Writes results/tying_check_periodic.json.
"""
import json

import numpy as np

from gva.ansatz import uccsd_pool
from gva.periodic import tr_partner, spin_partner, build_orbits

per = json.load(open("results/periodic.json"))
out = []
for s in per["systems"]:
    z = np.load(f"results/periodic_integrals_{s['system']}.npz")
    nocc, norb = int(z["nelec"][0]), z["h1"].shape[0]
    pool = uccsd_pool(nocc, norb)
    orbs, _ = build_orbits(pool, z["kidx"], s["n_k"], z["partner"], True, False, False)
    r = [x for x in per["runs"] if x["system"] == s["system"] and x["ansatz"] == "k only"
         and x["layers"] == 1 and x["seed"] == 0][0]
    th = np.asarray(r["theta"])
    amp = {o.members[0][0]: th[2 * i] + 1j * th[2 * i + 1] for i, o in enumerate(orbs)}
    tmax = max(abs(t) for t in amp.values())
    dtr = max(abs(amp[tr_partner(e, z["partner"])[0]] - tr_partner(e, z["partner"])[1] * np.conj(t))
              for e, t in amp.items())
    dsp = max(abs(amp[spin_partner(e)[0]] - t) for e, t in amp.items())
    gk = [x for x in per["runs"] if x["system"] == s["system"] and x["ansatz"] == "GVA-k"
          and x["layers"] == 1 and x["seed"] == 0][0]
    rec = {"system": s["system"], "n_amplitudes": len(amp), "max_abs_t": tmax,
           "max_dev_time_reversal": dtr, "max_dev_spin_flip": dsp,
           "rel_dev_time_reversal": dtr / tmax, "rel_dev_spin_flip": dsp / tmax,
           "E_k_only_minus_E_GVA_k_mHa": (r["E"] - gk["E"]) * 1e3}
    out.append(rec)
    print(f"{s['system']:28s} max|t| {tmax:.3e}  TR dev {dtr:.1e} ({dtr/tmax:.1%})  "
          f"spin dev {dsp:.1e} ({dsp/tmax:.1%})  E(k only) - E(GVA-k) "
          f"{rec['E_k_only_minus_E_GVA_k_mHa']:+.4f} mHa")
json.dump(out, open("results/tying_check_periodic.json", "w"), indent=1)
