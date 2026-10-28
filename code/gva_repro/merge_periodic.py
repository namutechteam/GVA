"""Merge the partial outputs of run_periodic.py (two parallel processes) into
results/periodic.json, ordered as in run_periodic.SYSTEMS."""
import glob
import json

order = [f"H2chain_r{r:.2f}_d{d:.2f}_Nk{nk}" for (r, d, nk) in
         [(0.74, 1.50, nk) for nk in (2, 3, 4, 5)] +
         [(R, R, nk) for R in (1.0, 1.5, 2.0) for nk in (3, 5)]]
systems, runs, parts = {}, [], sorted(glob.glob("results/periodic_part*.json"))
for p in parts:
    d = json.load(open(p))
    for s in d["systems"]:
        systems.setdefault(s["system"], s)
    runs += d["runs"]
missing = [t for t in order if t not in systems]
assert not missing, missing
runs.sort(key=lambda r: (order.index(r["system"]), r["finished"]))
json.dump({"systems": [systems[t] for t in order], "runs": runs, "merged_from": parts},
          open("results/periodic.json", "w"), indent=1)
print(f"merged {len(parts)} files: {len(systems)} systems, {len(runs)} runs")
