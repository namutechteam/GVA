"""Compare the v2 rerun (results/) with the original v1 results (results_v1/).

Runs are matched by (system, geometry, ansatz, layers, seed). Reports the
largest energy difference, differences in iteration counts and wall times,
and any run present in only one version. Writes results/logs/compare_v1_v2.txt.
"""
import json

import numpy as np

out = []


def compare(label, v1, v2, key):
    a = {key(r): r for r in v1}
    b = {key(r): r for r in v2}
    only1, only2 = sorted(set(a) - set(b), key=str), sorted(set(b) - set(a), key=str)
    common = [k for k in a if k in b]
    dE = np.array([abs(a[k]["E"] - b[k]["E"]) for k in common])
    dit = [(k, a[k].get("nit"), b[k].get("nit")) for k in common
           if a[k].get("nit") != b[k].get("nit")]
    out.append(f"== {label}: {len(common)} matched runs")
    out.append(f"   max |E_v2 - E_v1| = {dE.max():.2e} Ha, median {np.median(dE):.2e} Ha")
    worst = sorted(zip(dE, common), key=lambda x: -x[0])[:5]
    for d, k in worst:
        out.append(f"   {str(k):55s} |dE| = {d:.2e} Ha")
    out.append(f"   runs with a different iteration count: {len(dit)}")
    for k, i1, i2 in dit[:10]:
        out.append(f"     {str(k):53s} v1 {i1} -> v2 {i2}")
    if only1:
        out.append(f"   only in v1: {only1}")
    if only2:
        out.append(f"   only in v2: {only2}")
    return dE.max()


s1 = json.load(open("results_v1/small_molecules.json"))["runs"]
s2 = json.load(open("results/small_molecules.json"))["runs"]
m1 = compare("small molecules", s1, s2,
             lambda r: (r["system"], r["scale"], r["ansatz"], r["layers"], r["seed"]))
c1 = json.load(open("results_v1/cisplatin_eq.json"))["runs"]
c2 = json.load(open("results/cisplatin_eq.json"))["runs"]
m2 = compare("cisplatin C2v", c1, c2, lambda r: (r["ansatz"], r["layers"], r["seed"]))


def flat(p):
    rows = []
    for r in p:
        for k in ("GVA_1L", "GVA_2L", "UCCSD_1L"):
            rows.append({"R": round(r["R_PtCl"], 3), "k": k, "E": r[k]["E"], "nit": r[k]["nit"]})
        rows.append({"R": round(r["R_PtCl"], 3), "k": "CASCI", "E": r["E_CASCI"], "nit": None})
    return rows


p1 = flat(json.load(open("results_v1/pes.json")))
p2 = flat(json.load(open("results/pes.json")))
m3 = compare("Pt-Cl scan", p1, p2, lambda r: (r["R"], r["k"]))
txt = "\n".join(out)
open("results/logs/compare_v1_v2.txt", "w").write(txt + "\n")
print(txt)
