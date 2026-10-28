"""Assemble the deliverables: checksums, reproduction package, log archive.

  ../deliverables/GVA_reproduction_package.zip  all code, inputs, results, logs, docs
  ../deliverables/GVA_calculation_logs.zip      every calculation log (v1, v2, v3)
  ../deliverables/<manuscript, SI, source data>
"""
import hashlib
import os
import shutil
import zipfile

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(ROOT), "deliverables")
BUILT = {"GVA_npjQI_Manuscript.docx", "GVA_npjQI_Supplementary_Information.docx",
         "GVA_Source_Data.xlsx"}
os.makedirs(OUT, exist_ok=True)


def files(base, skip_built=True):
    for d, dirs, fs in os.walk(base):
        dirs[:] = sorted(x for x in dirs if x != "__pycache__")
        for f in sorted(fs):
            if f.endswith(".pyc") or f == "CHECKSUMS.sha256":
                continue
            if skip_built and f in BUILT and d == base:
                continue
            yield os.path.join(d, f)


# figure files under the submitted (npj Quantum Information) numbering
import numbering
FIGDIR = os.path.join(ROOT, "figures_npj")
os.makedirs(FIGDIR, exist_ok=True)
for old, new in numbering.FILES.items():
    for ext in ("png", "pdf"):
        shutil.copy(os.path.join(ROOT, "figures", f"{old}.{ext}"), os.path.join(FIGDIR, f"{new}.{ext}"))

lines = []
for p in files(ROOT):
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    lines.append(f"{h}  ./{os.path.relpath(p, ROOT)}")
open(os.path.join(ROOT, "CHECKSUMS.sha256"), "w").write("\n".join(lines) + "\n")

pkg = os.path.join(OUT, "GVA_reproduction_package.zip")
with zipfile.ZipFile(pkg, "w", zipfile.ZIP_DEFLATED) as z:
    for p in list(files(ROOT)) + [os.path.join(ROOT, "CHECKSUMS.sha256")]:
        z.write(p, os.path.join("gva_repro", os.path.relpath(p, ROOT)))

readme = """GVA calculation logs
====================
results/logs/                 v2 molecular runs and all v3 periodic runs
  run_small.log, run_cisplatin.log, run_pes.log, run_orbital_sensitivity.log
  run_periodic_A.log ... run_periodic_G.log   periodic benchmark, seven parts
      (A, B, D stopped at the 2-h job limit; E stopped by hand during ADAPT-VQE;
       every completed run is in results/periodic.json; see PROVENANCE.md, section 9)
  run_scaling_periodic.log, run_measure_periodic.log
  tests_periodic.log, verify_periodic.log, verify_periodic.stdout,
  analyze_tying_periodic.log, check_jw_cisplatin.log, verify_theta.log,
  compare_v1_v2.txt, periodic_start_time.txt, v2_start_time.txt
  pyscf/       full PySCF output (verbose 4): molecular runs and the CISD/CCSD/
               CCSD(T)/FCI references of every periodic system (<system>__mol1.log)
  pyscf_pbc/   PySCF output of the periodic KRHF (<system>__krhf.log) and of the
               independent explicit-supercell check (<system>__supercell.log)
results_v1/logs/              original (v1) run logs, validation and reference logs
"""
logs = os.path.join(OUT, "GVA_calculation_logs.zip")
with zipfile.ZipFile(logs, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("GVA_logs/LOGS_README.txt", readme)
    for sub in ("results/logs", "results_v1/logs"):
        for p in files(os.path.join(ROOT, sub), skip_built=False):
            z.write(p, os.path.join("GVA_logs", os.path.relpath(p, ROOT)))

figzip = os.path.join(OUT, "GVA_npjQI_figures.zip")
with zipfile.ZipFile(figzip, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(os.listdir(FIGDIR)):
        z.write(os.path.join(FIGDIR, f), os.path.join("figures", f))

for f in BUILT:
    shutil.copy(os.path.join(ROOT, f), OUT)
shutil.copy(os.path.join(ROOT, "uccsd_counts.py"), OUT)
for f in sorted(os.listdir(OUT)):
    print(f"{f:45s} {os.path.getsize(os.path.join(OUT, f)) / 1e6:8.2f} MB")
print(len(lines), "files checksummed")
