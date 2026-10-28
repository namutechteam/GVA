# gva_repro — reproduction package

Scripts that regenerate every number, table and figure in
"Symmetry-screened unitary coupled-cluster circuits for molecules and crystals"
and its Supplementary Information.

## Requirements

Python ≥ 3.11, NumPy ≥ 2.0 (uses `np.bitwise_count`), SciPy, PySCF 2.14,
OpenFermion 1.8, matplotlib, openpyxl, python-docx.

```
pip install pyscf==2.14.0 openfermion==1.8.1 numpy scipy matplotlib openpyxl python-docx
```

## Layout

| Path | Content |
|---|---|
| `gva/systems.py` | Geometries (incl. idealised C2v cisplatin), RHF, symmetry-adapted MP2 natural orbitals, CAS Hamiltonians, CISD/CCSD/CCSD(T)/FCI references |
| `gva/ansatz.py` | UCCSD pool, point-group screening, spin-flip tying, Jordan–Wigner CNOT counts |
| `gva/engine.py` | Exact determinant-space simulator (Givens-rotation form of each exponential), adjoint gradients, L-BFGS-B VQE, ADAPT-VQE |
| `gva/qubit.py` | Jordan–Wigner Hamiltonian, CI-to-qubit mapping, QWC/GC grouping (largest-first, DSATUR, sorted insertion), shot estimates |
| `gva/noise.py` | Global depolarizing error budget, Richardson ZNE, rescaling overhead |
| `gva/periodic.py` | Crystals: k-point integrals (gauge-fixed Bloch orbitals, real supercell basis), crystal-momentum screening, spin-flip and time-reversal tying, complex determinant-space simulator, closed-form counts for any mesh, exact JW builder |
| `uccsd_counts.py` | Closed-form parameter/CNOT formulas checked against enumeration |
| `run_small.py` | H2, H4, LiH, BeH2 benchmark (3 geometries, 5 starts) |
| `run_cisplatin.py` | Cisplatin CAS(10e,10o) benchmark (resumable) |
| `run_measure.py` | Measurement grouping and shot counts |
| `run_scaling.py` | Parameter/CNOT counts vs active-space size |
| `run_pes.py` | Pt–Cl bond-stretching scan (resumable) |
| `run_orbitals.py` | Active natural-orbital irreps and composition |
| `run_periodic.py` | H2-chain benchmark (dimerized chain N_k = 2-5; uniform chains 1.0/1.5/2.0 A, N_k = 3, 5); resumable; optional subset argument for parallel processes |
| `merge_periodic.py` | Merges the parallel parts into `results/periodic.json` |
| `run_scaling_periodic.py` | Counts for 1D and 3D meshes; closed forms checked by enumeration on 12 meshes |
| `run_measure_periodic.py` | Pauli terms, grouping and shots in the Bloch and real bases (20-qubit chains) |
| `verify_periodic.py` | Rebuilds every periodic state from stored parameters; exact JW energies |
| `analyze_tying_periodic.py` | Untied optimum against the imposed time-reversal and spin-flip relations |
| `tests_periodic.py` | Integral checks, three FCI routes, gradient checks, zero gradients of momentum-violating excitations, JW energy |
| `check_jw_cisplatin.py` | Cisplatin JW energy with and without OpenFermion's internal 1e-8 pruning |
| `make_figures_periodic.py` | Figures 9, 10 and Supplementary Fig. 2 |
| `ms_periodic.py`, `si_periodic.py` | Crystal sections of the manuscript and SI |
| `make_figures.py`, `fig1_schematic.py` | Figures 1–8 and Supplementary Fig. 1 |
| `make_workbook.py` | Supplementary Data 1 (GVA_Source_Data.xlsx) |
| `ms_build.py`, `si_build.py` | Manuscript and SI documents; every quoted number is read from `results/` via `facts.py` |
| `log_references.py` | Re-runs every reference calculation with full PySCF output and checks it against the stored energies |
| `tests_engine.py`, `tests_qubit.py` | Engine validation (gradients vs finite differences, convergence to FCI, Jordan–Wigner energy check) |
| `results/` | JSON outputs, cisplatin active-space integrals (`cisplatin_cas10_10.npz`) |
| `results/logs/` | v2 run logs, full PySCF output (`pyscf/`), `verify_theta.log`, `compare_v1_v2.txt` |
| `results_v1/` | Original (v1) results and logs, kept unchanged; v2 reproduces them exactly |
| `run_orbital_sensitivity.py` | Cisplatin with MP2 natural, canonical RHF and CASSCF orbitals (Supplementary Note 10) |
| `verify_theta.py` | Rebuilds every final state from the stored parameters and checks energy and ⟨S²⟩ |
| `compare_v1_v2.py` | Run-by-run comparison of v2 with v1 |
| `METHODS_AND_REFERENCES.md` | For each function: what it computes, equation, settings, literature (43 verified references) |
| `PROVENANCE.md` | Chronology of all runs, code changes during the work, aborted runs, what was not kept |
| `environment.txt` | Python and package versions, `pip freeze`, OS, CPU, memory |
| `CHECKSUMS.sha256` | SHA-256 of every file (`sha256sum -c CHECKSUMS.sha256`) |
| `DIFFERENCES_FROM_ORIGINAL_ko.md` | Internal note (Korean) on why the values differ from the earlier draft; remove before public deposition |

## Figure numbering (npj Quantum Information version)

Scripts and `figures/` keep historical file names; `numbering.py` maps them to the
submitted numbering, and `make_release.py` writes renamed copies to `figures_npj/`.

| figures/ | submitted | content |
|---|---|---|
| Fig1 | Fig. 1 | construction of GVA |
| Fig2 | Fig. 2 | small molecules |
| Fig4 | Fig. 3 | cisplatin CAS(10e,10o) |
| Fig9 | Fig. 4 | savings against symmetry-group order |
| Fig10 | Fig. 5 | periodic hydrogen chains |
| FigS1 | Supplementary Fig. 1 | spin contamination, small molecules |
| Fig6 | Supplementary Fig. 2 | resource scaling, cisplatin active spaces |
| Fig5 | Supplementary Fig. 3 | cisplatin accuracy against circuit resources |
| Fig3 | Supplementary Fig. 4 | measurement grouping and shots |
| Fig7 | Supplementary Fig. 5 | depolarizing error budget |
| Fig8 | Supplementary Fig. 6 | Pt-Cl bond stretching |
| FigS2 | Supplementary Fig. 7 | tying ablation, chains |

## Order of execution

```
# redirect stdout to results/logs/<script>.log to keep a run log
OMP_NUM_THREADS=1 python run_small.py          # ~50 min
python run_cisplatin.py      # ~3 h
python run_measure.py        # ~25 min
python run_scaling.py
python run_pes.py            # ~2 h
python run_orbitals.py
python run_orbital_sensitivity.py   # ~15 min
python verify_theta.py       # rebuilds all final states from stored theta (~5 min)
python run_scaling_periodic.py
# periodic benchmark: about 11 h of CPU time in total (each 20-qubit system 1-2 h);
# resumable, any subset of system indices, parts can run in parallel
OMP_NUM_THREADS=1 python run_periodic.py 0,1,2,3 C &
OMP_NUM_THREADS=1 python run_periodic.py 4,6,8,5 F &
OMP_NUM_THREADS=1 python run_periodic.py 7,9 G
python merge_periodic.py
python run_measure_periodic.py
python verify_periodic.py && python analyze_tying_periodic.py
python tests_periodic.py && python check_jw_cisplatin.py
python uccsd_counts.py       # exits non-zero if any formula disagrees with enumeration
python log_references.py     # optional: full PySCF reference logs + energy check (v1 output kept in results_v1/logs/reference)
python fig1_schematic.py && python make_figures.py && python make_figures_periodic.py
python make_workbook.py
python ms_build.py && python si_build.py
```

Timings refer to a 2-core cloud VM. No quantum hardware is used anywhere.
`run_cisplatin.py` and `run_pes.py` resume from existing JSON files; delete
`results/cisplatin_eq.json` and `results/pes.json` to recompute from scratch.

## Closed-form counts (n_o occupied, n_v virtual spatial orbitals)

```
N_UCCSD = 2 n_o n_v + 2 C(n_o,2) C(n_v,2) + (n_o n_v)^2
N_spin  = n_o n_v + C(n_o,2) C(n_v,2) + n_o n_v (n_o n_v + 1)/2
N_PG    = 2S + 2D_ss + D_ab
N_GVA   = S + D_ss + (D_ab + n_o n_v)/2
  S    = sum_G o_G v_G
  D_ab = sum_G X(G)^2,  X(G) = sum_{G1 x G2 = G} o_G1 v_G2
  D_ss = sum_G P_o(G) P_v(G)
CNOT (ladder, interleaved JW): single p<q: 4(q-p); double p<q<r<s: 16(q-p+s-r+1)

Crystals (any abelian group; inverse of k is -k):
  D_ab = sum_q X(q) X(-q),  X(q) = #(i,a) with k_a - k_i = q
  D_ss = sum_K P_o(K) P_v(K)
  n_self = #(i,a) with 2(k_a - k_i) = 0
  N_GVA-k = S + D_ss + (D_ab + n_self)/2      (real parameters)
  o/v bands per cell: N_GVA-k ~ (3/4) N_k^3 o^2 v^2 ; N_UCCSD(supercell) ~ 1.5 N_k^4 o^2 v^2
```
