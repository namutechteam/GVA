# Provenance of the results in this package

This file records how every number in the manuscript and Supplementary
Information was produced: which script and code version, when, which runs
were aborted or discarded, what was changed in the code during the work, and
what was not kept. Times are Korea Standard Time (UTC+09:00) on 2026-10-23
unless a date is given.

**Current state.** The results in `results/` come from a complete rerun (v2,
Section 8) that stores the final parameters of every run and full PySCF
output. All v2 energies and iteration counts are identical to the original
run (v1, Sections 1–7), which is kept unchanged in `results_v1/`.
Times marked `≈` are reconstructed from run durations; all others are file
modification times or log entries.

## 1. Chronology (v1; paths below are relative to `results_v1/`)

| Time | Step | Output | Used in paper |
|---|---|---|---|
| ≈17:20–17:30 | Engine validation (`tests_engine.py`, `tests_qubit.py`); exploratory cisplatin set-up and timing scripts (`tests_cis.py`, `tests_cis2.py`, not packaged) | console only | No (validation only) |
| ≈17:33 | `run_small.py`, attempt 1 | discarded | No |
| ≈17:35 | `run_cisplatin.py`, attempt 1 | discarded | No |
| ≈17:46 | Both attempts stopped; **code change 1** (Section 2) | — | — |
| 17:46:39 → 18:37:23 | `run_small.py`, attempt 2 (final) | `results/small_molecules.json`, `logs/run_small.log` | Table 1, Fig. 2, Fig. 3 (small molecules), Fig. 7, Tables S4, S8, Fig. S1 |
| ≈17:46:30 → 17:59:57 | `run_cisplatin.py`, attempt 2, part 1: GVA 1L seeds 0–2, GVA 2L seed 0 | records kept in `results/cisplatin_eq.json`, `logs/run_cisplatin_part1.log` | Table 2, Figs 4, 5, 7, Table S3 |
| ≈18:00 | Part 1 stopped; **code change 2** (Section 2) | — | — |
| 18:01:10 → 20:11:20 | `run_cisplatin.py`, part 2 (resumed; skips the four finished runs) | `results/cisplatin_eq.json`, `results/cisplatin_cas10_10.npz`, `logs/run_cisplatin_part2.log` | Table 2, Figs 1, 4, 5, 7, Tables S3, S5 |
| 18:37:23 → 19:01:07 | `run_measure.py` | `results/measurement.json`, `logs/run_measure.log` | Fig. 3, Table S6, Fig. 7 (identity coefficients) |
| 19:01:07 → 19:02:30 | `run_scaling.py` | `results/scaling.json`, `logs/run_scaling.log` | Fig. 6, Table S5 |
| 19:02:30 → 19:16:08 | `run_pes.py`, version 1 (12-point grid); point R = 2.32 Å completed | record kept in `results/pes.json`, `logs/run_pes_part1.log` | Fig. 8, Table S7 |
| 19:08:46 | `run_orbitals.py` (run alongside) | `results/cisplatin_orbitals.json` | Fig. 4a, Table S9 |
| ≈19:16 | Version 1 stopped (≈14 min per point); **code change 3** (Section 2) | — | — |
| ≈19:17 → 20:28:15 | `run_pes.py`, version 2 (9-point grid, resumed; reuses R = 2.32 Å) | `results/pes.json`, `logs/run_pes_part2.log` | Fig. 8, Table S7 |
| 20:28–20:34 | `fig1_schematic.py`, `make_figures.py`, `make_workbook.py`, `ms_build.py`, `si_build.py` | `figures/`, `results/noise_budget.json`, workbook, documents | all |
| 21:29–21:40 | This supplement: `environment.txt`, validation logs, `log_references.py` | `results/logs/validation/`, `results/logs/reference/` | verification |

## 2. Code changes made during the work

All final results were produced with the code in this package. The three
changes below were made between attempts; none alters the electronic-structure
model or the ansatz definitions.

**Change 1 (`gva/engine.py`, `run_small.py`, `run_cisplatin.py`; ≈17:46).**
Attempt 1 used an L-BFGS-B relative energy tolerance `ftol = 1e-15` and an
ADAPT-VQE loop without a stall test (limit 200 operators). In attempt 1 the
H4 ×1.0 and ×1.5 geometries took 68 s and 91 s, against 7 s and 5 s after the
change; a stack sample showed the time going into ADAPT-VQE re-optimizations.
The cisplatin two-layer GVA run was still converging after more than 10 min.
Changes:
`run_vqe` gained the argument `ftol` with default `1e-12` (≈1e-9 Ha absolute
for cisplatin, ≈1e-11 Ha for the small molecules); `adapt_vqe` gained a stall
test (energy gain below 1e-7 Ha over five additions), a `stop_reason` field,
and re-optimization settings `maxiter = 1500`, `ftol = 1e-12`; the operator
limit in the run scripts became 120. Both attempts were restarted from scratch;
nothing from attempt 1 is used (its stdout was overwritten and no JSON was written).

**Change 2 (`run_cisplatin.py`; ≈18:00).** Added resumption from
`results/cisplatin_eq.json` (finished runs are skipped) and reduced the
two-layer UCCSD random starts from seeds (0, 1, 2) to (0, 1) to fit the time
budget. The four runs finished before the change were computed by the same
`engine.py`, `ansatz.py` and `systems.py` as all later runs.

**Change 3 (`run_pes.py`; ≈19:16).** Version 1 computed each point fully
before moving on. Version 2 first builds all active-space Hamiltonians with
orbital tracking (phase 1), then runs references and VQE with resumption
(phase 2), on a 9-point grid (2.0, 2.2, 2.32, 2.5, 2.7, 2.9, 3.2, 3.5, 3.8 Å)
instead of 12 points. The R = 2.32 Å point from version 1 was kept; version 2
asserts that its CASCI energy agrees with the newly built Hamiltonian to 1e-8 Ha.
The per-point computation (functions, `maxiter = 2500`, Hartree–Fock start)
is identical in both versions.

**Other edits that do not affect stored results.**
`gva/systems.py`: the `symmetry` argument of `cisplatin_cas_hamiltonian` was
added (default unchanged) while part 2 of the cisplatin run was in progress
with the module already loaded; it is used only by `run_pes.py`. `gva/noise.py`: the rescaling
overhead was rewritten in logarithmic form to avoid floating-point overflow
(analysis stage, before `noise_budget.json` was written). Figure, workbook and
document scripts were edited repeatedly for layout and wording.

## 3. Runs that are not part of the results

| Run | Reason | Data kept |
|---|---|---|
| `run_small.py` attempt 1 | stopped for code change 1 | none |
| `run_cisplatin.py` attempt 1 | stopped for code change 1 | none |
| `run_pes.py` version 1 beyond R = 2.32 Å | stopped for code change 3 | none (no other point had finished) |
| Exploratory scripts `tests_cis.py`, `tests_cis2.py` | set-up checks and timing | none |

No completed run was dropped from any table or figure. Every run listed in
`results/small_molecules.json` and `results/cisplatin_eq.json` appears in
Supplementary Data 1.

## 4. What was not kept in v1 (resolved in v2, Section 8)

* Final optimized parameter vectors θ of the VQE and ADAPT-VQE runs, and the
  list of operators chosen by ADAPT-VQE (only its energy and gradient trace).
  The final states can be regenerated only by repeating the optimizations,
  which are deterministic for a given seed and thread count.
* Small-molecule Hamiltonians (regenerated deterministically by
  `gva.systems.small_molecule_hamiltonian`).
* PySCF output of the production runs (`verbose = 0`). The reference steps
  (RHF, MP2, active-space construction, CISD, CCSD, CCSD(T), FCI/CASCI) were
  re-run with full output by `log_references.py`; the logs are in
  `results/logs/reference/` and every energy agrees with the stored value to
  better than 1e-7 Ha (`results/logs/reference/SUMMARY.txt`).

## 5. Determinism

All production runs used `OMP_NUM_THREADS=1`. VQE starts are fixed by
`numpy.random.default_rng(1000 + k)` (small molecules) or `(2000 + k)`
(cisplatin). Re-running on other hardware or BLAS builds can change energies at
the 1e-10 Ha level and iteration counts by a few percent.

## 6. Logs

| File | Content |
|---|---|
| `results/logs/run_small.log`, `run_cisplatin.log`, `run_pes.log`, `run_orbital_sensitivity.log` | v2 run logs (one line per run or point; ADAPT-VQE energy and gradient per added operator) |
| `results/logs/pyscf/` | v2: full PySCF output (SCF, MP2, CISD, CCSD, CCSD(T)) for every molecule and active space |
| `results/logs/verify_theta.log`, `compare_v1_v2.txt` | v2: reconstruction from stored parameters; run-by-run comparison with v1 |
| `results_v1/logs/run_small.log`, `run_cisplatin_part1.log`, `_part2.log`, `run_measure.log`, `run_scaling.log`, `run_pes_part1.log`, `_part2.log` | v1 run logs |
| `results_v1/logs/validation/` | `tests_engine.py`, `tests_qubit.py`, `uccsd_counts.py` output |
| `results_v1/logs/reference/` | full PySCF output of all v1 reference calculations re-run by `log_references.py`, and `SUMMARY.txt` |
| `environment.txt` | Python, package versions, `pip freeze`, OS, CPU, memory |

Per-iteration energies of every VQE run are stored in the `history` field of
`results/small_molecules.json` and `results/cisplatin_eq.json`.

## 7. Checksums

`CHECKSUMS.sha256` lists SHA-256 hashes of every file in the package.
Verify with `sha256sum -c CHECKSUMS.sha256`.

## 8. Rerun v2 with parameter storage (23–24 October 2026)

**Code changes before v2 (no change to the calculation itself).**
`gva/systems.py`: optional PySCF logging (`set_pyscf_log`, `_mole`).
`gva/ansatz.py`: `excitation_label` and `parameter_map`.
`run_small.py`, `run_cisplatin.py`, `run_pes.py`: store for every run the final
parameters `theta`, start and finish times, the L-BFGS-B message, for ADAPT-VQE
the chosen operators and their parameters, and per system the generator list
behind each parameter (`parameter_maps_1L`); write PySCF output to
`results/logs/pyscf/`. Run settings, seeds and orbital choices are unchanged.

| Time | Step | Output |
|---|---|---|
| 2026-10-23 21:47:48 | v1 results copied to `results_v1/`; v2 runs started on two cores | — |
| 21:47:48 → 22:38:36 | `run_small.py` | `results/small_molecules.json`, `logs/run_small.log` |
| 21:47:48 → 2026-10-24 00:14:46 | `run_cisplatin.py` (8,818 s, from scratch) | `results/cisplatin_eq.json`, `results/cisplatin_cas10_10.npz`, `logs/run_cisplatin.log` |
| 22:38:36 → 00:08:18 | `run_pes.py` (9 points, from scratch) | `results/pes.json`, `logs/run_pes.log` |
| ≈22:40 → 22:54:56 | `run_orbital_sensitivity.py` (new; three orbital choices for cisplatin) | `results/orbital_sensitivity.json`, `logs/run_orbital_sensitivity.log` |
| 00:15 | `compare_v1_v2.py` | `logs/compare_v1_v2.txt` |
| 00:19 | `verify_theta.py` | `results/verify_theta.json`, `logs/verify_theta.log` |
| 00:20 onward | figures, workbook, manuscript and SI rebuilt from v2 | — |

**Reproduction of v1.** All 276 small-molecule runs, 16 cisplatin runs and
36 scan values (27 VQE runs and 9 CASCI energies) agree with v1 to 0 Ha, with
identical iteration counts (`logs/compare_v1_v2.txt`). `results/noise_budget.json`
is byte-identical to the v1 file.

**Reconstruction from stored parameters.** `verify_theta.py` rebuilt all
319 final states (276 + 16 + 27) from the stored parameters and ansatz
definitions; energy and ⟨Ŝ²⟩ agree with the stored values with zero
difference (tolerances 1e-9 Ha and 1e-8).

**Data that are still not stored.** Full final state vectors (they follow from
θ in less than a second each) and intermediate parameter vectors during
optimization (only the energy after each iteration is kept).

**Orbital sensitivity (new in v2).** Cisplatin CAS(10e,10o) was also run with
canonical RHF orbitals 32–41 (HOMO−4 to LUMO+4) and with CASSCF(10e,10o)
orbitals (Supplementary Note 10, Table S10). The MP2-natural-orbital entry
reproduces the main results (CASCI energy asserted equal to 1e-8 Ha).

**Text change, 2026-10-24.** A Discussion paragraph on molecules without spatial
symmetry was added (`ms_build.py`), with five new references (`refs.py`,
`references_verified.txt`), and the caption of Supplementary Table S5 now states
that tying alone leaves the CNOT count at the UCCSD value. All numbers in the new
paragraph are read from existing results (`cisplatin_eq.json` counts "spin only",
`scaling.json`, the cisplatin ADAPT-VQE and GVA runs); no calculation was rerun.

## 9. v3: extension to crystals (2026-10-24)

**Motivation.** For molecules the order of an abelian point group is at most eight,
so the saving from screening is bounded and vanishes for C1 molecules. The
translation group of a crystal has order N_k and grows with the k-point mesh.
v3 adds crystal-momentum screening with spin-flip and time-reversal tying (GVA-k),
closed-form counts for any abelian group, and a VQE benchmark on periodic H2 chains.

**New code.** `gva/periodic.py` (integrals, gauge fixing, orbit construction,
complex simulator, counts, exact JW builder); run, test, verification, figure and
document scripts listed in README. The molecular code (`gva/engine.py`,
`gva/ansatz.py`, `gva/systems.py`, `gva/qubit.py`) is unchanged; the molecular
results of v2 were not recomputed.

**Development runs not kept.** Before the production runs, the periodic code was
exercised on N_k = 2-5 test systems (gradient checks, counts, first VQE runs). A
first production launch at 02:47 was stopped after about one minute (it had been
started as detached processes, which the environment may stop; it was restarted
as tracked jobs) and its partial logs were deleted. A one-system test run of
`run_periodic.py` (N_k = 2) was also deleted before the production runs.

**Findings during development.**
* For N_k = 2, every ansatz from the zero start stops 0.29 uHa above FCI; full
  BFGS from a random start reaches FCI to 1e-14 Ha. Left as is (Supplementary Note 14).
* Two-layer runs from the zero start can stall at a symmetric stationary point
  (seen once in testing at r = d = 1.0 A, N_k = 3); three starts are therefore
  used for all two-layer settings.
* OpenFermion removes coefficients below 1e-8 while accumulating operators. With
  this default, Jordan-Wigner energies differ from determinant-space energies by
  ~1e-9 Ha for the uniform chains (many integrals below 1e-8). `qubit_hamiltonian_c`
  disables the pruning; agreement is then ~1e-15 Ha. For cisplatin the difference is
  2.5e-12 Ha with or without pruning (`check_jw_cisplatin.py`); the "1e-8 Ha" quoted
  in the v2 manuscript was the tolerance asserted by run_measure.py, not a measured
  difference, and the v3 text gives the measured value. The run-time JW values stored in `results/periodic.json` were
  computed with the earlier builder (pruned for the Bloch basis); the exact values
  are in `results/verify_periodic.json`.
* Run-to-run differences in iteration counts were seen in tests with
  multithreaded BLAS; production runs use OMP_NUM_THREADS=1.

**Production runs and the cost reduction.** Parts A (systems 3,5,0,4) and B
(7,9,1,2,6,8) were started at 02:47 as two tracked jobs and stopped by the
environment's 2-hour limit at 04:47. By then the two 20-qubit systems processed
first (dimerized N_k = 5; uniform 1.5 A N_k = 5) were complete except for one
GVA-k 2L start and ADAPT-VQE; the interrupted run left no record. Two-layer
supercell UCCSD from random starts had needed 1,453-1,972 s at the 3,000-iteration
cap. The remaining work was relaunched as parts C (0,1,2,3), D (7,4,6,8,5) and E (9),
which resume from all existing part files. For the two remaining 20-qubit uniform
chains (1.0 and 2.0 A) two-layer supercell UCCSD was run from the Hartree-Fock start
only; every other setting is unchanged.

**ADAPT-VQE wall cap.** Part D (04:47-06:47) spent its whole 2-h limit inside
ADAPT-VQE for the uniform 1.5 A chain (N_k = 5) and was stopped before writing that
run; part E completed every run for the 2.0 A chain except ADAPT-VQE, which was
stopped by hand at 06:49 for the same reason. `adapt_vqe_c` received a `max_wall`
argument (stop before adding a further operator once 1,800 s have elapsed; the state
reached is recorded with stop_reason 'time limit'). The remaining work was run as
parts F (4,6,8,5) and G (ADAPT-VQE for 7 and 9). ADAPT-VQE runs from parts C and the
smaller systems finished well inside the cap.

**Completion of the periodic runs.** Part C finished at 05:17, part G at 07:52 and
part F at 08:06 (all within their limits). `merge_periodic.py` combined the six part
files into `results/periodic.json` (10 systems, 176 runs, no duplicates).

**Verification and measurement.** A first `verify_periodic.py` run, which built the
exact (unpruned) Jordan-Wigner operator for every system, was stopped after reaching
the 20-qubit systems (each 20-qubit build takes more than 10 min). The exact check is
now restricted to systems with up to 16 qubits; all 176 states are still rebuilt
from the stored parameters (zero difference in energy and <S^2>), and the 20-qubit
states keep the JW check made during the runs (max 5.6e-13 Ha). A first
`run_measure_periodic.py` run with the exact builder was stopped after 21 min for
the same reason and restarted with OpenFermion's default pruning during
construction (the counts prune at 1e-8 in either case, as for the molecular counts).

## 10. npj Quantum Information version (2026-10-24)

Text and layout only; no calculation was rerun.
* Abstract cut to 150 words; Introduction heading added; Code availability moved into
  Methods and Data availability placed after Methods (npj order).
* Five molecular figures moved to the SI (measurement grouping, cisplatin accuracy
  against resources, active-space scaling, depolarizing error budget, Pt-Cl scan);
  their paragraphs are now at the start of Supplementary Notes 4, 6, 7, 8 and 9, and a
  one-paragraph summary remains in the main text. Main display items: 5 figures, 3 tables.
* Figure numbers are mapped by `numbering.py`; workbook sheet names follow the map;
  `figures_npj/` holds renamed figure files.
* The SI has its own reference list (Supplementary References); the main reference
  list contains only references cited in the main text (50).
* Wording changes from a final style pass (no change of content).
* Every cited reference was re-checked against Crossref/OpenAlex/arXiv metadata
  (results in references_verified.txt).
