# GVA: Symmetry-screened unitary coupled-cluster circuits for molecules and crystals

Code and data for the manuscript *"Symmetry-screened unitary coupled-cluster
circuits for molecules and crystals"*.

The geometric variational ansatz (GVA) uses symmetry to reduce the number of
independent amplitudes and two-qubit gates in UCCSD circuits. Molecular
calculations combine point-group screening with spin-flip tying; the periodic
extension, GVA-k, uses crystal-momentum screening and time-reversal tying.
The benchmarks cover H₂, H₄, LiH, BeH₂, cisplatin CAS(10e,10o), and periodic
hydrogen chains. All calculations use classical simulation.

## Workflow

```text
01  Reference problem   geometry, orbitals, integrals and FCI/CASCI reference
02  Ansatz             UCCSD excitation pool, symmetry screening and tying
03  Optimization       determinant-space simulation and adjoint gradients
04  Resource analysis  parameter/CNOT counts, measurement groups and shot estimates
05  Figures            stored source data → manuscript and supplementary figures
```

## Contents

| Path | Content |
| --- | --- |
| [`code/gva_repro/gva/`](code/gva_repro/gva/) | Hamiltonians, ansatz construction, simulation and resource analysis |
| [`code/gva_repro/results/`](code/gva_repro/results/) | Stored runs, energy histories, parameters and active-space integrals |
| [`data/GVA_Source_Data.xlsx`](data/GVA_Source_Data.xlsx) | Supplementary Data 1: source values for figures and tables |
| [`plotting/`](plotting/) | Scripts for main Figures 1–5 and Supplementary Figures 1–7 |
| [`figures/`](figures/) | Manuscript figure exports in PDF, SVG and PNG |
| [`archives/GVA_calculation_logs.zip`](archives/GVA_calculation_logs.zip) | Calculation and reference logs |
| [`verification/`](verification/) | File hashes and local preparation checks |

## 1. Set up the environment

Use Python 3.11 or newer. To regenerate figures from the included workbook:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-figures.txt
```

Scientific calculations additionally require PySCF, OpenFermion and SciPy:

```bash
python -m pip install -r requirements.txt
```

The recorded calculation environment is in
[`code/gva_repro/environment.txt`](code/gva_repro/environment.txt).
The requirements files are installation inputs; they are not a validated
lockfile for a clean installation.

## 2. Reproduce the figures

From the repository root:

```bash
python scripts/reproduce_figures.py
```

This reads the included workbook and writes 12 figures, each as PDF, SVG and
PNG, to `output/figures/`. No new quantum-chemistry calculation or GPU is needed.
Use `--output /path/to/new-output` to choose another output directory.

| Scripts | Manuscript items |
| --- | --- |
| `plotting/make_fig1.py` | Fig. 1: construction of GVA |
| `plotting/make_fig2.py` | Fig. 2: small-molecule benchmarks |
| `plotting/make_fig3.py` | Fig. 3: cisplatin |
| `plotting/make_fig4.py` | Fig. 4: savings versus symmetry-group order |
| `plotting/make_fig5.py` | Fig. 5: periodic hydrogen chains |
| `plotting/make_suppfig1.py`–`make_suppfig7.py` | Supplementary Figs. 1–7 |

Figure 1 includes schematic text and selected illustrative values in the script;
the other scripts load workbook sheets. Plotting stored values does not rerun
the underlying optimizations.

The regeneration scripts use Matplotlib's DejaVu Sans. The included manuscript
exports use Arimo, so text spacing and line wrapping may differ. See
[`docs/REPRODUCTION_NOTES.md`](docs/REPRODUCTION_NOTES.md) for font details and
historical figure numbering.

## 3. Verify stored states and run calculations

After installing the scientific dependencies:

```bash
cd code/gva_repro
python tests_engine.py
python tests_qubit.py
python tests_periodic.py
python verify_theta.py
python verify_periodic.py
```

The `verify_*` scripts rebuild final states from stored parameters and compare
energies and spin expectations with the recorded results. The test scripts
check the simulation engine, qubit mapping and periodic calculations; see their
source for the checks performed.

To rerun the benchmarks, follow the execution order in
[`code/gva_repro/README.md`](code/gva_repro/README.md). It lists the molecular,
cisplatin and periodic runners, resource estimates and resume behaviour.
Use a separate copy for new calculations to preserve the included reference
results. Full calculations can take hours.

## Data and reproducibility

`data/GVA_Source_Data.xlsx` contains the manuscript source tables, optimization
histories and variational parameters. `code/gva_repro/results/` also includes
cisplatin and periodic-system integrals. Full reference and execution logs are
included in the code package and the separate log archive.

- [`PROVENANCE.md`](code/gva_repro/PROVENANCE.md): calculation history and code changes.
- [`METHODS_AND_REFERENCES.md`](code/gva_repro/METHODS_AND_REFERENCES.md): implementation details and literature.
- [`PREPARATION_REPORT.md`](verification/PREPARATION_REPORT.md): checks executed during repository preparation and their limits.
- [`PAPER_DATA_AUDIT.md`](verification/PAPER_DATA_AUDIT.md): manuscript/SI consistency checks against the stored data.

Check the integrity of the included files from the repository root:

```bash
python scripts/verify_manifest.py
```

This checks SHA-256 hashes; it does not validate scientific claims.

## Citation

Citation metadata will be finalized after the remaining author details are
confirmed. The repository is public; no archival DOI has been assigned.

## License

The research code and associated software documentation are released under the
[MIT License](LICENSE). Copyright © 2026 NamuICT.

Original research datasets, numerical outputs, calculation logs and project-created
figure exports are licensed under [CC BY 4.0](LICENSE-DATA.md), with scope and
attribution described in that notice. Third-party materials retain their original
terms.
