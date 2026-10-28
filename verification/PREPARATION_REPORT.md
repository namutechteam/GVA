# Local preparation validation

Date: 2026-10-28 (Asia/Seoul).

## Passed

- 273 scientific package files match the supplied reproduction ZIP byte for byte.
  The internal Korean comparison note is intentionally omitted per the source README.
- 58 Python files compiled for syntax (without running heavy calculations).
- 29 JSON result files parsed.
- 12 NPZ files opened and every array loaded with allow_pickle=False.
- Source workbook opened: 31 worksheets, including stored molecular/periodic runs,
  histories, variational parameters, figure values, and enumeration checks.
- All 12 final plotting scripts completed successfully and generated 36 nonempty
  PDF/SVG/PNG exports in output/figures/.
- Cisplatin stored orbital irreps: formula counts agree with explicit enumeration:
  UCCSD 875; spin-only 450; point-group-only 297; GVA 161 parameters.
  Enumerated CNOT counts: UCCSD 108800; GVA 36496 per layer.
- UCCSD parameter formula agrees with full excitation-pool enumeration for all
  36 combinations of 1–6 occupied and 1–6 virtual spatial orbitals.

## Limits and remaining checks

- PySCF and OpenFermion are absent from the current host Python environment.
  Molecular/periodic energy reconstruction, gradients and full optimization were
  not rerun. Existing verification logs are supplied provenance, not new validation.
- Workbook and stored results were inventoried and read; this is not a cell-by-cell
  audit of all manuscript tables or numerical claims.
- Figure executions check functionality; no new visual review of all panels was
  performed. Original exports are preserved. DejaVu Sans fallback can alter text layout.
- Nonfatal font fallback and existing invalid-escape warnings appeared during
  rendering. All scripts exited successfully.
- Original requirements versions are supplied metadata; a clean installation of
  the complete scientific environment has not been tested here.
- Code licence was subsequently set to MIT by the project owner (2026 NamuICT).
  Data licensing, author metadata, DOI and GitHub publication remain pending.
  No external upload took place.

## Preparation environment

Python: 3.12.12

Packages: {"matplotlib": "3.10.8", "numpy": "2.2.6", "openpyxl": "3.1.5", "scipy": "1.18.0"}
