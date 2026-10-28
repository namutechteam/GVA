# Preparation provenance

Prepared on 2026-10-28 (Asia/Seoul).

## Inputs

- `GVA_final/GVA_reproduction_package.zip`: final supplied scientific package.
- `GVA_final/GVA_Source_Data.xlsx`: final supplied source workbook.
- `GVA_final/GVA_calculation_logs.zip`: supplied logs.
- `GVA_final_work/figures/*.py` and `out/`: final redraws referenced by the current
  manuscript and SI in `paper/`.

## Preparation changes

- Scientific code and results are copied without edits, except that the internal
  `DIFFERENCES_FROM_ORIGINAL_ko.md` is omitted as explicitly requested by the
  supplied package README. It remains available in the source archive.
- Original `CHECKSUMS.sha256` is retained as historical source metadata. It is not
  the root manifest: its omitted internal note no longer exists in this copy.
- Final plotting code is copied; shared font configuration is changed from Arimo
  to Matplotlib-bundled DejaVu Sans. Supplied figure exports remain unchanged.
- Font binaries without accompanying redistribution notices are not copied.
- Portable figure and integrity commands, requirements, and English README added.
- No research result is recalculated or replaced during packaging.
- No GitHub repository, public release, Zenodo deposition, or licence is created.

See `verification/source_manifest.json` for hashes of the inputs and
`verification/files.sha256` for hashes of the prepared files.

## README revision

Reviewed MGQC-BT, TriSPD and Keap1-Nrf2-peptide-design README layouts.
The root README now leads with the study, workflow, contents and numbered
execution steps. Detailed preparation notes are in REPRODUCTION_NOTES.md.

## MIT licence selection

On 2026-10-28 the project owner selected MIT for the code. Added LICENSE with
Copyright (c) 2026 NamuICT and updated the README. Data, logs and figure exports
remain outside this code-licence declaration until their terms are specified.
The earlier preparation status above records the state before this selection.

## Data and figure licence selection

On 2026-10-28 the project owner selected CC BY 4.0 for original data, outputs,
logs and project-created figures. LICENSE-DATA.md states the scope and links
to the official legal code. Code remains MIT; third-party rights are preserved.
The initial private GitHub upload was confirmed by the owner. These licence changes are committed from the server while retaining the
initial GitHub history. The repository remains private.

## Public repository

The project owner authorized public visibility on 2026-10-28. GitHub confirms
that namutechteam/GVA is public. A pre-publication scan of tracked files and
log archive entries found no known GitHub token, AWS key ID or private-key
patterns; this was a targeted check, not an exhaustive security audit.
The working manuscript now cites source-data commit
f65a20ed047628c767385dab5690079d86c2977a, MIT code and CC BY 4.0 original data.
No DOI or tagged release is claimed.
