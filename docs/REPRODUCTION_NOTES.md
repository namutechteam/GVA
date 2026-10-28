# Reproduction notes

## Font handling

The final redraw scripts originally used Arimo. The supplied font binaries
had no accompanying redistribution notice, so they are not bundled here.
The portable plotting configuration uses Matplotlib-bundled DejaVu Sans;
`figures/` preserves the original exports. Regenerated labels can differ in
spacing and line wrapping. Add the appropriate redistribution notice before
bundling Arimo for exact typography.

## Scientific package and final redraws

`code/gva_repro/` is copied from the final supplied reproduction package.
Its original `make_figures*.py` and `fig1_schematic.py` build historical figure
designs; `numbering.py` maps those names to the manuscript numbering.
`plotting/` holds the final redraw scripts for main Figures 1–5 and
Supplementary Figures 1–7. Use these for the current designs.
The original `ms_build.py` and `si_build.py` build historical documents and
do not incorporate all subsequent figure and layout edits.

## Dataset provenance

Scientific inputs come from `GVA_final`, supplied after the earlier local v2.4
work. The historical `results_v1/` directory within the final package is retained
for the package's own provenance and is not the separate local v2.4 release.
The source archives and their hashes are listed in
`../verification/source_manifest.json`.

## Release status

This folder is a local preparation for a repository named `namutechteam/GVA`.
The owner confirmed the initial private GitHub upload. No public release or
DOI has been created as part of preparation.
Author spellings, final author-order approval, corresponding emails
and archival metadata remain in `RELEASE_CHECKLIST.md`.
The manuscript and SI in the working `paper/` directory still contain metadata
placeholders and are not included in this preparation.

## README references

The README layout was reviewed against the company repositories
`namutechteam/MGQC-BT`, `namutechteam/TriSPD` and
`namutechteam/Keap1-Nrf2-peptide-design` on 2026-10-28.
The GVA commands and reported checks refer to this package only.

The code and software documentation now use MIT (2026 NamuICT), as selected
by the project owner. Original datasets, logs and project-created figures are separately licensed
under CC BY 4.0; see LICENSE-DATA.md.
