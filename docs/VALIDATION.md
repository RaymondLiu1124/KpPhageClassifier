# KpPhageClassifier 0.1.2 release validation

[Back to overview](../README.md)

Local checks were performed on **18 November 2026**, using **Linux and Python 3.13**.

## Software and installation

**28 unit and regression tests passed.** Coverage includes per-genome reciprocal
matching, family ambiguity, order orientation/start invariance, rearrangements,
evidence thresholds, paralog/MTase exclusion, input IDs and coordinates,
existing-output preservation, installer integrity and transfer recovery,
portable paths, packaged release metadata and optional gene-caller checks.

A fresh virtual environment installed the package and DNA prediction dependencies.
The packaged CLI ran outside the source checkout. The wheel includes both the
reference-release registry and prediction helper. Source and wheel metadata passed
build/Twine checks; Ruff's `F` and `E9` checks passed. Dependency checks found no conflicts.

The offline installer verified the real pinned release archive and installed it
successfully into a path containing spaces. A previously verified live Figshare
download was used for the full classification run below.

## Gene calling and full classification

The packaged Pyrodigal-gv helper produced identical normalized protein sequences
and gene coordinates to the original Pharokka-based helper for three references.
Raw prediction files differ only in internal sequence counters, which do not
change protein sequences, coordinates or downstream identifiers.

The full workflow was **DNA → gene prediction → bidirectional MMseqs2 search →
assignment**, using Pyrodigal-gv **0.3.2**, Pyrodigal **3.7.1** and MMseqs2
**18.8cc5c**. References `G0053`, `G0524` and `G0119` were held out of scoring,
family assignment and order comparisons.

| Query | W group | Subgroup | wGRR |
| --- | --- | --- | ---: |
| `OQ830459.1` | W005 | W005-S02 | 0.890482 |
| `MK422451.1` | W020 | W020-S01 | 0.826011 |
| `LN866626.1` | W001 | W001-S05 | 0.935410 |
| `N_ONLY_CONTROL` | Unassigned | Unassigned | 0.000000 |

All computed assignment fields match the previous validated run. The previous run
explicitly labelled input completeness as `complete`; this run used the default
`unknown`. This user-supplied field differs; the inferred results do not.

Per-reference similarities, family assignments, gene-order comparisons,
normalized proteins and coordinates are unchanged. The N-only control is
unassigned with no predicted proteins. Reference core-file hashes were checked
and remain unchanged.

## Scope

These are engineering checks with known reference genomes. Original family
construction retained the held-out sequences; these checks are not independent
biological validation or calibrated accuracy estimates.

Linux and Python 3.13 were tested locally. After upload, the
[GitHub test and packaging matrix](https://github.com/lxd98/KpPhageClassifier/actions/runs/44263954708)
passed on **Python 3.10, 3.11, 3.12 and 3.13**. A new conda solve, macOS and
Windows/WSL were not tested locally.
