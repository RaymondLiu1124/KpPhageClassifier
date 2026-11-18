# Changelog

[Back to overview](../README.md)

## 0.1.2 — GitHub release preparation

- Added installable Python packaging, the `kp-phage` console command, optional DNA prediction dependencies, a conda environment and a GitHub test workflow.
- Replaced machine-specific executable and database defaults with portable options.
- Preserved development settings in an ignored local configuration file.
- Moved gene calling into the package and use Pyrodigal-gv directly, retaining the same calling settings without requiring full Pharokka.
- Improved `doctor` to check reference availability and optional gene-calling modules.
- Preserved numeric protein and genome IDs, including leading zeros, in input tables.
- Added the MIT software license, release documentation and repository exclusions.
- Formatted the README and documentation as Markdown, with linked guides, command examples and tables.

The reference library, protein similarity thresholds and group labels are unchanged.

## 0.1.1 — Reference installer

- Added pinned Figshare and offline archive installation.
- Added archive/file checksums, restricted extraction and interrupted-transfer retry/resume.
- Added installer tests.

## 0.1.0 — Initial local classifier

- Added reference-relative W/S assignment, family mapping and signed gene-order checks.
