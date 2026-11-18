# KpPhageClassifier usage guide

[Back to overview](../README.md)

KpPhageClassifier compares a new phage genome with a fixed Klebsiella pneumoniae
phage reference library. It reports candidate whole-proteome wGRR groups (W),
existing similarity subgroups (S), and agreement with reference gene order.

This research tool projects queries onto existing groups; it does not recluster
the reference collection or automatically create new subgroups. Group labels are
local research labels, not taxonomic or host predictions. See the [methods and interpretation guide](METHODS.md)
for thresholds, confidence rules and the limits of gene-order comparisons.

## Installation

Download or clone this repository, then enter its top-level directory.
The recommended Linux environment includes the tested MMseqs2 version:

```bash
conda env create -f environment.yml
conda activate kpphageclassifier
```

Alternatively, with Python >=3.10 and MMseqs2 already installed on PATH:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[dna]'
```

The dna extra installs Pyrodigal-gv 0.3.2 and Pyrodigal 3.7.1, matching the tested
reference gene-calling versions. The remaining Python dependencies are NumPy,
pandas and Biopython. MMseqs2 is a separate executable; version 18.8cc5c was used
for the reference and validation. The conda environment file pins that version.
Linux is the validated platform; macOS and Windows/WSL have not been tested.

For users supplying existing proteins and gene coordinates, install with

```bash
python -m pip install .
```

Pyrodigal-gv is then optional. Full Pharokka and its annotation databases are not
required. Functional annotations in results are derived from reference families;
this tool does not provide comprehensive de novo functional annotation.

The source-checkout launcher `./kp-phage` is also available. It uses `.venv/bin/python`
when present, otherwise python3; `KP_PHAGE_PYTHON` can override this interpreter.
The installed kp-phage command uses the environment into which pip installed it.

## Install the reference library

```bash
kp-phage install-db
kp-phage doctor --require-gene-calling
```

The default database directory is `~/.cache/kpphageclassifier/reference_db`, or
`$XDG_CACHE_HOME/kpphageclassifier/reference_db` if `XDG_CACHE_HOME` is set.
Use `KP_PHAGE_DB` or `--db` to select a different directory:

```bash
kp-phage install-db --db references/kp_phage_v1.0.0
kp-phage doctor --db references/kp_phage_v1.0.0 --require-gene-calling
```

The [reference library](https://doi.org/10.6084/m9.figshare.35048395.v2) is released
separately on Figshare.

Library version 1.0.0 contains 1,687 unique reference genomes represented by
173,570 proteins, genome/family metadata, 30 W groups and 125 existing wGRR
subgroup labels. Whole-genome DNA is not included in this library. The 106
separately downloaded literature phage-plasmids and chromosomal proviruses are
not reference members. Figshare item version 2 and library version 1.0.0 refer
to different versioning schemes.

The installer downloads approximately 48 MB, checks the pinned SHA256, validates
all packaged files, and installs approximately 106 MB. Allow at least 160 MB plus
spare disk space during installation. An existing database is never replaced.
Interrupted downloads retry up to three times, resuming when supported. If a
command ultimately fails, temporary files are removed; rerun the command to retry.
No Figshare account or token is required.

To install an archive downloaded separately:

```bash
kp-phage install-db --archive /path/to/kp_phage_reference_v1.0.0.tar.gz
```

Other supported source options are `--figshare` (public item ID, article link or
DOI) and `--url` (direct HTTPS URL). They must provide the exact registered archive.
Installed libraries contain `installation.json`, `release_documentation/` and
`INSTALLED_SHA256SUMS.txt`. From inside the database directory, verify with:

```bash
sha256sum -c INSTALLED_SHA256SUMS.txt
```

## Classify a genome

```bash
kp-phage classify --genomes new_phage.fasta --out new_results --threads 8
```

Add `--db references/kp_phage_v1.0.0` if you chose a custom database directory.
A classification run does not need network access. Use a NEW output directory;
existing outputs are never overwritten. MMseqs2 may require several GB of RAM
and temporary disk space; `--memory-limit` defaults to 8G for its split-memory
setting and is not a hard limit on total process memory.

Each DNA FASTA record is treated as a separate query genome. Aggregation of
multiple contigs into one genome is not implemented. Use `--genome-status complete`
only when completeness is independently supported; the default is unknown.
partial is also supported. These are user-supplied flags, not completeness tests.

For an external MMseqs2 or gene-calling interpreter:

```bash
kp-phage classify --genomes new_phage.fasta --out new_results \
  --mmseqs /path/to/mmseqs --caller-python /path/to/python
```

Equivalent defaults can be set with `KP_PHAGE_MMSEQS` and `KP_PHAGE_CALLER_PYTHON`.
The selected caller interpreter must have Biopython, Pyrodigal-gv and Pyrodigal.
In a source checkout only, an ignored `.kp-phage.local.json` can provide db, mmseqs
and `caller_python` defaults. Environment variables override that file; command
options override both. Local settings are not part of the distributed package.

## Use existing proteins

```bash
kp-phage classify --proteins phage.faa --genes genes.tsv --out protein_results
```

The tab-delimited gene table requires these columns:

```text
genome_id	gene_id	start	end	strand
```

Coordinates are 1-based and inclusive; strand is + or -. Gene IDs must be unique
across genomes and match the FASTA protein IDs exactly. A terminal stop symbol is
removed; internal stop symbols are rejected. Use a consistent prediction method
for comparability with the references. Order differences from precalled proteins
remain subject to quality review because no query DNA is available.

## Outputs

| File | Contents |
| --- | --- |
| `REPORT.txt` | Readable query summaries |
| `assignments.tsv` / `assignments.json` | Candidate W/S labels, evidence, closest reference, wGRR, matched fractions and gene-order status |
| `reference_similarity.tsv` | Scores against every eligible reference |
| `group_and_subgroup_rankings.tsv` | W rankings and subgroup ranks within the best W group |
| `query_family_assignments.tsv` | Family matches, ambiguity and coordinates |
| `gene_order_comparisons.tsv` | Signed shared-family order and anchor coverage |
| `reciprocal_best_matches.tsv` | Protein pairs contributing to wGRR |
| `query_genomes.tsv` / `query_genes.tsv` / `query_proteins.faa` | Normalized inputs |
| `run_manifest.json` | Library/input hashes, parameters, versions and exclusions |
| `run_status.json` | Running, complete or failed; detailed errors are logged |

Interpret group confidence as heuristic evidence, not a calibrated probability.
An exact shared-family gene-order match does not imply identical genomes. Known
MTase families are excluded from order comparison; whole-proteome wGRR includes
MTases. Low-support or ambiguous queries can remain unassigned.

## Reference maintenance

Normal users should use `install-db`. The developer-only `build-reference` command
requires the original analysis project's tables, proteins and genome FASTA,
which are not distributed with this software repository:

```bash
kp-phage build-reference --project-root /path/to/original_analysis \
  --db /path/to/new_reference_db
```

Build into a new directory and select it explicitly with `--db`. To hold references
out of a query comparison, repeat `--exclude-reference` with a genome ID or exact
accession. Exclusions apply to scoring, family assignment and order comparisons;
family definitions remain those of the original reference construction.

## Tests and validation

```bash
python -m pip install '.[dna,dev]'
python -m unittest discover -s tests -v
ruff check src tests
python -m build
python -m twine check dist/*
```

Unit tests use synthetic fixtures without downloading the reference library.
See the [validation summary](VALIDATION.md) for local release checks and their scope. The GitHub
workflow runs unit tests and packaging checks; it does not download the full
library or run whole-reference MMseqs searches. Existing held-out tests are
engineering checks, not independent biological validation of all thresholds.

## License and citation

The software is released under the [MIT License](../LICENSE). Dependencies retain
their own licenses. The [reference library on Figshare](https://doi.org/10.6084/m9.figshare.35048395.v2)
is labelled **CC BY 4.0**. The software license does not relicense reference data.

When reporting results, record the software version, reference DOI/version,
MMseqs2 and gene-caller versions and any changed parameters. Cite the reference
library using its Figshare citation. See the [methods and interpretation guide](METHODS.md) for methodological
resources. No software paper or separate software DOI is claimed for this release.
