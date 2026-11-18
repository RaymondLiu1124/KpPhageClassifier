# KpPhageClassifier

**KpPhageClassifier** assigns candidate groups and subgroups to new phage genomes
by comparing their proteins with a *Klebsiella pneumoniae* phage reference library.
It uses whole-proteome weighted gene repertoire relatedness (**wGRR**) for group
assignment, then compares the order of shared protein families.

**Version 0.1.2** · **Python 3.10+** · [Documentation](docs/USAGE.md) · [Release notes](docs/CHANGELOG.md)

| Command | Function |
| --- | --- |
| `classify` | Assign candidate W groups and S subgroups, and compare reference gene order |
| `install-db` | Download and verify the versioned reference library from Figshare |
| `doctor` | Check the reference library, MMseqs2 and optional gene-calling dependencies |
| `build-reference` | Build a new library from the original analysis project; for reference maintainers |

## Installation

Clone the repository and create the Linux environment:

```bash
git clone https://github.com/lxd98/KpPhageClassifier.git
cd KpPhageClassifier
conda env create -f environment.yml
conda activate kpphageclassifier
kp-phage --help
```

The environment includes **MMseqs2 18.8cc5c** and installs the tool with
**Pyrodigal-gv 0.3.2** and **Pyrodigal 3.7.1** for DNA input. Full Pharokka and its
annotation databases are not required.

If MMseqs2 is already available on your `PATH`, you can use a Python virtual
environment instead:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[dna]'
```

For existing proteins and coordinates, `python -m pip install .` is sufficient;
MMseqs2 is still required, while the gene caller is optional. See the
[usage guide](docs/USAGE.md) for executable paths and environment overrides.

## Reference library

Install the library once, then check the environment:

```bash
kp-phage install-db
kp-phage doctor --require-gene-calling
```

The [reference library v1.0.0](https://doi.org/10.6084/m9.figshare.35048395.v2)
contains **1,687 reference genomes represented by 173,570 proteins**, with
**30 W groups** and **125 existing wGRR subgroup labels**. It includes protein
sequences, family membership, gene coordinates and group assignments; genome DNA
is not bundled. The library is hosted separately from the software on Figshare.

The download is approximately **48 MB** and expands to approximately **106 MB**.
The installer checks the archive and file checksums and preserves existing
libraries. No Figshare account is required.

The default location is `~/.cache/kpphageclassifier/reference_db`, or under
`$XDG_CACHE_HOME` when set. To choose another location, use
`kp-phage install-db --db /path/to/reference_db`, then pass the same `--db` to
`classify` and `doctor`. Offline archive installation is also supported; see the
[usage guide](docs/USAGE.md).

## Quick usage

Classify a DNA FASTA:

```bash
kp-phage classify \
  --genomes /path/to/new_phage.fasta \
  --out phage_results \
  --threads 8
```

Each FASTA record is treated as **one query genome**. Multiple contigs belonging
to a single genome are not combined automatically. The output directory must
be new. Classification runs locally after the reference library is installed.

Use `--genome-status complete` only when completeness is independently supported.
The default is `unknown`; `partial` is also accepted. This flag supplies context
for gene-order interpretation and does not estimate genome completeness.

To use existing proteins and gene coordinates:

```bash
kp-phage classify \
  --proteins /path/to/phage.faa \
  --genes /path/to/genes.tsv \
  --out protein_results \
  --threads 8
```

The tab-delimited gene table requires `genome_id`, `gene_id`, `start`, `end` and
`strand`. Coordinates are **1-based, inclusive**, with `+` or `-` for strand.
Each `gene_id` must match a protein FASTA ID exactly and be unique across genomes.
Use consistent gene predictions for comparability with the reference library.

## Output

Start with `REPORT.txt` for a readable summary or `assignments.tsv` for one row
per query. The latter reports candidate groups, the closest reference, wGRR,
matched-protein fractions and gene-order status.

| File | Contents |
| --- | --- |
| `assignments.tsv` / `assignments.json` | Candidate W/S labels and supporting evidence |
| `reference_similarity.tsv` | Similarity scores against each eligible reference |
| `group_and_subgroup_rankings.tsv` | W-group ranks and subgroup ranks within the best W group |
| `query_family_assignments.tsv` | Protein-family matches, ambiguity and gene coordinates |
| `gene_order_comparisons.tsv` | Shared-family order similarity and anchor coverage |
| `reciprocal_best_matches.tsv` | Protein pairs contributing to wGRR |
| `run_manifest.json` / `run_status.json` | Input/library provenance, parameters and run status |

**W and S labels are research groups, not formal taxonomy or host predictions.**
The tool projects a query onto existing groups; low-support or ambiguous queries
can remain unassigned. It does not recluster the references or automatically
create a new subgroup.

Gene order is checked using shared, single-copy protein families. An exact order
match describes those shared families, not identical genomes or gene content.
Known methyltransferase families are excluded from the order comparison;
whole-proteome wGRR still includes methyltransferases. See the
[methods and interpretation guide](docs/METHODS.md) for thresholds and order statuses.

## Documentation

The [usage guide](docs/USAGE.md) covers installation alternatives, input
requirements, database management, outputs, configuration and reference maintenance.

- [Similarity thresholds and gene-order interpretation](docs/METHODS.md)
- [Validation and its scope](docs/VALIDATION.md)
- [Release notes](docs/CHANGELOG.md)
- [Report an issue](https://github.com/lxd98/KpPhageClassifier/issues)

Use `kp-phage classify --help`, `kp-phage install-db --help`, or
`kp-phage doctor --help` for command-specific options.

## Citation

When reporting results, cite the
[reference library](https://doi.org/10.6084/m9.figshare.35048395.v2) using its
Figshare citation and record the software version, library version, MMseqs2 and
gene-caller versions, and any changed parameters.

## License

KpPhageClassifier software and documentation are available under the
[MIT License](LICENSE). The reference library is released separately under
[CC BY 4.0 on Figshare](https://doi.org/10.6084/m9.figshare.35048395.v2).
Dependencies retain their own licenses.
