# KpPhageClassifier

**KpPhageClassifier** assigns a query phage genome to one of 30 reference **W groups**
defined for *Klebsiella pneumoniae*-associated phages. The W groups were identified
by Louvain clustering of a whole-proteome similarity network. To place a new genome,
the classifier compares its proteins with the reference proteins using weighted gene
repertoire relatedness (**wGRR**). It ranks each W group by its best-scoring reference
and reports a candidate assignment when the protein matches support it.

**Version 0.1.2** · **Python 3.10+** · [Documentation](docs/USAGE.md) · [Release notes](docs/CHANGELOG.md)

[Open the W-group explorer](https://raymondliu1124.github.io/KpPhageClassifier/) to browse
the 1,687 reference genomes. Search by phage name or accession to find its W group,
view its whole-proteome similarities and inspect the group's representative gene framework.

| Command | Function |
| --- | --- |
| `classify` | Report a candidate W group, best reference, wGRR and assignment evidence |
| `install-db` | Download and verify the versioned reference library from Figshare |
| `doctor` | Check the reference library, MMseqs2 and optional gene-calling dependencies |
| `build-reference` | Build a new library from the original analysis project; for reference maintainers |

## Installation

Clone the repository and create the Linux environment:

```bash
git clone https://github.com/RaymondLiu1124/KpPhageClassifier.git
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
contains **173,570 predicted proteins from 1,687 reference genomes across 30 W groups**.
It provides reference protein sequences, family memberships, gene coordinates and
W-group assignments. The reference library is distributed through Figshare.

The download is approximately **48 MB** and expands to approximately **106 MB**.
The installer checks the archive and file checksums and preserves existing
libraries. No Figshare account is required.

The default location is `~/.cache/kpphageclassifier/reference_db`, or under
`$XDG_CACHE_HOME` when set. To choose another location, use
`kp-phage install-db --db /path/to/reference_db`, then pass the same `--db` to
`classify` and `doctor`. Offline archive installation is also supported; see the
[usage guide](docs/USAGE.md).

## Test your installation with two genomes

From the repository root, after installing the tool and reference library, classify
the [example FASTA](examples/phage_plasmids/phage_plasmids.fna) and check the output:

```bash
kp-phage classify --genomes examples/phage_plasmids/phage_plasmids.fna --out example_results --threads 8
python examples/phage_plasmids/check_results.py example_results/assignments.tsv
```

The FASTA contains one record for each phage–plasmid. With classifier 0.1.2,
reference library 1.0.0 and default thresholds, the expected results are:

| Query | Accession | W group | Candidate subgroup |
| --- | --- | --- | --- |
| phiKO2 | NC_005857.1 | W005 | W005-S02 |
| pKPHS1 | NC_016838.1 | W020 | W020-S01 |

The checker prints the observed assignments and `PASS` when both rows match the
[expected results](examples/phage_plasmids/expected.tsv). It exits with an error
for a missing query or a different assignment. For a second run, use a fresh
output directory in both commands.

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

Use **`candidate_W_group`** in `assignments.tsv` as the project classification,
and read **`W_evidence`** alongside it. A blank assignment means the query remains
unresolved. A query can still have a closest reference in `best_reference` in that
case. The W labels describe protein-repertoire similarity among the project references.

| File | What to look for |
| --- | --- |
| `REPORT.txt` | Readable summary for each query |
| `assignments.tsv` / `assignments.json` | Candidate W group, evidence, best reference, wGRR and matched-protein fractions |
| `reference_similarity.tsv` | Pairwise scores for checking how well other references match the query |
| `reciprocal_best_matches.tsv` | Protein pairs contributing to each wGRR score |
| `run_manifest.json` / `run_status.json` | Inputs, library version, run parameters and status |

The [methods guide](docs/METHODS.md) documents the W-group support thresholds and
the additional comparisons produced by the software.

## Documentation

The [usage guide](docs/USAGE.md) covers installation alternatives, input
requirements, database management, outputs, configuration and reference maintenance.

- [W-group assignment thresholds and methods](docs/METHODS.md)
- [Validation and its scope](docs/VALIDATION.md)
- [Release notes](docs/CHANGELOG.md)
- [Report an issue](https://github.com/RaymondLiu1124/KpPhageClassifier/issues)

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
