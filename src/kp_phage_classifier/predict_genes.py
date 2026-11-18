"""Predict viral genes without requiring a full Pharokka installation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from importlib.metadata import version
import json
from pathlib import Path

from Bio import SeqIO
import pyrodigal_gv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fasta', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    args.out.mkdir(parents=True, exist_ok=True)
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)

    def predict(record):
        return record.id, finder.find_genes(str(record.seq))

    with ExitStack() as stack:
        source = stack.enter_context(args.fasta.open())
        gff = stack.enter_context((args.out / 'prodigal-gv_out.gff').open('w'))
        dna = stack.enter_context((args.out / 'prodigal-gv_out_tmp.fasta').open('w'))
        proteins = stack.enter_context((args.out / 'prodigal-gv_out_aas_tmp.fasta').open('w'))
        pool = stack.enter_context(ThreadPoolExecutor(max_workers=args.threads))
        for record_id, genes in pool.map(predict, SeqIO.parse(source, 'fasta')):
            genes.write_gff(gff, sequence_id=record_id, include_translation_table=True)
            genes.write_genes(dna, sequence_id=record_id)
            genes.write_translations(proteins, sequence_id=record_id)
    versions = {name: version(name) for name in ('pyrodigal-gv', 'pyrodigal', 'biopython')}
    (args.out / 'prediction_versions.json').write_text(json.dumps(versions, indent=2) + '\n')


if __name__ == '__main__':
    main()
