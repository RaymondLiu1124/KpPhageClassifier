"""Check the W-group and subgroup assignments for the two example genomes."""

import argparse
import csv
from pathlib import Path


def read_rows(path, required):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        rows = {}
        for row in reader:
            accession = row['query_name']
            if accession in rows:
                raise ValueError(f'{path}: duplicate query {accession}')
            rows[accession] = row
        return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('assignments', type=Path, help='assignments.tsv from the example run')
    args = parser.parse_args()
    columns = {'query_name', 'candidate_W_group', 'candidate_subgroup'}
    try:
        expected = read_rows(Path(__file__).with_name('expected.tsv'), columns | {'phage_name'})
        observed = read_rows(args.assignments, columns)
    except (OSError, ValueError) as error:
        parser.exit(1, f'FAIL: {error}\n')

    errors = []
    print('phage_name\taccession\tcandidate_W_group\tcandidate_subgroup')
    for accession, ref in expected.items():
        row = observed.get(accession)
        if row is None:
            errors.append(f'Missing query: {accession}')
            continue
        group, subgroup = row['candidate_W_group'], row['candidate_subgroup']
        print(f"{ref['phage_name']}\t{accession}\t{group or 'unresolved'}\t{subgroup or 'unresolved'}")
        if (group, subgroup) != (ref['candidate_W_group'], ref['candidate_subgroup']):
            errors.append(f"Expected {accession} to have {ref['candidate_W_group']} and {ref['candidate_subgroup']}.")
    for accession in sorted(observed.keys() - expected.keys()):
        errors.append(f'Unexpected query: {accession}')
    if errors:
        parser.exit(1, 'FAIL: ' + '\n'.join(errors) + '\n')
    print('PASS: both W-group and subgroup assignments match.')


if __name__ == '__main__':
    main()
