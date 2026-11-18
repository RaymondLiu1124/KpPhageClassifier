# KpPhageClassifier methods and interpretation

[Back to overview](../README.md)

## Workflow

1. Call proteins using Pyrodigal-gv `ViralGeneFinder(meta=True)`, matching the
   original Pharokka gene-calling settings, or use the supplied proteins and gene coordinates. Full Pharokka
   database annotation is not rerun: grouping uses sequence alignments, and known
   functional labels come from existing reference protein families.

2. Search proteins in both directions with MMseqs2. Identity >=35%, coverage >=50%
   of both proteins, E<1e-4, sensitivity 7.5, explicit identity mode. Best hits are
   selected separately for each genome pair (not globally across all references).
   Ties: descending bits, ascending E, descending identity, coverage, alignment
   length, then lexical target ID. Internal IDs are zero-padded for stable ties.

3. Keep reciprocal best matches. wGRR is the sum of mean forward/reverse identity
   fractions divided by the smaller predicted proteome size. Both matched-protein
   fractions and match counts are reported to expose short-genome effects.

4. Rank W groups by their maximum reference wGRR, as in the earlier assignment
   pipeline. Score <=0.05: no reference connection. Relative best/runner-up lead
   <10%: ambiguous. Fewer than 3 reciprocal matches or <20% of either proteome
   matched: limited. Score >=0.50 and >=50% of BOTH proteomes matched: strong.
   Other supported cases: moderate. Zero-protein input is explicitly unresolved.

5. Rank the existing S subgroups within the best W group using the same maximum-
   reference approach and checks. This is a candidate nearest-reference projection,
   not insertion into the original hierarchical tree or a rerun of Louvain.

6. Map query proteins to families within the best W group using hits to ALL family
   members (more inclusive than a single representative). Require >=50% identity,
   >=80% coverage of both proteins, E<1e-5. Keep the best hit per family. If the best
   family bitscore leads the next family by <5%, flag an ambiguous family rather
   than force it into an order anchor. The 5% margin is exploratory and configurable.

7. Compare observed family order against every reference in that W group. Use
   shared, single-copy families and exclude families with ANY member already
   annotated as a methyltransferase. Hypothetical families remain included.
   Unknown MTases may remain; the original W/wGRR comparison still includes all
   proteins, including MTases.

8. For each reference pair, require >=10 shared anchors and >=50% of each genome's
   single-copy non-MTase repertoire. Project orders onto shared anchors; compare
   signed neighbour pairs after allowing reverse complement and start-position
   differences. The last-to-first edge is computational closure, not proof that
   a genome is circular. Extra genes/gene losses omitted from this projection do
   not themselves count as rearrangements. Report gene content and coverage too.

## Strict order comparison

An exact score of 1.0 is required to call a match of shared-family order. This
does not imply identical whole genomes, identical total gene content or
identical DNA.
It means the assessed shared single-copy families have the same signed adjacency.

### Order status

- **`matches_existing_arrangement`**: at least one quality-eligible reference in the
  candidate subgroup has exactly matching shared-family order, at adequate coverage.

- **`candidate_different_arrangement`**: assessable eligible subgroup references exist,
  none match exactly, and a DNA query was explicitly marked complete, with no
  ambiguous bases or interior partial gene calls. This is a candidate for review,
  not proof of a new biological lineage or absence of an unsampled arrangement.

- **`order_difference_needs_quality_review`**: comparable differences exist but those
  quality conditions were not established.

- **`insufficient_order_evidence`**: too few comparable anchors/eligible references.

- **`unresolved_subgroup`**: no supported W/S assignment; do not force order assignment.

Order scores against all references in the best W group are retained, including
other S subgroups. A perfect order match can occur in several S subgroups; it does
not override the protein-similarity evidence. Ambiguous family calls and missing
proteins can reduce comparison coverage. Apparent order differences can also
reflect gene-calling or assembly errors and require inspection.

No strict gene-order subgroup catalogue has yet been accepted for the reference
collection. Consequently the tool reports the closest observed arrangements and
potential differences, but does not manufacture new S or O labels. The original
reference groups remain untouched. This version supports reference-relative comparison; gene-order differences
require biological review before new subgroup labels are accepted.

## Implementation resources

- [MMseqs2 documentation and source](https://github.com/soedinglab/MMseqs2)
- [Pyrodigal-gv documentation and source](https://github.com/althonos/pyrodigal-gv)
- [Pyrodigal gene-caller API](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html)

Protein thresholds, family ambiguity margins and reference-assignment rules above
are this project's operational settings, not universal biological boundaries.
