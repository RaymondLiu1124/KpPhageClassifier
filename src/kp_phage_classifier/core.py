"""Protein matching and reference-relative assignment. Scores are not taxonomy."""
import numpy as np
import pandas as pd

ALIGN_COLS=['query','target','fident','alnlen','qlen','tlen','qcov','tcov','evalue','bits']

def reciprocal_matches(forward,reverse,queries,references):
    def prepare(frame,rev=False):
        f=frame[(frame.fident>=.35)&(frame.qcov>=.5)&(frame.tcov>=.5)&(frame.evalue<1e-4)].copy()
        f=f.rename(columns={'target' if rev else 'query':'q','query' if rev else 'target':'r'})
        f=f[f.q.isin(queries.index)&f.r.isin(references.index)]
        f['qg']=f.q.map(queries.genome_id);f['rg']=f.r.map(references.genome_id)
        f['mincov']=f[['qcov','tcov']].min(axis=1)
        keys=['r','qg'] if rev else ['q','rg'];tie='q' if rev else 'r'
        return f.sort_values(['bits','evalue','fident','mincov','alnlen',tie],ascending=[False,True,False,False,False,True]).drop_duplicates(keys)
    a=prepare(forward);b=prepare(reverse,True)
    out=a.merge(b[['q','r','fident']],on=['q','r'],suffixes=('_forward','_reverse'))
    out['identity']=(out.fident_forward+out.fident_reverse)/2
    assert not out.duplicated(['q','rg']).any() and not out.duplicated(['r','qg']).any()
    return out[['q','r','qg','rg','identity']]

def classify_ranked(best,runner):
    score=float(best.wgrr);n=int(best.reciprocal_pairs)
    qf=float(best.query_matched_fraction);rf=float(best.reference_matched_fraction)
    margin=(score-runner)/score if score else 0.
    if score<=.05:status='no_reference_edge'
    elif margin<.10:status='ambiguous'
    elif n<3 or min(qf,rf)<.20:status='limited_support'
    elif score>=.50 and min(qf,rf)>=.50:status='strong_similarity'
    else:status='moderate_similarity'
    return status,margin

def rank_groups(pairs,column):
    ranks=pairs.sort_values(['wgrr','reference_genome_id'],ascending=[False,True]).drop_duplicates(column).copy()
    ranks['rank']=range(1,len(ranks)+1)
    return ranks

def adjacency(seq):
    if len(seq)<2:return set()
    return {min((a,sa,b,sb),(b,-sb,a,-sa)) for (a,sa),(b,sb) in zip(seq,seq[1:]+seq[:1])}

def order_score(a,b,min_anchors=10,min_fraction=.5):
    # Inputs must already contain only one copy of each included family.
    assert len({f for f,s in a})==len(a) and len({f for f,s in b})==len(b)
    common={f for f,s in a}&{f for f,s in b};k=len(common)
    fa=k/len(a) if a else 0.;fb=k/len(b) if b else 0.
    if k<min_anchors or min(fa,fb)<min_fraction:return k,fa,fb,np.nan
    aa=[x for x in a if x[0] in common];bb=[x for x in b if x[0] in common]
    return k,fa,fb,len(adjacency(aa)&adjacency(bb))/k

def single_copy_order(frame):
    frame=frame.reset_index(drop=True)
    f=frame[(frame.family_status=='assigned')&(~frame.family_excluded_MTase.astype(bool))].copy()
    counts=f.family_id.value_counts();f=f[f.family_id.map(counts)==1].sort_values(['start','end','gene_id'])
    return [(r.family_id,1 if r.strand=='+' else -1) for r in f.itertuples()]

def assign_families(hits,reference_proteins,gene_ids,bits_margin=.05):
    """Match against members of the candidate W group, retaining family ambiguity."""
    h=hits[(hits.fident>=.50)&(hits.qcov>=.80)&(hits.tcov>=.80)&(hits.evalue<1e-5)&hits.target.isin(reference_proteins.index)].copy()
    h['family_id']=h.target.map(reference_proteins.family_id)
    h=h.sort_values(['bits','evalue','fident','target'],ascending=[False,True,False,True]).drop_duplicates(['query','family_id'])
    lookup={q:f for q,f in h.groupby('query',sort=False)};rows=[]
    for gene in gene_ids:
        f=lookup.get(gene)
        r=dict(gene_id=gene,family_id='',family_status='no_family_match',best_matching_protein='',family_excluded_MTase=False,family_bits_margin=np.nan,family_identity=np.nan)
        if f is not None and len(f):
            best=f.iloc[0];margin=(best.bits-f.iloc[1].bits)/best.bits if len(f)>1 else 1.
            r.update(family_id=best.family_id,family_status='assigned' if margin>=bits_margin else 'ambiguous_family',
                best_matching_protein=best.target,family_excluded_MTase=bool(reference_proteins.loc[best.target,'family_excluded_MTase']),family_bits_margin=margin,family_identity=best.fident)
        rows.append(r)
    return pd.DataFrame(rows)
