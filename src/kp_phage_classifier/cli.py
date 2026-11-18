"""Reference-relative K. pneumoniae phage classification CLI."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,shutil,subprocess,sys,traceback
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from . import __version__
from .database import install_database
from .core import ALIGN_COLS,reciprocal_matches,rank_groups,classify_ranked,assign_families,single_copy_order,order_score

from .config import DEFAULT_DB, DEFAULT_MMSEQS, DEFAULT_CALLER

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def canonical_dna(seq):
    s=str(seq).upper();rc=str(Seq(s).reverse_complement())
    return hashlib.sha256(min(s,rc).encode()).hexdigest()

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def table(df,path):df.to_csv(path,sep='\t',index=False,na_rep='NA')

def load(path):return pd.read_csv(path,sep='\t',keep_default_na=False)

def run_command(command,log):
    with Path(log).open('a') as f:
        f.write(json.dumps([str(c) for c in command])+'\n');f.flush()
        subprocess.run([str(c) for c in command],stdout=f,stderr=subprocess.STDOUT,check=True)

def require_program(path):
    p=shutil.which(str(path))
    if not p:raise ValueError(f'Executable not found: {path}')
    return p

def build_reference(a):
    project=a.project_root.resolve();db=a.db.resolve()
    if db.exists():raise ValueError(f'Reference directory already exists: {db}. Use a new directory to rebuild.')
    src=project/'07_whole_proteome_wgrr';framework=project/'15_W_group_family_frameworks'
    genomes=load(src/'genome_clusters.tsv')[['genome_id','original_genome_id','protein_count','cluster']].rename(columns={'original_genome_id':'accessions','cluster':'W_group'})
    proteins=[];groups=[];quality=[];inputs=[src/'genome_clusters.tsv',src/'all_predicted_proteins.faa']
    for w in sorted(genomes.W_group.unique()):
        base=framework/w
        op=base/'08_similarity_review'/f'{w}_heatmap_genome_order.tsv'
        fp=base/'03_families/gene_family_membership.tsv';qp=base/'02_quality/genome_quality_and_lineages.tsv'
        inputs.extend([op,fp,qp]);sub=load(op)[['genome_id','candidate_subgroup_at_distance_0.40']].rename(columns={'candidate_subgroup_at_distance_0.40':'subgroup'})
        groups.append(sub);quality.append(load(qp)[['genome_id','primary_eligible','quality_flags']])
        f=load(fp);excluded=set(f.loc[f.annotated_methyltransferase.astype(int)==1,'family_id'])
        f['family_excluded_MTase']=f.family_id.isin(excluded);f['family_status']='assigned';f['W_group']=w
        proteins.append(f[['gene_id','genome_id','W_group','family_id','representative','start','end','strand','annotation','family_excluded_MTase','family_status']])
    genomes=genomes.merge(pd.concat(groups),on='genome_id',validate='one_to_one').merge(pd.concat(quality),on='genome_id',validate='one_to_one')
    proteins=pd.concat(proteins,ignore_index=True)
    assert not proteins.gene_id.duplicated().any()
    counts=proteins.groupby('genome_id').size()
    assert np.array_equal(genomes.protein_count.to_numpy(int),genomes.genome_id.map(counts).to_numpy(int))
    fa_ids=[r.id for r in SeqIO.parse(src/'all_predicted_proteins.faa','fasta')]
    assert len(set(fa_ids))==len(fa_ids) and set(fa_ids)==set(proteins.gene_id)
    dna_path=project/'02_kp_phage_genomes/kp_phages_unique_1687.fasta';inputs.append(dna_path)
    hashes={r.id:canonical_dna(r.seq) for r in SeqIO.parse(dna_path,'fasta')}
    genomes['canonical_dna_sha256']=genomes.accessions.map(hashes)
    if genomes.canonical_dna_sha256.isna().any():raise ValueError('Reference DNA headers do not match accession manifest.')
    db.mkdir(parents=True);(db/'mmseqs').mkdir()
    table(genomes,db/'genomes.tsv');table(proteins,db/'proteins.tsv')
    shutil.copyfile(src/'all_predicted_proteins.faa',db/'proteins.faa')
    reps=set(proteins.representative)
    SeqIO.write((r for r in SeqIO.parse(db/'proteins.faa','fasta') if r.id in reps),db/'family_representatives.faa','fasta')
    mmseqs=require_program(a.mmseqs)
    run_command([mmseqs,'createdb',db/'proteins.faa',db/'mmseqs/reference','--dbtype','1'],db/'build.log')
    files={str(p.relative_to(project)):digest(p) for p in inputs}
    manifest=dict(schema_version=1,tool_version=__version__,created_utc=datetime.now(timezone.utc).isoformat(),
        reference_genomes=len(genomes),reference_proteins=len(proteins),W_groups=genomes.W_group.nunique(),subgroups=genomes.subgroup.nunique(),
        mmseqs_version=subprocess.check_output([mmseqs,'version'],text=True).strip(),source_project=str(project),source_hashes=files,
        files={n:digest(db/n) for n in ['genomes.tsv','proteins.tsv','proteins.faa','family_representatives.faa']})
    write_json(db/'reference.json',manifest)
    print(f'Built reference library: {len(genomes)} genomes, {len(proteins)} proteins, {manifest["subgroups"]} subgroups.\n{db}')

def read_fasta(path):
    with Path(path).open() as handle:
        return list(SeqIO.parse(handle,'fasta'))

def normalize_dna(path):
    records=read_fasta(path)
    if not records:raise ValueError('Genome FASTA contains no sequences.')
    if len({r.id for r in records})!=len(records):raise ValueError('Genome FASTA IDs must be unique.')
    normalized=[];rows=[]
    for i,r in enumerate(records,1):
        seq=str(r.seq).upper()
        if not seq or set(seq)-set('ACGTRYSWKMBDHVNX'):raise ValueError(f'Invalid DNA sequence in {r.id}. Use ungapped DNA FASTA.')
        gid=f'QG{i:05d}';normalized.append(SeqRecord(Seq(seq),id=gid,description=''))
        rows.append(dict(genome_id=gid,query_name=r.id,length_nt=len(seq),ambiguous_bases=sum(x not in 'ACGT' for x in seq),canonical_dna_sha256=canonical_dna(seq)))
    return normalized,pd.DataFrame(rows)

def validate_proteins(path,genes_path):
    records=read_fasta(path);genes=pd.read_csv(genes_path,sep='\t',keep_default_na=False,dtype={'genome_id':str,'gene_id':str})
    required={'genome_id','gene_id','start','end','strand'}
    if not required<=set(genes):raise ValueError('Gene table requires: genome_id, gene_id, start, end, strand.')
    if not records or genes.empty:raise ValueError('Protein FASTA and gene table must be nonempty.')
    ids=[r.id for r in records]
    if len(set(ids))!=len(ids) or genes.gene_id.duplicated().any():raise ValueError('Gene/protein IDs must be globally unique.')
    if set(ids)!=set(genes.gene_id):raise ValueError('Protein IDs and gene-table gene_id values must match exactly.')
    for c in ['start','end']:
        v=pd.to_numeric(genes[c],errors='raise')
        if (v%1!=0).any():raise ValueError('Coordinates must be integers.')
        genes[c]=v.astype(int)
    if (genes.start<1).any() or (genes.end<genes.start).any() or not set(genes.strand)<={'+','-'}:raise ValueError('Use 1-based inclusive start/end coordinates and + or - strand.')
    if (genes.genome_id.astype(str).str.strip()=='').any():raise ValueError('genome_id cannot be blank.')
    for r in records:
        seq=str(r.seq).upper().rstrip('*')
        if not seq or '*' in seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):raise ValueError(f'Invalid protein sequence: {r.id}')
        r.seq=Seq(seq)
    return records,genes

def prepare_queries(a,out,validated):
    if a.genomes:
        records,genomes=validated
        SeqIO.write(records,out/'queries.fna','fasta');calls=out/'gene_calls';calls.mkdir()
        informative=[r for r in records if set(str(r.seq))&set('ACGT')]
        SeqIO.write(informative,calls/'informative_queries.fna','fasta')
        if informative:
            run_command([require_program(a.caller_python),Path(__file__).with_name('predict_genes.py'),calls/'informative_queries.fna',calls,'--threads',a.threads],out/'prediction.log')
            called=read_fasta(calls/'prodigal-gv_out_aas_tmp.fasta')
        else:
            called=[]
            (out/'prediction.log').write_text('No informative DNA bases; no proteins predicted.\n')
        rows=[];seqs=[]
        lengths=genomes.set_index('genome_id').length_nt
        for i,r in enumerate(called,1):
            parts=r.description.split(' # ');attrs=dict(x.split('=',1) for x in parts[4].split(';') if '=' in x)
            gid=r.id.rsplit('_',1)[0];pid=f'QP{i:07d}';start,end=int(parts[1]),int(parts[2]);partial=attrs.get('partial','')
            rows.append(dict(genome_id=gid,gene_id=pid,original_gene_id=r.id,start=start,end=end,strand='+' if parts[3]=='1' else '-',partial=partial,
                interior_partial=int(partial!='00' and start>3 and end<int(lengths[gid])-2)))
            seqs.append(SeqRecord(Seq(str(r.seq).rstrip('*')),id=pid,description=''))
        genes=pd.DataFrame(rows,columns=['genome_id','gene_id','original_gene_id','start','end','strand','partial','interior_partial'])
    else:
        records,raw=validated
        names=list(dict.fromkeys(raw.genome_id.astype(str)));mapping={g:f'QG{i:05d}' for i,g in enumerate(names,1)}
        genomes=pd.DataFrame([dict(genome_id=q,query_name=g,length_nt='',ambiguous_bases='',canonical_dna_sha256='') for g,q in mapping.items()])
        pmap={r.id:f'QP{i:07d}' for i,r in enumerate(records,1)}
        genes=raw.copy();genes['original_gene_id']=genes.gene_id;genes['gene_id']=genes.gene_id.map(pmap);genes['genome_id']=genes.genome_id.astype(str).map(mapping)
        genes['partial']='unknown';genes['interior_partial']=0
        seqs=[SeqRecord(r.seq,id=pmap[r.id],description='') for r in records]
    SeqIO.write(seqs,out/'query_proteins.faa','fasta')
    counts=genes.groupby('genome_id').size();genomes['protein_count']=genomes.genome_id.map(counts).fillna(0).astype(int)
    genomes['genome_status']=a.genome_status
    table(genes,out/'query_genes.tsv');table(genomes,out/'query_genomes.tsv')
    return genes,genomes

def read_hits(path):
    try:return pd.read_csv(path,sep='\t',dtype={'query':str,'target':str})
    except pd.errors.EmptyDataError:return pd.DataFrame(columns=ALIGN_COLS)

def classify(a):
    if a.threads<1:raise ValueError('--threads must be positive.')
    if not 0<=a.family_bits_margin<=1:raise ValueError('--family-bits-margin must be between 0 and 1.')
    if a.proteins and not a.genes:raise ValueError('--proteins requires --genes.')
    if a.genomes and a.genes:raise ValueError('--genes is only used with --proteins.')
    out=a.out.resolve();db=a.db.resolve()
    if out.exists():raise ValueError(f'Output directory already exists: {out}. Choose a new output directory.')
    meta=json.loads((db/'reference.json').read_text())
    for name,h in meta['files'].items():
        if digest(db/name)!=h:raise ValueError(f'Reference file changed: {name}. Rebuild a separate reference library.')
    refs=load(db/'genomes.tsv');rp=load(db/'proteins.tsv').set_index('gene_id',drop=False)
    excludes=set()
    for x in a.exclude_reference:
        selected=refs[(refs.genome_id==x)|refs.accessions.map(lambda y:x in y.split('|'))]
        if selected.empty:raise ValueError(f'Excluded reference not found: {x}')
        excludes.update(selected.genome_id)
    refs=refs[~refs.genome_id.isin(excludes)].copy();rp=rp[rp.genome_id.isin(refs.genome_id)].copy()
    if refs.empty:raise ValueError('No reference genomes remain.')
    validated=normalize_dna(a.genomes) if a.genomes else validate_proteins(a.proteins,a.genes)
    mmseqs=require_program(a.mmseqs)
    out.mkdir(parents=True)
    a._created_output=True
    write_json(out/'run_status.json',{'status':'running'})
    print('Preparing query proteins...',flush=True)
    genes,queries=prepare_queries(a,out,validated)
    qp=genes.set_index('gene_id',drop=False)
    if len(genes):
        common=['-s','7.5','--max-seqs','200000','--min-seq-id','0.35','-c','0.5','--cov-mode','0','-e','1e-4',
            '--alignment-mode','3','--seq-id-mode','0','--threads',str(a.threads),'--split-memory-limit',a.memory_limit,
            '--format-output',','.join(ALIGN_COLS),'--format-mode','4','--remove-tmp-files','1']
        print('Searching proteins in both directions...',flush=True)
        run_command([mmseqs,'easy-search',out/'query_proteins.faa',db/'mmseqs/reference',out/'query_to_reference.tsv',out/'tmp_forward']+common,out/'search.log')
        run_command([mmseqs,'easy-search',db/'proteins.faa',out/'query_proteins.faa',out/'reference_to_query.tsv',out/'tmp_reverse']+common,out/'search.log')
        forward=read_hits(out/'query_to_reference.tsv');reverse=read_hits(out/'reference_to_query.tsv')
    else:
        forward=reverse=pd.DataFrame(columns=ALIGN_COLS)
    print('Assigning groups and comparing gene order...',flush=True)
    rbh=reciprocal_matches(forward,reverse,qp,rp)
    table(rbh,out/'reciprocal_best_matches.tsv')
    aggregates=rbh.groupby(['qg','rg']).agg(identity_sum=('identity','sum'),reciprocal_pairs=('identity','size'))
    ref_orders={g:single_copy_order(f) for g,f in rp.groupby('genome_id',sort=False)}
    pair_rows=[];rank_rows=[];family_frames=[];order_rows=[];summaries=[]
    for q in queries.itertuples():
        pairs=refs[['genome_id','accessions','protein_count','W_group','subgroup','primary_eligible']].rename(columns={'genome_id':'reference_genome_id','accessions':'reference_accessions','protein_count':'reference_protein_count'}).copy()
        pairs['query_genome_id']=q.genome_id;pairs['query_name']=q.query_name;pairs['query_protein_count']=q.protein_count
        vals=[aggregates.loc[(q.genome_id,g)] if (q.genome_id,g) in aggregates.index else None for g in pairs.reference_genome_id]
        pairs['reciprocal_pairs']=[int(x.reciprocal_pairs) if x is not None else 0 for x in vals]
        sums=np.array([float(x.identity_sum) if x is not None else 0 for x in vals]);den=np.minimum(q.protein_count,pairs.reference_protein_count.to_numpy(int))
        pairs['wgrr']=np.divide(sums,den,out=np.zeros(len(pairs)),where=den>0)
        pairs['query_matched_fraction']=pairs.reciprocal_pairs/q.protein_count if q.protein_count else 0.
        pairs['reference_matched_fraction']=pairs.reciprocal_pairs/pairs.reference_protein_count
        assert (pairs.reciprocal_pairs<=den).all() and (pairs.wgrr<=1.000001).all()
        pair_rows.append(pairs)
        ranks=rank_groups(pairs,'W_group');best=ranks.iloc[0];runner=float(ranks.iloc[1].wgrr) if len(ranks)>1 else 0.
        status,margin=classify_ranked(best,runner)
        if not q.protein_count:status='no_predicted_proteins'
        ranks['ranking_level']='W_group';rank_rows.append(ranks)
        group=best.W_group;subset=pairs[pairs.W_group==group];sranks=rank_groups(subset,'subgroup')
        sbest=sranks.iloc[0];srunner=float(sranks.iloc[1].wgrr) if len(sranks)>1 else 0.
        ss,smargin=classify_ranked(sbest,srunner);sranks['ranking_level']='subgroup';rank_rows.append(sranks)
        supported=status in ['strong_similarity','moderate_similarity']
        subgroup=sbest.subgroup if supported and ss in ['strong_similarity','moderate_similarity'] else ''
        qgenes=genes[genes.genome_id==q.genome_id]
        emptyfam=pd.DataFrame(columns=['gene_id','family_id','family_status','best_matching_protein','family_excluded_MTase','family_bits_margin','family_identity'])
        fam=assign_families(forward[forward['query'].isin(qgenes.gene_id)],rp[rp.W_group==group],qgenes.gene_id,a.family_bits_margin) if len(qgenes) else emptyfam
        fam=qgenes.merge(fam,on='gene_id',validate='one_to_one');fam['candidate_W_context']=group;family_frames.append(fam)
        seq=single_copy_order(fam);comparisons=[]
        if best.wgrr>.05 and len(seq):
            for r in subset.itertuples():
                k,fa,fb,score=order_score(seq,ref_orders[r.reference_genome_id])
                rec=dict(query_genome_id=q.genome_id,query_name=q.query_name,reference_genome_id=r.reference_genome_id,reference_accessions=r.reference_accessions,
                    W_group=group,subgroup=r.subgroup,reference_primary_eligible=int(r.primary_eligible),shared_families=k,query_anchor_fraction=fa,reference_anchor_fraction=fb,
                    order_similarity=score,assessable=int(np.isfinite(score)),wgrr=r.wgrr)
                comparisons.append(rec);order_rows.append(rec)
        within=[r for r in comparisons if r['subgroup']==subgroup and r['assessable'] and r['reference_primary_eligible']]
        exact=[r for r in within if r['order_similarity']>=1-1e-12]
        best_order=max(within,key=lambda r:(r['order_similarity'],r['shared_families'],r['wgrr'])) if within else None
        if not subgroup:order_status='unresolved_subgroup'
        elif not within:order_status='insufficient_order_evidence'
        elif exact:order_status='matches_existing_arrangement'
        elif a.genome_status=='complete' and a.genomes and q.ambiguous_bases==0 and qgenes.interior_partial.sum()==0:order_status='candidate_different_arrangement'
        else:order_status='order_difference_needs_quality_review'
        duplicates=refs.loc[refs.canonical_dna_sha256==q.canonical_dna_sha256,'accessions'].to_list() if q.canonical_dna_sha256 else []
        summary=dict(query_name=q.query_name,query_genome_id=q.genome_id,query_proteins=q.protein_count,genome_status=a.genome_status,
            candidate_W_group=group if supported else '',W_evidence=status,best_W_group=group if best.wgrr>0 else '',best_reference=best.reference_accessions if best.wgrr>0 else '',wgrr=float(best.wgrr),reciprocal_pairs=int(best.reciprocal_pairs),
            query_matched_fraction=float(best.query_matched_fraction),reference_matched_fraction=float(best.reference_matched_fraction),W_relative_margin=margin,
            candidate_subgroup=subgroup,subgroup_evidence=ss if supported else 'unresolved_W_group',subgroup_relative_margin=smargin,
            assigned_family_proteins=int((fam.family_status=='assigned').sum()),ambiguous_family_proteins=int((fam.family_status=='ambiguous_family').sum()),
            single_copy_non_MTase_anchors=len(seq),order_status=order_status,assessed_subgroup_references=len(within),exact_order_matching_references=len(exact),
            best_order_reference=best_order['reference_accessions'] if best_order else '',best_order_similarity=best_order['order_similarity'] if best_order else None,
            shared_order_families=best_order['shared_families'] if best_order else 0,
            exact_DNA_reference_matches=';'.join(duplicates),new_subgroup_created=False)
        summaries.append(summary)
    summarydf=pd.DataFrame(summaries);table(summarydf,out/'assignments.tsv');write_json(out/'assignments.json',summaries)
    table(pd.concat(pair_rows,ignore_index=True),out/'reference_similarity.tsv')
    table(pd.concat(rank_rows,ignore_index=True),out/'group_and_subgroup_rankings.tsv')
    table(pd.concat(family_frames,ignore_index=True),out/'query_family_assignments.tsv')
    table(pd.DataFrame(order_rows,columns=['query_genome_id','query_name','reference_genome_id','reference_accessions','W_group','subgroup','reference_primary_eligible','shared_families','query_anchor_fraction','reference_anchor_fraction','order_similarity','assessable','wgrr']),out/'gene_order_comparisons.tsv')
    commands=dict(tool_version=__version__,reference_manifest_sha256=digest(db/'reference.json'),reference_library=str(db),
        input_sha256=digest(a.genomes or a.proteins),genes_sha256=digest(a.genes) if a.genes else None,
        input_mode='DNA_with_Pyrodigal_gv_gene_calling' if a.genomes else 'precalled_proteins_and_coordinates',
        mmseqs_version=subprocess.check_output([mmseqs,'version'],text=True).strip(),excluded_reference_ids=sorted(excludes),
        parameters=dict(identity=.35,coverage_both=.50,evalue=1e-4,family_identity=.50,family_coverage_both=.80,family_evalue=1e-5,
            family_bits_margin=a.family_bits_margin,minimum_order_anchors=10,minimum_order_anchor_fraction=.50,exact_order_score=1.0,
            group_score='maximum_reference_wGRR',reference_subgroups='existing_wGRR_heatmap_labels',gene_order='shared_single_copy_non_MTase_families',threads=a.threads),
        created_utc=datetime.now(timezone.utc).isoformat())
    write_json(out/'run_manifest.json',commands)
    lines=['Kp phage reference-relative classification','',f'Reference library: {meta["reference_genomes"]} genomes; {len(refs)} searched after exclusions.',
        'W/subgroup assignments are candidates. Exact shared-family order does not imply identical genomes.',
        'No reference labels or new subgroup IDs are created. Family matching is sequence-based; known MTase families are excluded from order only.','']
    for r in summaries:
        lines.extend([r['query_name'],f'  W group: {r["candidate_W_group"] or "unresolved"} ({r["W_evidence"]})',
            f'  Best reference: {r["best_reference"] or "no match"}; wGRR {r["wgrr"]:.4f}',
            f'  Candidate subgroup: {r["candidate_subgroup"] or "unresolved"}',
            f'  Gene order: {r["order_status"]}; {r["exact_order_matching_references"]} exact shared-family order matches',
            f'  Best order reference: {r["best_order_reference"] or "unassessable"}; shared anchors {r["shared_order_families"]}',''])
    (out/'REPORT.txt').write_text('\n'.join(lines)+'\n')
    write_json(out/'run_status.json',{'status':'complete','queries':len(queries),'reference_groups_modified':False})
    print(summarydf[['query_name','candidate_W_group','W_evidence','candidate_subgroup','order_status']].to_string(index=False))
    print(f'\nReport: {out/"REPORT.txt"}')

def parser():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--version',action='version',version=__version__)
    sub=p.add_subparsers(dest='command',required=True)
    install=sub.add_parser('install-db',help='Download, verify and install a versioned reference library.')
    source=install.add_mutually_exclusive_group()
    source.add_argument('--figshare',help='Public Figshare item link, ID, or DOI.')
    source.add_argument('--url',help='Direct HTTPS archive download URL.')
    source.add_argument('--archive',type=Path,help='Install an already downloaded release archive.')
    install.add_argument('--reference-version',default='1.0.0')
    install.add_argument('--db',type=Path,default=DEFAULT_DB,help='New database directory; never overwrites existing data.')
    install.set_defaults(func=install_database)
    b=sub.add_parser('build-reference',help='Snapshot the current W groups and protein families.')
    b.add_argument('--project-root',type=Path,required=True);b.add_argument('--db',type=Path,default=DEFAULT_DB);b.add_argument('--mmseqs',default=DEFAULT_MMSEQS);b.set_defaults(func=build_reference)
    c=sub.add_parser('classify',help='Predict proteins, compare references, and check gene order.')
    inp=c.add_mutually_exclusive_group(required=True);inp.add_argument('--genomes',type=Path,help='DNA FASTA; each record is a separate query genome.');inp.add_argument('--proteins',type=Path,help='Precalled protein FASTA; requires --genes.')
    c.add_argument('--genes',type=Path,help='TSV: genome_id, gene_id, start, end, strand (1-based inclusive).')
    c.add_argument('--out',type=Path,required=True);c.add_argument('--db',type=Path,default=DEFAULT_DB);c.add_argument('--threads',type=int,default=8)
    c.add_argument('--memory-limit',default='8G');c.add_argument('--mmseqs',default=DEFAULT_MMSEQS);c.add_argument('--caller-python',default=DEFAULT_CALLER)
    c.add_argument('--genome-status',choices=['unknown','complete','partial'],default='unknown',help='User-supplied completeness information, not a completeness estimate.')
    c.add_argument('--family-bits-margin',type=float,default=.05,help='Minimum relative bitscore lead over a different family; default 0.05, exploratory.')
    c.add_argument('--exclude-reference',action='append',default=[],help='Genome ID or exact accession; repeat to exclude held-out references from all scoring.')
    c.set_defaults(func=classify)
    d=sub.add_parser('doctor',help='Check local tools and reference availability.')
    d.add_argument('--db',type=Path,default=DEFAULT_DB)
    d.add_argument('--mmseqs',default=DEFAULT_MMSEQS)
    d.add_argument('--caller-python',default=DEFAULT_CALLER)
    d.add_argument('--require-gene-calling',action='store_true',help='Also require Pyrodigal-gv for DNA input.')
    d.set_defaults(func=doctor)
    return p

def doctor(a):
    caller=shutil.which(str(a.caller_python))
    probe=None
    if caller:
        probe=subprocess.run([caller,'-c',
            'import pyrodigal_gv,pyrodigal; import importlib.metadata as m; '
            'print(m.version("pyrodigal-gv") + " / " + m.version("pyrodigal"))'],
            capture_output=True,text=True,timeout=30)
    manifest=a.db.expanduser()/'reference.json'
    reference_valid=False
    if manifest.is_file():
        try:
            meta=json.loads(manifest.read_text())
            reference_valid=meta.get('schema_version')==1 and all(
                (manifest.parent/name).is_file() for name in
                ['genomes.tsv','proteins.tsv','proteins.faa','family_representatives.faa','mmseqs/reference'])
        except (ValueError,OSError):
            pass
    checks={'python':sys.executable,'mmseqs':shutil.which(str(a.mmseqs)),
        'caller_python':caller,'gene_calling_available':bool(probe and probe.returncode==0),
        'gene_caller_versions':probe.stdout.strip() if probe and probe.returncode==0 else None,
        'reference_directory':str(a.db),'reference_available':reference_valid}
    print(json.dumps(checks,indent=2))
    missing=[]
    if not checks['mmseqs']:missing.append('MMseqs2 (--mmseqs or KP_PHAGE_MMSEQS)')
    if not reference_valid:missing.append('reference library (run install-db or specify --db)')
    if a.require_gene_calling and not checks['gene_calling_available']:
        missing.append('Pyrodigal-gv (install the dna extra or specify --caller-python)')
    if missing:raise ValueError('Unavailable: '+', '.join(missing))

def main():
    a=parser().parse_args()
    try:a.func(a)
    except Exception as e:
        if getattr(a,'_created_output',False) and (a.out/'run_status.json').exists():
            status=json.loads((a.out/'run_status.json').read_text())
            if status.get('status')=='running':write_json(a.out/'run_status.json',{'status':'failed','error':str(e)})
        if getattr(a,'_created_output',False):
            (a.out/'error.log').write_text(traceback.format_exc())
        print(f'ERROR: {e}',file=sys.stderr);sys.exit(1)

if __name__=='__main__':main()
