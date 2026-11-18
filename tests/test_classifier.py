import unittest
import pandas as pd
import numpy as np
from kp_phage_classifier.core import ALIGN_COLS,reciprocal_matches,order_score,assign_families,classify_ranked

class ClassifierTests(unittest.TestCase):
    def test_order_invariance_and_rearrangement(self):
        a=[(f'F{i}',1 if i%3 else -1) for i in range(20)]
        self.assertEqual(order_score(a,a[8:]+a[:8])[3],1)
        self.assertEqual(order_score(a,[(f,-s) for f,s in reversed(a)])[3],1)
        self.assertEqual(order_score(a,a[:4]+[('extra',1)]+a[4:])[3],1)
        b=a[:4]+a[10:15]+a[4:10]+a[15:]
        self.assertLess(order_score(a,b)[3],1)
        self.assertTrue(np.isnan(order_score(a,a[:9])[3]))
        self.assertTrue(np.isnan(order_score(a,a[:10]+[(f'other{i}',1) for i in range(20)])[3]))
    def test_reciprocal_matching_does_not_double_count(self):
        q=pd.DataFrame({'genome_id':['qg','qg']},index=['q1','q2'])
        r=pd.DataFrame({'genome_id':['rg','rg']},index=['r1','r2'])
        def row(a,b,identity,bits):return [a,b,identity,100,100,100,1.,1.,1e-20,bits]
        f=pd.DataFrame([row('q1','r1',.9,100),row('q2','r1',.8,90),row('q2','r2',.6,80)],columns=ALIGN_COLS)
        rev=pd.DataFrame([row('r1','q1',.8,100),row('r1','q2',.8,90),row('r2','q2',.6,80)],columns=ALIGN_COLS)
        got=reciprocal_matches(f,rev,q,r)
        self.assertEqual(len(got),1);self.assertEqual(got.iloc[0].q,'q1');self.assertAlmostEqual(got.iloc[0].identity,.85)
    def test_best_match_is_per_reference_genome(self):
        q=pd.DataFrame({'genome_id':['qg']},index=['q1']);r=pd.DataFrame({'genome_id':['a','b']},index=['r1','r2'])
        f=pd.DataFrame([['q1','r1',.9,100,100,100,1,1,1e-20,100],['q1','r2',.7,100,100,100,1,1,1e-20,80]],columns=ALIGN_COLS)
        rev=f.rename(columns={'query':'target','target':'query'})[ALIGN_COLS]
        self.assertEqual(len(reciprocal_matches(f,rev,q,r)),2)
    def test_ambiguous_family_is_not_forced(self):
        refs=pd.DataFrame({'family_id':['F1','F2'],'family_excluded_MTase':[False,True]},index=['r1','r2'])
        f=pd.DataFrame([['q1','r1',.8,100,100,100,1,1,1e-20,100],['q1','r2',.8,100,100,100,1,1,1e-20,99]],columns=ALIGN_COLS)
        a=assign_families(f,refs,['q1','q2'])
        self.assertEqual(a.iloc[0].family_status,'ambiguous_family');self.assertEqual(a.iloc[1].family_status,'no_family_match')
    def test_group_margin_and_coverage(self):
        r=pd.Series({'wgrr':.8,'reciprocal_pairs':50,'query_matched_fraction':.7,'reference_matched_fraction':.7})
        self.assertEqual(classify_ranked(r,.3)[0],'strong_similarity')
        self.assertEqual(classify_ranked(r,.79)[0],'ambiguous')
        r['reference_matched_fraction']=.1
        self.assertEqual(classify_ranked(r,.3)[0],'limited_support')
    def test_empty_alignments(self):
        empty=pd.DataFrame(columns=ALIGN_COLS)
        q=pd.DataFrame({'genome_id':['qg']},index=['q1']);r=pd.DataFrame({'genome_id':['rg']},index=['r1'])
        self.assertTrue(reciprocal_matches(empty,empty,q,r).empty)

if __name__=='__main__':unittest.main()
