import unittest
import pandas as pd
from kp_phage_classifier.core import single_copy_order

class RegressionTests(unittest.TestCase):
    def test_indexed_reference_gene_table(self):
        f=pd.DataFrame({'gene_id':['p2','p1'],'family_id':['F2','F1'],'family_status':['assigned','assigned'],
            'family_excluded_MTase':[False,False],'start':[20,1],'end':[30,10],'strand':['-','+']}).set_index('gene_id',drop=False)
        self.assertEqual(single_copy_order(f),[('F1',1),('F2',-1)])
    def test_multicopy_and_mtase_exclusion(self):
        f=pd.DataFrame({'gene_id':['p1','p2','p3','p4'],'family_id':['F1','F1','MT','F2'],'family_status':['assigned']*4,
            'family_excluded_MTase':[False,False,True,False],'start':[1,20,40,60],'end':[10,30,50,70],'strand':['+']*4})
        self.assertEqual(single_copy_order(f),[('F2',1)])

if __name__=='__main__':unittest.main()
