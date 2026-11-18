import unittest,tempfile,subprocess,sys
from pathlib import Path
from kp_phage_classifier.cli import normalize_dna,validate_proteins

class InputTests(unittest.TestCase):
    def test_invalid_dna_and_duplicate_names(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'dna.fa';p.write_text('>a\nACGT\n>a\nACGT\n')
            with self.assertRaises(ValueError):normalize_dna(p)
            p.write_text('>a\nPROTEIN\n')
            with self.assertRaises(ValueError):normalize_dna(p)
    def test_precalled_coordinate_and_id_validation(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'a.faa';p.write_text('>gene1\nMALW*\n');t=Path(td)/'a.tsv'
            t.write_text('genome_id\tgene_id\tstart\tend\tstrand\ngenome\tgene1\t1\t15\t+\n')
            r,g=validate_proteins(p,t);self.assertEqual(str(r[0].seq),'MALW')
            t.write_text('genome_id\tgene_id\tstart\tend\tstrand\ngenome\tgene2\t1\t15\t+\n')
            with self.assertRaises(ValueError):validate_proteins(p,t)
            t.write_text('genome_id\tgene_id\tstart\tend\tstrand\ngenome\tgene1\t0\t15\t+\n')
            with self.assertRaises(ValueError):validate_proteins(p,t)
    def test_numeric_identifiers_are_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'a.faa';p.write_text('>001\nMALW\n')
            t=Path(td)/'a.tsv';t.write_text('genome_id\tgene_id\tstart\tend\tstrand\n002\t001\t1\t12\t+\n')
            _,genes=validate_proteins(p,t)
            self.assertEqual(genes.iloc[0].gene_id,'001')
            self.assertEqual(genes.iloc[0].genome_id,'002')
    def test_existing_output_is_untouched(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'existing';out.mkdir();status=out/'run_status.json';before='{"status":"running","owner":"someone_else"}\n';status.write_text(before)
            r=subprocess.run([sys.executable,'-m','kp_phage_classifier.cli','classify','--genomes','missing.fa','--out',str(out)],capture_output=True,text=True)
            self.assertNotEqual(r.returncode,0);self.assertIn('already exists',r.stderr);self.assertEqual(status.read_text(),before)

if __name__=='__main__':unittest.main()
