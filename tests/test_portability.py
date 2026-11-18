import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from kp_phage_classifier import config, database
from kp_phage_classifier.cli import doctor


class PortabilityTests(unittest.TestCase):
    def test_installed_cache_and_environment_overrides(self):
        with tempfile.TemporaryDirectory() as td:
            # Simulate an installed package with no source-checkout config.
            module = Path(td) / 'lib' / 'site-packages' / 'kp_phage_classifier' / 'config.py'
            with patch.object(config, '__file__', str(module)), patch.dict(
                os.environ, {'XDG_CACHE_HOME': td}, clear=True
            ):
                db, mmseqs, caller = config.runtime_defaults()
                self.assertEqual(db, Path(td) / 'kpphageclassifier' / 'reference_db')
                self.assertEqual(mmseqs, 'mmseqs')
                self.assertTrue(caller)
                with patch.dict(os.environ, {'KP_PHAGE_DB': td + '/custom', 'KP_PHAGE_MMSEQS': 'custom-mmseqs'}):
                    db, mmseqs, _ = config.runtime_defaults()
                    self.assertEqual(db, Path(td) / 'custom')
                    self.assertEqual(mmseqs, 'custom-mmseqs')

    def test_release_registry_is_packaged(self):
        release = json.loads(database.REGISTRY.read_text())['1.0.0']
        self.assertEqual(release['reference_genomes'], 1687)
        self.assertEqual(len(release['sha256']), 64)

    def test_doctor_requires_reference_and_optional_caller(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            args = argparse.Namespace(db=root, mmseqs='mmseqs', caller_python='python', require_gene_calling=False)
            failed_probe = subprocess.CompletedProcess([], 1, '', 'missing module')
            with patch('kp_phage_classifier.cli.shutil.which', side_effect=lambda x: x), patch(
                'kp_phage_classifier.cli.subprocess.run', return_value=failed_probe
            ), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'reference library'):
                    doctor(args)
                (root / 'reference.json').write_text('{"schema_version": 1}')
                (root / 'mmseqs').mkdir()
                for name in ['genomes.tsv', 'proteins.tsv', 'proteins.faa', 'family_representatives.faa', 'mmseqs/reference']:
                    (root / name).touch()
                doctor(args)  # Precalled proteins do not need a gene caller.
                args.require_gene_calling = True
                with self.assertRaisesRegex(ValueError, 'Pyrodigal-gv'):
                    doctor(args)
