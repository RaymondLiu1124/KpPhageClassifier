import argparse,io,json,tarfile,tempfile,unittest
from contextlib import contextmanager,redirect_stdout
from pathlib import Path
from unittest.mock import patch
from kp_phage_classifier import database as db

class DatabaseTests(unittest.TestCase):
    def fixture(self,base):
        name='kp_phage_reference_v1.0.0';package=base/name;ref=package/'reference_db';ref.mkdir(parents=True)
        for n in db.CORE_FILES:(ref/n).write_text('fixture '+n+'\n')
        (ref/'mmseqs').mkdir()
        for n in db.MMSEQS_FILES:(ref/'mmseqs'/n).write_bytes(b'fixture')
        meta={'schema_version':1,'reference_release':'1.0.0','reference_genomes':2,'reference_proteins':3,'files':{n:db.sha256(ref/n) for n in db.CORE_FILES}}
        (ref/'reference.json').write_text(json.dumps(meta))
        for n in ['README.txt','DATA_DICTIONARY.txt','RELEASE.json']:(package/n).write_text('fixture')
        (package/'SHA256SUMS.txt').write_text(''.join(f'{db.sha256(p)}  {p.relative_to(package).as_posix()}\n' for p in sorted(package.rglob('*')) if p.is_file()))
        archive=base/(name+'.tar.gz')
        with tarfile.open(archive,'w:gz') as tar:tar.add(package,arcname=name)
        release={'filename':archive.name,'package_root':name,'archive_bytes':archive.stat().st_size,'expanded_bytes':sum(p.stat().st_size for p in package.rglob('*') if p.is_file()),'sha256':db.sha256(archive),'reference_genomes':2,'reference_proteins':3,'download_url':None,'figshare_item':None}
        registry=base/'registry.json';registry.write_text(json.dumps({'1.0.0':release}))
        return archive,release,registry
    def args(self,base,archive=None,**kw):
        values=dict(reference_version='1.0.0',db=base/'installed',archive=archive,url=None,figshare=None);values.update(kw)
        return argparse.Namespace(**values)
    def test_figshare_identifiers(self):
        cases={'35048395':'https://api.figshare.com/v2/articles/35048395',
            '10.6084/m9.figshare.35048395.v2':'https://api.figshare.com/v2/articles/35048395/versions/2',
            'https://doi.org/10.6084/m9.figshare.35048395.v2':'https://api.figshare.com/v2/articles/35048395/versions/2',
            'https://figshare.com/articles/dataset/title/35048395':'https://api.figshare.com/v2/articles/35048395',
            'https://figshare.com/articles/dataset/title/35048395/2':'https://api.figshare.com/v2/articles/35048395/versions/2',
            'https://api.figshare.com/v2/articles/35048395/versions/2':'https://api.figshare.com/v2/articles/35048395/versions/2'}
        for value,expected in cases.items():self.assertEqual(db.figshare_api_url(value),expected)
        for value in ['https://figshare.com/account/articles/123','https://figshare.com/projects/title/123','https://example.com/articles/123']:
            with self.assertRaises(ValueError):db.figshare_api_url(value)
    def test_https_required(self):
        for value in ['http://example.com/a','file:///tmp/a','https://user:secret@example.com/a']:
            with self.assertRaises(ValueError):db.check_https(value)
    def test_figshare_selects_expected_archive(self):
        release={'filename':'library.tar.gz','archive_bytes':10}
        result={'files':[{'name':'README.txt','size':3,'download_url':'https://example.org/readme'}, {'name':'library.tar.gz','size':10,'download_url':'https://example.org/archive'}]}
        with patch.object(db,'json_get',return_value=result):self.assertEqual(db.resolve_figshare('123',release),'https://example.org/archive')
        result['files'][1]['size']=11
        with patch.object(db,'json_get',return_value=result),self.assertRaises(ValueError):db.resolve_figshare('123',release)
    def test_offline_install_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);archive,release,registry=self.fixture(base);a=self.args(base,archive)
            with patch.object(db,'REGISTRY',registry),redirect_stdout(io.StringIO()):db.install_database(a)
            self.assertTrue((a.db/'reference.json').is_file())
            self.assertTrue((a.db/'INSTALLED_SHA256SUMS.txt').is_file())
            self.assertTrue(json.loads((a.db/'installation.json').read_text())['all_package_checksums_verified'])
            original=db.sha256(a.db/'reference.json')
            with patch.object(db,'REGISTRY',registry),self.assertRaises(ValueError):db.install_database(a)
            self.assertEqual(original,db.sha256(a.db/'reference.json'))
            self.assertFalse(list(base.glob('.installed.install*')))
    def test_corrupt_archive_does_not_install(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);archive,release,registry=self.fixture(base);archive.write_bytes(archive.read_bytes()+b'bad');a=self.args(base,archive)
            with patch.object(db,'REGISTRY',registry),self.assertRaises(ValueError):db.install_database(a)
            self.assertFalse(a.db.exists());self.assertFalse(list(base.glob('.installed.install*')))
    def test_reject_links_and_traversal(self):
        for name,kind in [('root/../../escape','file'),('/tmp/escape','file'),('root/link','link')]:
            with tempfile.TemporaryDirectory() as td:
                base=Path(td);archive=base/'bad.tar.gz'
                with tarfile.open(archive,'w:gz') as tar:
                    item=tarfile.TarInfo(name)
                    if kind=='link':item.type=tarfile.SYMTYPE;item.linkname='/tmp'
                    else:item.size=1
                    tar.addfile(item,None if kind=='link' else io.BytesIO(b'x'))
                with self.assertRaises(ValueError):db.unpack(archive,base/'out',{'package_root':'root','expanded_bytes':10})
    def test_checksums_cover_every_file(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);archive,release,registry=self.fixture(base);package=base/release['package_root']
            (package/'extra.txt').write_text('not checksummed')
            with self.assertRaises(ValueError):db.verify_package(package,'1.0.0',release)
    def test_download_then_install(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);archive,release,registry=self.fixture(base);data=archive.read_bytes()
            class Response(io.BytesIO):headers={'Content-Length':str(len(data))}
            @contextmanager
            def response(url,**kwargs):yield Response(data)
            a=self.args(base,url='https://example.org/library.tar.gz')
            with patch.object(db,'REGISTRY',registry),patch.object(db,'open_https',response),redirect_stdout(io.StringIO()):db.install_database(a)
            self.assertEqual(db.sha256(a.db/'proteins.faa'),db.sha256(base/release['package_root']/'reference_db/proteins.faa'))
    def test_truncated_download(self):
        with tempfile.TemporaryDirectory() as td:
            class Response(io.BytesIO):headers={}
            @contextmanager
            def response(url,**kwargs):yield Response(b'partial')
            with patch.object(db,'open_https',response),redirect_stdout(io.StringIO()),self.assertRaises(ValueError):
                db.download('https://example.org/archive',Path(td)/'x',{'archive_bytes':100,'filename':'archive'})
    def test_resume_after_interrupted_connection(self):
        data=b'abcdefghijklmnopqrstuvwxyz';calls=[]
        class Interrupted:
            headers={'Content-Length':str(len(data))};status=200
            def __init__(self):self.first=True
            def read(self,n):
                if self.first:self.first=False;return data[:7]
                raise TimeoutError('simulated timeout')
        class Resumed(io.BytesIO):
            status=206;headers={'Content-Length':str(len(data)-7),'Content-Range':f'bytes 7-{len(data)-1}/{len(data)}'}
        @contextmanager
        def response(url,extra_headers=None):
            calls.append(extra_headers)
            yield Interrupted() if len(calls)==1 else Resumed(data[7:])
        with tempfile.TemporaryDirectory() as td:
            target=Path(td)/'download'
            with patch.object(db,'open_https',response),redirect_stdout(io.StringIO()):
                db.download('https://example.org/archive',target,{'filename':'archive','archive_bytes':len(data)})
            self.assertEqual(target.read_bytes(),data)
            self.assertEqual(calls,[{}, {'Range':'bytes=7-'}])
    def test_server_ignoring_range_restarts_safely(self):
        data=b'correct_full_archive'
        class Response(io.BytesIO):status=200;headers={'Content-Length':str(len(data))}
        @contextmanager
        def response(url,extra_headers=None):yield Response(data)
        with tempfile.TemporaryDirectory() as td:
            target=Path(td)/'download';target.write_bytes(b'partial')
            with patch.object(db,'open_https',response),redirect_stdout(io.StringIO()):
                db.download('https://example.org/archive',target,{'filename':'archive','archive_bytes':len(data)})
            self.assertEqual(target.read_bytes(),data)
    def test_persistent_timeout_cleans_staging(self):
        calls=[]
        @contextmanager
        def response(url,**kwargs):
            calls.append(url)
            raise TimeoutError('simulated timeout')
            yield
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);archive,release,registry=self.fixture(base);a=self.args(base,url='https://example.org/archive')
            with patch.object(db,'REGISTRY',registry),patch.object(db,'open_https',response),redirect_stdout(io.StringIO()),self.assertRaises(ValueError):
                db.install_database(a)
            self.assertEqual(len(calls),3);self.assertFalse(a.db.exists());self.assertFalse(list(base.glob('.installed.install*')))
    def test_default_source_is_pinned(self):
        r=json.loads(db.REGISTRY.read_text())['1.0.0']
        self.assertEqual(r['figshare_item'],'10.6084/m9.figshare.35048395.v2')
        self.assertEqual(r['sha256'],'e926539849ab76a2b494ab33a8b31ae0d229726792ea5b1d466ca8bb5d0e7044')

if __name__=='__main__':unittest.main()
