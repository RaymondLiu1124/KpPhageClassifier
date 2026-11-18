"""Download and install pinned reference releases without modifying existing data."""
from contextlib import contextmanager
from datetime import datetime,timezone
from pathlib import Path,PurePosixPath
import hashlib,json,os,re,shutil,tarfile,tempfile
from http.client import IncompleteRead
from urllib.parse import urlparse
from urllib.request import Request,HTTPRedirectHandler,build_opener
from urllib.error import HTTPError,URLError
from . import __version__

REGISTRY=Path(__file__).with_name('reference_releases.json')
CORE_FILES={'genomes.tsv','proteins.tsv','proteins.faa','family_representatives.faa'}
MMSEQS_FILES={'reference','reference.dbtype','reference.index','reference.lookup','reference.source','reference_h','reference_h.dbtype','reference_h.index'}

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def check_https(url):
    u=urlparse(url)
    if u.scheme!='https' or not u.hostname or u.username or u.password:
        raise ValueError('Downloads require an HTTPS URL without embedded credentials.')
    return url

class HTTPSRedirects(HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        check_https(newurl)
        return super().redirect_request(req,fp,code,msg,headers,newurl)

@contextmanager
def open_https(url,extra_headers=None):
    check_https(url)
    headers={'User-Agent':f'KpPhageClassifier/{__version__}','Accept-Encoding':'identity'}
    headers.update(extra_headers or {})
    request=Request(url,headers=headers)
    try:
        with build_opener(HTTPSRedirects()).open(request,timeout=120) as response:
            check_https(response.geturl())
            yield response
    except HTTPError as e:
        if e.code in (429,500,502,503,504):
            raise ConnectionError(f'Temporary server error HTTP {e.code}') from e
        if e.code in (401,403,404):
            raise ValueError(f'Download returned HTTP {e.code}. Use a public Figshare item or file URL; private drafts are not supported.') from e
        raise ValueError(f'Download returned HTTP {e.code}; retry later.') from e
    except URLError as e:
        raise ConnectionError(f'Could not download the reference library: {e.reason}') from e

def json_get(url):
    with open_https(url) as response:
        data=response.read(10*1024*1024+1)
    if len(data)>10*1024*1024:raise ValueError('Figshare metadata response exceeds the allowed size.')
    return json.loads(data)

def figshare_api_url(value):
    """Public item ID, m9 Figshare DOI, public item URL, or API article URL."""
    value=str(value).strip()
    if re.fullmatch(r'[1-9][0-9]*',value):return f'https://api.figshare.com/v2/articles/{value}'
    doi=re.fullmatch(r'(?:https://doi\.org/)?10\.6084/m9\.figshare\.([0-9]+)(?:\.v([0-9]+))?',value,re.I)
    if doi:
        article,version=doi.groups()
    else:
        parsed=urlparse(check_https(value));host=parsed.hostname.lower()
        if host=='api.figshare.com':
            m=re.fullmatch(r'/v2/articles/([0-9]+)(?:/versions/([0-9]+))?/?',parsed.path)
        elif host=='figshare.com' or host.endswith('.figshare.com'):
            m=re.fullmatch(r'/articles/(?:[^/]+/)*([0-9]+)(?:/([0-9]+))?/?',parsed.path)
            # Greedy optional path sections must not consume the article ID as a slug.
            if m:
                parts=parsed.path.strip('/').split('/')
                if len(parts)>2 and parts[-1].isdigit() and parts[-2].isdigit():
                    article,version=parts[-2:];m=None
                    return f'https://api.figshare.com/v2/articles/{article}/versions/{version}'
        else:m=None
        if not m:raise ValueError('Provide a public Figshare item ID, article link, or 10.6084/m9.figshare DOI (not a project/account link).')
        article,version=m.groups()
    return f'https://api.figshare.com/v2/articles/{article}'+(f'/versions/{version}' if version else '')

def resolve_figshare(value,release):
    metadata=json_get(figshare_api_url(value))
    matches=[f for f in metadata.get('files',[]) if f.get('name')==release['filename']]
    if len(matches)!=1:
        raise ValueError(f'Expected exactly one {release["filename"]} in the public Figshare item; found {len(matches)}.')
    record=matches[0]
    if record.get('is_link_only') or not record.get('download_url'):
        raise ValueError('The archive must be an uploaded downloadable file, not a link-only item.')
    if int(record.get('size',release['archive_bytes']))!=release['archive_bytes']:
        raise ValueError('Figshare file size differs from the pinned reference release.')
    return check_https(record['download_url'])

def download_once(url,destination,release):
    expected=release['archive_bytes'];destination=Path(destination)
    total=destination.stat().st_size if destination.exists() else 0
    if total==expected:return
    if total>expected:raise ValueError('Partial download exceeds the expected archive size.')
    headers={'Range':f'bytes={total}-'} if total else {}
    with open_https(url,extra_headers=headers) as response:
        status=getattr(response,'status',200)
        if status==206:
            content_range=response.headers.get('Content-Range','')
            match=re.fullmatch(r'bytes ([0-9]+)-([0-9]+)/([0-9]+)',content_range)
            if not match or tuple(map(int,match.groups()))!=(total,expected-1,expected):
                raise ValueError('Unexpected Content-Range while resuming the archive.')
        elif status==200:
            total=0  # Server ignored Range: safely restart, never append a full response.
        else:raise ValueError(f'Unexpected archive download HTTP status {status}.')
        length=response.headers.get('Content-Length')
        if length and int(length)!=expected-total:raise ValueError('Download Content-Length does not match the reference release.')
        next_progress=(total//(10*1024*1024)+1)*(10*1024*1024)
        with destination.open('ab' if total else 'wb') as output:
            while True:
                block=response.read(1<<18)
                if not block:break
                total+=len(block)
                if total>expected:raise ValueError('Download is larger than the pinned reference archive.')
                output.write(block)
                if total>=next_progress:
                    print(f'  {total/1e6:.1f} / {expected/1e6:.1f} MB',flush=True);next_progress+=10*1024*1024
    if total!=expected:raise ConnectionError('Download ended before the complete archive was received.')

def download(url,destination,release):
    print(f'Downloading {release["filename"]} ({release["archive_bytes"]/1e6:.1f} MB)...',flush=True)
    for attempt in range(1,4):
        try:
            download_once(url,destination,release)
            return
        except (TimeoutError,ConnectionError,IncompleteRead) as e:
            if attempt==3:raise ValueError(f'Download failed after 3 attempts: {e}. No database was installed; rerun to retry.') from e
            print(f'Connection interrupted; retrying/resuming (attempt {attempt+1}/3)...',flush=True)

def relative_path(name):
    p=PurePosixPath(name)
    if not name or p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name or not p.parts:
        raise ValueError(f'Unsafe path in reference package: {name!r}')
    return p

def unpack(archive,stage,release):
    seen=set();size=0;root=release['package_root']
    with tarfile.open(archive,'r:gz') as tar:
        for i,member in enumerate(tar,1):
            if i>10000:raise ValueError('Reference archive contains too many entries.')
            p=relative_path(member.name)
            if p.parts[0]!=root:raise ValueError('Unexpected top-level directory in reference archive.')
            key=p.as_posix()
            if key in seen:raise ValueError('Duplicate archive path.')
            seen.add(key)
            if not member.isdir() and not member.isfile():raise ValueError('Links and special files are not allowed in reference archives.')
            size+=member.size
            if size>release['expanded_bytes']+65536:raise ValueError('Expanded archive exceeds the expected size.')
            target=Path(stage).joinpath(*p.parts)
            if member.isdir():target.mkdir(parents=True,exist_ok=True)
            else:
                target.parent.mkdir(parents=True,exist_ok=True)
                with tar.extractfile(member) as src,target.open('xb') as dst:shutil.copyfileobj(src,dst,1<<20)
                target.chmod(0o644)
    return Path(stage)/root

def verify_package(package,version,release):
    sums=package/'SHA256SUMS.txt'
    if not sums.is_file():raise ValueError('Package is missing SHA256SUMS.txt.')
    covered=set()
    for line in sums.read_text().splitlines():
        match=re.fullmatch(r'([0-9a-fA-F]{64})  (.+)',line)
        if not match:raise ValueError('Malformed internal checksum manifest.')
        expected,name=match.groups();path=relative_path(name);name=path.as_posix()
        if name in covered:raise ValueError('Duplicate entry in checksum manifest.')
        covered.add(name);f=package.joinpath(*path.parts)
        if not f.is_file() or sha256(f)!=expected.lower():raise ValueError(f'Package checksum failed: {name}')
    actual={p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file()}-{'SHA256SUMS.txt'}
    if actual!=covered:raise ValueError('Checksum manifest does not cover every package file.')
    db=package/'reference_db';meta=json.loads((db/'reference.json').read_text())
    if meta.get('schema_version')!=1 or meta.get('reference_release')!=version:
        raise ValueError('Unsupported reference schema or unexpected release version.')
    if set(meta.get('files',{}))!=CORE_FILES:raise ValueError('Unexpected core-file manifest.')
    for name,h in meta['files'].items():
        if sha256(db/name)!=h:raise ValueError(f'Classifier core-file checksum failed: {name}')
    for name in MMSEQS_FILES:
        if not (db/'mmseqs'/name).is_file():raise ValueError(f'Missing MMseqs database file: {name}')
    if meta.get('reference_genomes')!=release['reference_genomes'] or meta.get('reference_proteins')!=release['reference_proteins']:
        raise ValueError('Reference counts differ from the pinned release.')
    return db

def install_database(a):
    releases=json.loads(REGISTRY.read_text());version=a.reference_version
    if version not in releases:raise ValueError(f'Unknown reference version {version}. Available: {", ".join(releases)}')
    release=releases[version];destination=a.db.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError(f'Database destination already exists: {destination}. Choose a new --db directory; existing libraries are not replaced.')
    source_archive=Path(a.archive).expanduser().resolve() if a.archive else None
    if source_archive and not source_archive.is_file():raise ValueError(f'Archive not found: {source_archive}')
    url=None
    if not source_archive:
        if a.url:url=check_https(a.url)
        elif a.figshare:url=resolve_figshare(a.figshare,release)
        elif release.get('download_url'):url=check_https(release['download_url'])
        elif release.get('figshare_item'):url=resolve_figshare(release['figshare_item'],release)
        else:raise ValueError('The default public Figshare URL has not been configured yet. Supply --figshare ITEM_LINK_OR_DOI, --url DIRECT_DOWNLOAD_URL, or --archive LOCAL_TAR_GZ.')
    destination.parent.mkdir(parents=True,exist_ok=True)
    lock=destination.parent/f'.{destination.name}.install.lock'
    try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError:raise ValueError(f'An installation lock already exists: {lock}. Check whether another install is running.')
    os.close(fd)
    try:
        with tempfile.TemporaryDirectory(prefix=f'.{destination.name}.install-',dir=destination.parent) as tmp:
            tmp=Path(tmp);archive=source_archive or tmp/release['filename']
            if not source_archive:download(url,archive,release)
            if archive.stat().st_size!=release['archive_bytes'] or sha256(archive)!=release['sha256']:
                raise ValueError('Archive SHA256/size verification failed; no database was installed.')
            print('Verifying and extracting the reference library...',flush=True)
            package=unpack(archive,tmp/'unpacked',release);payload=verify_package(package,version,release)
            docs=payload/'release_documentation';docs.mkdir()
            for name in ['README.txt','DATA_DICTIONARY.txt','RELEASE.json','SHA256SUMS.txt']:
                shutil.copyfile(package/name,docs/('ARCHIVE_SHA256SUMS.txt' if name=='SHA256SUMS.txt' else name))
            receipt={'reference_version':version,'archive_sha256':release['sha256'],'source':str(source_archive) if source_archive else url,
                'figshare_item':a.figshare or release.get('figshare_item'),'tool_version':__version__,
                'installed_utc':datetime.now(timezone.utc).isoformat(),'all_package_checksums_verified':True}
            (payload/'installation.json').write_text(json.dumps(receipt,indent=2)+'\n')
            (payload/'INSTALLED_SHA256SUMS.txt').write_text(''.join(
                f'{sha256(p)}  {p.relative_to(payload).as_posix()}\n'
                for p in sorted(payload.rglob('*')) if p.is_file()))
            if destination.exists() or destination.is_symlink():raise ValueError('Database destination appeared during installation; refusing to replace it.')
            payload.rename(destination)
    finally:lock.unlink(missing_ok=True)
    print(f'Installed reference library {version}: {destination}\nVerified {release["reference_genomes"]:,} genomes and {release["reference_proteins"]:,} proteins.\nUse classify --db {destination}',flush=True)
