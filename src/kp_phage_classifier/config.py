"""Portable runtime defaults; explicit CLI options override these settings."""
import json
import os
import sys
from pathlib import Path


def runtime_defaults():
    # An optional ignored source-checkout file preserves developer settings.
    # Installed wheels never depend on the original analysis directory.
    root = Path(__file__).resolve().parents[2]
    local = root / '.kp-phage.local.json'
    settings = {}
    if (root / 'pyproject.toml').is_file() and local.is_file():
        settings = json.loads(local.read_text())
    cache = Path(os.environ.get('XDG_CACHE_HOME', Path.home() / '.cache'))
    return (
        Path(os.environ.get('KP_PHAGE_DB', settings.get('db', str(cache / 'kpphageclassifier' / 'reference_db')))).expanduser(),
        os.environ.get('KP_PHAGE_MMSEQS', settings.get('mmseqs', 'mmseqs')),
        os.environ.get('KP_PHAGE_CALLER_PYTHON', settings.get('caller_python', sys.executable)),
    )


DEFAULT_DB, DEFAULT_MMSEQS, DEFAULT_CALLER = runtime_defaults()
