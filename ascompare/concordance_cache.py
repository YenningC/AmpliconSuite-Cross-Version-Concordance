"""Integrity-checked, local component artifacts; never a native AA partial cache.

An entry becomes visible only after a successful stage and atomic publication.
The producer and local views share immutable output inodes when possible.  Every
use checks all bytes again; raw BAMs and user resources are never linked here.
"""
from contextlib import contextmanager
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tempfile
import uuid

from .common import digest, now

SCHEMA = 1


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def link_or_copy(source, destination):
    """Only call for completed program outputs, never for an input BAM."""
    source, destination = Path(source), Path(destination)
    if source.is_symlink() or not source.is_file():
        raise ValueError('Component artifact must be a regular file: ' + str(source))
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, destination)
    except OSError as exc:
        if exc.errno not in (errno.EXDEV, errno.EPERM, errno.EACCES,
                             errno.EMLINK, errno.ENOTSUP):
            raise
        shutil.copy2(source, destination)


def _relative(value):
    p = PurePosixPath(value)
    if p.is_absolute() or '..' in p.parts or str(p) in ('', '.'):
        raise ValueError('Unsafe component cache artifact path')
    return p


def read_entry(root, component, key):
    """Malformed, partial, or changed entries are misses, never successes."""
    entry = Path(root) / component / key
    try:
        if entry.is_symlink() or not entry.is_dir():
            return None
        manifest_path = entry / 'manifest.json'
        if manifest_path.is_symlink() or not manifest_path.is_file():
            return None
        manifest = json.loads(manifest_path.read_text())
        if not isinstance(manifest, dict):
            return None
        if (manifest.get('schema') != SCHEMA or manifest.get('component') != component
                or manifest.get('key') != key or manifest.get('status') != 'complete'
                or signature(manifest['identity']) != key):
            return None
        if (not isinstance(manifest.get('stage'), dict)
                or manifest['stage'].get('status') != 'complete'
                or manifest['stage'].get('returncode', 0) != 0
                or not isinstance(manifest.get('producer_run_dir'), str)):
            return None
        files = manifest['files']
        if not isinstance(files, dict) or not files:
            return None
        payload = entry / 'payload'
        if payload.is_symlink() or not payload.is_dir():
            return None
        for rel, expected in files.items():
            path = payload.joinpath(*_relative(rel).parts)
            # Reject symlinks in every path component, including directories.
            if any(p.is_symlink() for p in (path, *path.parents) if p != entry.parent):
                return None
            if not path.is_file() or digest(path) != expected:
                return None
        actual = {p.relative_to(payload).as_posix() for p in payload.rglob('*') if p.is_file()}
        if actual != set(files):
            return None
        return dict(path=str(entry), root=str(Path(root).resolve()), **manifest)
    except (OSError, ValueError, TypeError, KeyError):
        return None


@contextmanager
def stage_lock(root, component, key):
    """Serialize producers of one exact signature, including standalone callers."""
    directory = Path(root) / component
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / (key + '.lock')).open('a+') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def publish(root, component, identity, producer, files, stage):
    """Publish a successful stage under its lock; failed evidence stays in attempt."""
    root, producer = Path(root), Path(producer)
    key = signature(identity)
    current = read_entry(root, component, key)
    if current:
        return current
    directory = root / component
    directory.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.' + key + '.', dir=directory))
    try:
        hashes = {}
        for rel in sorted(set(map(str, files))):
            rel = _relative(rel).as_posix()
            source = producer / rel
            destination = temporary / 'payload' / rel
            link_or_copy(source, destination)
            hashes[rel] = digest(destination)
        if not hashes:
            raise ValueError('An empty component cache cannot be published')
        manifest = dict(schema=SCHEMA, component=component, key=key, status='complete',
                        identity=identity, producer_run_dir=str(producer.resolve()),
                        created_utc=now(), files=hashes, stage=stage)
        (temporary / 'manifest.json').write_text(json.dumps(manifest, indent=2,
                                                         allow_nan=False) + '\n')
        # Kept as a small explicit provenance pointer, not used as trusted state.
        (temporary / 'producer.json').write_text(json.dumps(dict(
            producer_run_dir=manifest['producer_run_dir'], component=component,
            key=key, created_utc=manifest['created_utc']), indent=2) + '\n')
        target = directory / key
        if target.exists():
            # Retain changed evidence instead of deleting or repairing it in place.
            target.rename(directory / ('.invalid-' + key + '-' + uuid.uuid4().hex))
        os.replace(temporary, target)
        result = read_entry(root, component, key)
        if result is None:
            raise ValueError('Component artifact changed while its cache was published')
        return result
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def materialize(entry, destination):
    """Create a per-version view after integrity verification, without overwrites."""
    source = Path(entry['path']) / 'payload'
    destination = Path(destination)
    created = []
    try:
        for rel, expected in entry['files'].items():
            path = source.joinpath(*_relative(rel).parts)
            target = destination / rel
            link_or_copy(path, target)
            created.append(target)
            if digest(target) != expected:
                raise ValueError('Component artifact changed during materialization')
    except BaseException:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        raise


def metadata(entry, hit):
    return dict(hit=bool(hit), key=entry['key'], producer_run_dir=entry['producer_run_dir'],
                source_output=str(Path(entry['root']).parent),
                artifact_paths='local_hardlinks_or_cross_filesystem_copies',
                producer_paths_in_logs='retained_verbatim')
