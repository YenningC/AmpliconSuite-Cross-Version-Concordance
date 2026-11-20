"""Memoize only successful, identical native pysam bedcov calls in WGS batches.

Executed inside each CNVkit interpreter (including spawned workers). The caller
checks the BAM and index identity. Unsupported argument shapes pass straight
through. Each consumer's own CNVkit parses the unmodified native output.
"""

import fcntl
from functools import lru_cache
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()


@lru_cache(maxsize=1)
def engine():
    import pysam
    import platform

    root = Path(pysam.__file__).resolve().parent
    # Include native extensions, Python dispatcher code and wheel-vendored
    # shared libraries. Restrict cache sharing to the same host and ABI.
    files = list(root.glob('*.py')) + list(root.glob('*.so')) + list(root.glob('*.dylib'))
    for folder in (root.parent / 'pysam.libs', root / '.dylibs'):
        if folder.exists():
            files.extend(p for p in folder.rglob('*') if p.is_file())
    return [
        platform.node(),
        platform.platform(),
        sys.implementation.cache_tag,
        pysam.__version__,
        pysam.__samtools_version__,
        {str(p.relative_to(root.parent)): sha(p) for p in sorted(files)},
    ]


def install(spec):
    import pysam

    native = pysam.bedcov
    root = Path(spec['cache'])

    def covered(*args, **kwargs):
        # Both locked versions use bedcov(BED, BAM, ...) with options either
        # before or after those paths. Keep option order/values exactly.
        if kwargs != {'split_lines': False} or not all(isinstance(x, str) for x in args):
            return native(*args, **kwargs)
        positional, options = [], []
        it = iter(args)
        for token in it:
            if token in ('--reference', '-Q', '-d'):
                val = next(it, None)
                if val is None:
                    return native(*args, **kwargs)
                if token == '--reference':
                    if os.path.realpath(val) != spec.get('reference'):
                        return native(*args, **kwargs)
                    val = spec['reference_sha256']
                options.append([token, val])
            elif token.startswith('-'):
                return native(*args, **kwargs)
            else:
                positional.append(token)
        if len(positional) != 2 or positional[1] not in spec['alignments']:
            return native(*args, **kwargs)
        key_parts = [
            1,
            spec['environment'],
            engine(),
            spec['alignments'][positional[1]],
            sha(positional[0]),
            options,
        ]
        key = hashlib.sha256(json.dumps(key_parts, sort_keys=True).encode()).hexdigest()
        directory = root / key[:2]
        directory.mkdir(parents=True, exist_ok=True)
        payload, checksum = directory / (key + '.gz'), directory / (key + '.sha256')
        with (directory / (key + '.lock')).open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if checksum.is_file() and payload.is_file() and sha(payload) == checksum.read_text():
                with gzip.open(payload, 'rt', encoding='utf8', newline='') as f:
                    result = f.read()
                hit = True
            else:
                result = native(*args, **kwargs)
                if not isinstance(result, str) or not result:
                    return result
                with tempfile.NamedTemporaryFile(dir=directory, delete=False) as tmp:
                    tmpname = tmp.name
                try:
                    with gzip.open(tmpname, 'wt', encoding='utf8', newline='', compresslevel=1) as f:
                        f.write(result)
                    os.replace(tmpname, payload)
                    checksum.write_text(sha(payload))
                finally:
                    Path(tmpname).unlink(missing_ok=True)
                hit = False
        with open(spec['events'], 'a') as f:
            f.write(json.dumps({'key': key, 'hit': hit, 'alignment': Path(positional[1]).stem}) + '\n')
        return result

    pysam.bedcov = covered


if __name__ in ('__main__', '__mp_main__') and os.environ.get('AS_GLOBAL_DEPTH_SPEC'):
    install(json.loads(Path(os.environ['AS_GLOBAL_DEPTH_SPEC']).read_text()))
if __name__ == '__main__':
    sys.argv = sys.argv[1:]
    # Keep this driver registered as the actual main module while executing the
    # native console script. Spawned CNVkit workers then install the same raw
    # counter cache before running their native function (important on macOS).
    script = Path(sys.argv[0])
    exec(
        compile(script.read_bytes(), str(script), 'exec'),
        {'__name__': '__main__', '__file__': str(script), '__package__': None},
    )
