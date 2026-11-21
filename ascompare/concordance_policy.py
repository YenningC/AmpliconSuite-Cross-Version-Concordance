"""Bounded, indexed BAM access for the workflow's locked Python callers.

This is an execution contract for audited scientific tools, not a security
sandbox for arbitrary hostile programs. No BAM copying, indexing, hashing in
full, whole-contig read requests, or unbounded read iteration is allowed.
"""
import argparse
from bisect import bisect_left
import atexit
import json
import os
from pathlib import Path
import runpy
import sys

DEFAULT_MAX_UNIQUE_BP = 100_000_000


class ScanForbidden(RuntimeError):
    pass


def guarded_command(python, script, args, policy_path):
    if not Path(policy_path).is_file():
        raise ValueError('A bounded BAM access policy is required')
    return [str(python), str(Path(__file__).resolve()), '--policy', str(policy_path),
            '--script', str(script), '--', *map(str, args)]


def install_guard(policy, audit_path=None):
    import pysam
    alignments = {str(Path(p).resolve()): str(Path(i).resolve())
                  for p, i in policy.get('alignments', {}).items()}
    if not alignments:
        raise ValueError('No approved indexed BAMs in policy')
    limits = {k: int(policy.get(k, v)) for k, v in (
        ('max_fetch_bp', 50000000), ('max_total_fetch_bp', 500000000),
        ('max_total_reads', 20000000), ('max_unique_bp', DEFAULT_MAX_UNIQUE_BP))}
    if any(v <= 0 for v in limits.values()):
        raise ValueError('BAM access budgets must be positive')
    state = dict(policy='indexed_regions_only', queries=0, requested_bp=0,
                 yielded_reads=0, blocked=[], opened=[], native_bedcov_calls=0,
                 read_count_limit_scope='Python fetch only; native bedcov is bounded by BED bp and process timeout', **limits)
    visited = {}
    state['unique_requested_bp'] = 0
    original = pysam.AlignmentFile

    def save():
        if audit_path:
            Path(audit_path).write_text(json.dumps(state, indent=2) + '\n')

    def reject(reason):
        state['blocked'].append(reason)
        save()
        raise ScanForbidden('Whole-genome/unchecked BAM access prohibited: ' + reason)

    def bam_path(value):
        if isinstance(value, bytes):
            value = os.fsdecode(value)
        if not isinstance(value, (str, os.PathLike)):
            reject('BAM must be an approved local file path')
        p = str(Path(value).resolve())
        if p not in alignments or not p.lower().endswith('.bam'):
            reject('Only the explicitly approved indexed BAMs can be opened')
        if not Path(alignments[p]).is_file():
            reject('Existing BAM index missing; automatic indexing is disabled')
        return p

    def charge(chrom, start, end, lengths):
        if chrom is None or start is None or end is None or chrom not in lengths:
            reject('Every read operation requires chromosome, start and end')
        start, end = int(start), int(end)
        # AA's inclusive-end adapter sometimes requests one extra terminal base.
        end = min(end, lengths[chrom])
        if not 0 <= start < end <= lengths[chrom]:
            reject('Invalid region bounds')
        if start == 0 and end == lengths[chrom]:
            reject('Whole-contig read request')
        size = end - start
        if size > limits['max_fetch_bp']:
            reject('Single-region budget exceeded')
        if state['requested_bp'] + size > limits['max_total_fetch_bp']:
            reject('Cumulative region budget exceeded; result must be incomplete')
        previous = visited.setdefault(chrom, [])
        first = bisect_left(previous, (start,))
        if first and previous[first - 1][1] >= start:
            first -= 1
        last = first
        left, right, replaced_bp = start, end, 0
        while last < len(previous) and previous[last][0] <= right:
            a, b = previous[last]
            left, right = min(left, a), max(right, b)
            replaced_bp += b - a
            last += 1
        new_bp = right - left - replaced_bp
        unique = state['unique_requested_bp'] + new_bp
        unique_limit = min(limits['max_unique_bp'], max(1, sum(lengths.values()) * 9 // 10))
        state['effective_max_unique_bp'] = unique_limit
        if unique > unique_limit:
            state['first_rejected_unique_request'] = state.get('first_rejected_unique_request',
                dict(chrom=chrom, start=start, end=end, new_unique_bp=new_bp,
                     projected_unique_bp=unique, accepted_unique_bp=state['unique_requested_bp'],
                     limit_bp=unique_limit))
            reject('Unique genomic coverage budget exceeded; genome tiling is prohibited')
        previous[first:last] = [(left, right)]
        state['unique_requested_bp'] = unique
        state['queries'] += 1
        state['requested_bp'] += size
        return chrom, start, end

    def counted(iterator):
        for read in iterator:
            state['yielded_reads'] += 1
            if state['yielded_reads'] > limits['max_total_reads']:
                reject('Read-count budget exceeded; result must be incomplete')
            yield read

    class IndexedBam:
        def __init__(self, filename, mode='rb', *args, **kwargs):
            p = bam_path(filename)
            if mode not in ('rb', 'r'):
                reject('Alignment writes and CRAM conversion are disabled')
            kwargs['index_filename'] = alignments[p]
            self._file = original(p, mode, *args, **kwargs)
            # Present the approved caller's alias for historical AA cache matching;
            # I/O and the index always use validated canonical paths.
            self.filename = os.fsencode(os.path.abspath(os.fsdecode(filename)))
            if not self._file.has_index():
                reject('BAM index cannot be read')
            self._lengths = dict(zip(self._file.references, self._file.lengths))
            state['opened'].append(p)

        def __getattr__(self, name):
            if name in ('mate', 'head', 'reset', 'seek', 'next', 'write'):
                reject('Unchecked alignment operation: ' + name)
            return getattr(self._file, name)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self._file.close()

        def __iter__(self):
            reject('Unbounded BAM iteration')

        def __next__(self):
            reject('Unbounded BAM iteration')

        def _region(self, args, kwargs):
            if kwargs.get('until_eof'):
                reject('until_eof scanning')
            if kwargs.get('region') is not None:
                reject('Use explicit finite chromosome/start/end arguments')
            chrom = args[0] if args else kwargs.get('contig', kwargs.get('reference'))
            start = args[1] if len(args) > 1 else kwargs.get('start')
            end = args[2] if len(args) > 2 else kwargs.get('stop', kwargs.get('end'))
            if len(args) > 3:
                reject('Unchecked positional fetch options')
            return charge(chrom, start, end, self._lengths)

        def fetch(self, *args, **kwargs):
            chrom, start, end = self._region(args, kwargs)
            return counted(self._file.fetch(chrom, start, end,
                           multiple_iterators=kwargs.get('multiple_iterators', False)))

        def count(self, *args, **kwargs):
            chrom, start, end = self._region(args, kwargs)
            return self._file.count(chrom, start, end,
                                    read_callback=kwargs.get('read_callback', 'all'))

        def count_coverage(self, *args, **kwargs):
            chrom, start, end = self._region(args, kwargs)
            return self._file.count_coverage(chrom, start, end,
                quality_threshold=kwargs.get('quality_threshold', 15),
                read_callback=kwargs.get('read_callback', 'all'))

        def pileup(self, *args, **kwargs):
            chrom, start, end = self._region(args, kwargs)
            allowed = {k: v for k, v in kwargs.items() if k not in
                       ('contig', 'reference', 'start', 'stop', 'end', 'region', 'until_eof')}
            allowed['truncate'] = True
            return self._file.pileup(chrom, start, end, **allowed)

    pysam.AlignmentFile = pysam.Samfile = IndexedBam
    # The shipped CNVkit versions use samtools bedcov through pysam. Validate
    # its finite BED before entering native code; retain the official algorithm.
    original_bedcov = pysam.bedcov

    def bedcov(*args, **kwargs):
        arguments = list(map(str, args))
        bams = [bam_path(a) for a in arguments if a.lower().endswith('.bam')]
        if not bams:
            reject('bedcov lacks an approved BAM')
        beds = [Path(a) for a in arguments if not a.startswith('-') and
                Path(a).is_file() and not a.lower().endswith(('.bam', '.bai', '.csi'))]
        if len(beds) != 1:
            reject('bedcov requires exactly one finite local BED')
        for bam in bams:
            with original(bam, 'rb', index_filename=alignments[bam]) as handle:
                lengths = dict(zip(handle.references, handle.lengths))
                count = 0
                with beds[0].open() as stream:
                    for line in stream:
                        if not line.strip() or line.startswith(('#', 'track', 'browser')):
                            continue
                        fields = line.split()
                        if len(fields) < 3:
                            reject('Invalid bedcov BED')
                        charge(fields[0], int(fields[1]), int(fields[2]), lengths)
                        count += 1
                if not count:
                    reject('Empty bedcov BED')
        state['native_bedcov_calls'] += 1
        return original_bedcov(*args, **kwargs)

    pysam.bedcov = bedcov
    pysam.samtools.bedcov = bedcov

    def prohibited(*args, **kwargs):
        reject('Unrestricted samtools operation/index creation')

    for name in ('view', 'coverage', 'depth', 'mpileup', 'flagstat', 'stats',
                 'index', 'sort', 'merge', 'cat', 'collate', 'reheader', 'quickcheck'):
        if hasattr(pysam, name):
            setattr(pysam, name, prohibited)
        if hasattr(pysam.samtools, name):
            setattr(pysam.samtools, name, prohibited)
    original_idxstats = pysam.idxstats

    def idxstats(*args, **kwargs):
        candidates = [a for a in args if str(a).lower().endswith('.bam')]
        if len(candidates) != 1:
            reject('idxstats requires an approved indexed BAM')
        bam_path(candidates[0])
        return original_idxstats(*args, **kwargs)

    pysam.idxstats = pysam.samtools.idxstats = idxstats
    import subprocess
    original_popen = subprocess.Popen

    def bounded_popen(args, *pos, **kwargs):
        if kwargs.get('shell') or isinstance(args, str):
            reject('Shell execution is not permitted in a BAM-reading caller')
        if any(str(x).lower().endswith(('.bam', '.cram')) for x in args):
            reject('External BAM-reading subprocess is not audited')
        if Path(str(args[0])).name in ('samtools', 'sambamba'):
            reject('External alignment command is not audited')
        return original_popen(args, *pos, **kwargs)

    subprocess.Popen = bounded_popen
    atexit.register(save)
    return state


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--policy', required=True)
    parser.add_argument('--script', required=True)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    # Executing a module file directly otherwise puts ascompare/statistics.py
    # ahead of the standard library statistics module in child tools.
    own_dir = Path(__file__).resolve().parent
    sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != own_dir]
    policy_path = Path(args.policy).resolve()
    policy = json.loads(policy_path.read_text())
    audit = Path.cwd() / ('bam_access.' + str(os.getpid()) + '.json')
    state = install_guard(policy, audit)
    arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
    sys.argv = [args.script, *arguments]
    sys.path.insert(0, str(Path(args.script).resolve().parent))
    try:
        runpy.run_path(args.script, run_name='__main__')
    finally:
        if state['blocked']:
            raise ScanForbidden('Caller caught a forbidden access; output is not complete')


if __name__ == '__main__':
    main()
