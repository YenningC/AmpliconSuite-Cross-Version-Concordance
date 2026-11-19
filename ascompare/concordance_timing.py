"""Measured invocation timing and explicitly qualified native AA log intervals.

The workflow total is elapsed wall time. Nested wall/subprocess/native rows are
diagnostics and must not be summed with that total. Reused producer intervals
are reference rows, never current execution time.
"""
from contextvars import ContextVar
from datetime import datetime, timezone
from functools import wraps
import json
from pathlib import Path
import re
import time
import uuid
import warnings

from .common import table


_SESSION = ContextVar('concordance_timing_session', default=None)
_PARENT = ContextVar('concordance_timing_parent', default=None)
FIELDS = ('invocation_id', 'event_id', 'parent_event_id', 'sample_id', 'version',
          'stage', 'measurement_kind', 'accounting', 'duration_seconds', 'status',
          'cache_status', 'producer_id', 'started_utc', 'ended_utc', 'log_path',
          'method', 'approximate', 'notes')


def _utc():
    return datetime.now(timezone.utc).isoformat()


class TimingSession:
    """Append each invocation and retain previous measurements on resume."""
    def __init__(self, out):
        self.out = Path(out)
        self.invocation_id = uuid.uuid4().hex

    def __enter__(self):
        self.out.mkdir(parents=True, exist_ok=True)
        # A killed writer may leave a partial final JSON line. Separate it from
        # this invocation, preserving all new events for the next TSV rebuild.
        with (self.out / '.stage_timing_events.jsonl').open('a+b') as handle:
            if handle.tell():
                handle.seek(-1, 2)
                if handle.read(1) != b'\n':
                    handle.write(b'\n')
        self.token = _SESSION.set(self)
        self.total = measure('workflow', accounting='invocation_total')
        self.total.__enter__()
        return self

    def emit(self, row):
        row = {key: str(value) if isinstance(value, Path) else value for key, value in row.items()}
        row['invocation_id'] = self.invocation_id
        path = self.out / '.stage_timing_events.jsonl'
        with path.open('a', encoding='utf-8') as handle:
            handle.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n')

    def flush(self):
        path = self.out / '.stage_timing_events.jsonl'
        rows = []
        with path.open(encoding='utf-8') as handle:
            for number, line in enumerate(handle, 1):
                try:
                    rows.append(json.loads(line))
                except ValueError:
                    warnings.warn(f'Ignoring incomplete timing event at line {number}: {path}')
        destination = self.out / 'stage_timing.tsv'
        temporary = destination.with_suffix('.tsv.tmp')
        table(temporary, rows, FIELDS)
        temporary.replace(destination)

    def __exit__(self, kind, value, traceback):
        try:
            self.total.__exit__(kind, value, traceback)
            self.flush()
        except OSError:
            if kind is None:
                raise
            warnings.warn('Timing could not be saved while handling a workflow failure')
        finally:
            _SESSION.reset(self.token)
        return False


class measure:
    """Measure monotonic elapsed time, including failed and interrupted calls."""
    def __init__(self, stage, measurement_kind='wall', **metadata):
        parent = _PARENT.get() or {}
        self.row = dict(event_id=uuid.uuid4().hex,
            parent_event_id=parent.get('event_id', ''), sample_id=parent.get('sample_id', ''),
            version=parent.get('version', ''), stage=stage, measurement_kind=measurement_kind,
            accounting='nested', status='complete', cache_status='not_applicable',
            producer_id='', log_path='', method='monotonic elapsed', approximate=False, notes='')
        self.row.update(metadata)
        self.session = _SESSION.get()

    @property
    def status(self):
        return self.row['status']

    @status.setter
    def status(self, status):
        self.row['status'] = status

    def mark(self, **metadata):
        self.row.update(metadata)

    def __enter__(self):
        self.row['started_utc'] = _utc()
        self.started = time.perf_counter()
        self.token = _PARENT.set(self.row)
        return self

    def __exit__(self, kind, value, traceback):
        self.row.update(ended_utc=_utc(), duration_seconds=max(0., time.perf_counter() - self.started))
        if kind is not None:
            self.row['status'] = 'interrupted' if issubclass(kind, KeyboardInterrupt) else 'failed'
            self.row['notes'] = (self.row.get('notes', '') + '; ' + kind.__name__).strip('; ')
        _PARENT.reset(self.token)
        if self.session is not None:
            self.session.emit(self.row)
        return False


def timed_call(stage, function, *args, timing=None, **kwargs):
    """Instrument a normal workflow call without changing its result."""
    with measure(stage, **(timing or {})) as clock:
        result = function(*args, **kwargs)
        records = result if isinstance(result, list) else [result]
        statuses = {r.get('status') for r in records if isinstance(r, dict)}
        if statuses & {'failed', 'invalid', 'interrupted', 'blocked', 'missing', 'unsupported'}:
            clock.status = 'failed' if statuses & {'failed', 'invalid', 'blocked'} else 'incomplete'
        elif isinstance(result, dict) and result.get('ready') is False:
            clock.status = 'incomplete'
        if isinstance(result, dict) and 'cache_hit' in result:
            clock.mark(cache_status='hit' if result['cache_hit'] else 'miss',
                producer_id=result.get('run_dir', '') if result['cache_hit'] else '',
                notes='Current invocation elapsed only; cached producer execution time is not repeated')
        return result


def timed_upstream(function):
    @wraps(function)
    def wrapped(sample, runtime, *args, **kwargs):
        return timed_call('local_upstream', function, sample, runtime, *args,
            timing=dict(sample_id=sample['sample_id'], version=runtime['id']), **kwargs)
    return wrapped


def timed_cnv(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        with measure('cnv_cache_or_compute') as clock:
            result = function(*args, **kwargs)
            cache = result[1]
            clock.mark(cache_status='hit' if cache['cache_hit'] else 'miss',
                       producer_id=str(cache.get('path', '')),
                       notes='Current validation/lookup elapsed; producer subprocess rows are recorded only when executed')
            return result
    return wrapped


_TIME = re.compile(r'#TIME\s+([0-9]+(?:\.[0-9]+)?)\s+(.*)')
_COVERAGE = re.compile(r'Time spent counting coverage \(slow, fast\):\s*([\d.]+) s,\s*([\d.]+) s')


def native_aa_intervals(log_path):
    """Parse marker-to-marker intervals without estimating from file mtimes.

    Plot completion has no dedicated upstream marker, so the next amplicon or
    native total is an approximate endpoint. Coverage calls overlap phases.
    Missing/unfinished phases retain empty durations.
    """
    # Stream in one pass: debug coverage markers can number in the millions.
    # Only top-level phases and a single coverage summary are retained.
    rows = []

    def interval(stage, start, end, approximate=False, notes=''):
        rows.append(dict(stage=stage, duration_seconds=end[0]-start[0] if start and end else None,
            status='complete' if start and end else 'unavailable', approximate=approximate,
            notes=notes, method='native #TIME marker difference',
            marker_start_line=start[2] if start else None, marker_end_line=end[2] if end else None))

    background = extension = determined = reconstruction = plotting = None
    previous = None
    monotonic = True
    coverage_counter = None
    coverage_pending = {}
    coverage_seconds = 0.
    coverage_pairs = coverage_starts = coverage_unmatched = 0
    with Path(log_path).open(errors='replace') as handle:
        for number, line in enumerate(handle, 1):
            counter = _COVERAGE.search(line)
            if counter:
                coverage_counter = (float(counter[1]), float(counter[2]), number)
            match = _TIME.search(line)
            if not match:
                continue
            event = (float(match[1]), match[2].strip(), number)
            elapsed, text, _ = event
            if previous is not None and elapsed < previous:
                monotonic = False
            previous = elapsed
            if text.startswith('Initiating bam_to_breakpoint object') and background is None:
                background = event
            elif text.startswith(('Exploring interval:', 'Finding integration sites:')) and extension is None:
                extension = event
            elif text.startswith('Interval sets for amplicons determined') and determined is None:
                determined = event
            is_reconstruction = text.startswith('Reconstructing amplicon')
            is_plot = text.startswith(('Plotting SV View for amplicon', 'Plotting viral view for interval'))
            if plotting is not None and (is_reconstruction or is_plot or text == 'Total Runtime'):
                interval('aa_plotting', plotting, event, True,
                         'Ends at next top-level phase or native total; includes intervening bookkeeping')
                plotting = None
            if is_reconstruction:
                if reconstruction is not None:
                    interval('aa_reconstruction', reconstruction, None, notes='Native completion marker missing')
                reconstruction = event
            elif text.startswith('Completed decomposition') and reconstruction is not None:
                interval('aa_reconstruction', reconstruction, event)
                reconstruction = None
            if is_plot:
                plotting = event
            if text.startswith('Fetching coverage for '):
                name = text[len('Fetching coverage for '):]
                if name in coverage_pending:
                    coverage_unmatched += 1
                coverage_pending[name] = event
                coverage_starts += 1
            elif text.startswith('Obtained coverage list for '):
                name = text[len('Obtained coverage list for '):]
                start = coverage_pending.pop(name, None)
                if start is not None:
                    coverage_seconds += max(0., elapsed - start[0])
                    coverage_pairs += 1
    interval('aa_background_and_input_setup', background, extension or determined, True,
             'Includes input setup between native markers; background-only time is not separately logged')
    interval('aa_seed_extension', extension, determined)
    if reconstruction is not None:
        interval('aa_reconstruction', reconstruction, None, notes='Native completion marker missing')
    if not any(row['stage'] == 'aa_reconstruction' for row in rows):
        interval('aa_reconstruction', None, None, notes='No native reconstruction start/end markers')
    if plotting is not None:
        interval('aa_plotting', plotting, None, True, 'Native plot endpoint missing')
    if not any(row['stage'] == 'aa_plotting' for row in rows):
        interval('aa_plotting', None, None, notes='Plotting absent, skipped, or not logged')
    if not monotonic:
        for row in rows:
            row.update(duration_seconds=None, status='unavailable', notes='Native elapsed markers are not monotonic')
    if coverage_counter is not None:
        slow, fast, line = coverage_counter
        rows.append(dict(stage='aa_coverage_counting', duration_seconds=slow + fast,
            status='complete', approximate=False, method='native cumulative coverage counters',
            notes=f'Overlaps AA phases; diagnostic only; slow_seconds={slow}; fast_seconds={fast}; native log line {line}'))
    elif coverage_pairs and monotonic:
        unfinished = coverage_unmatched + len(coverage_pending)
        rows.append(dict(stage='aa_coverage_counting', duration_seconds=coverage_seconds,
            status='partial' if unfinished else 'complete', approximate=True,
            method='sum of single-pass matched native #TIME coverage intervals',
            notes=f'Overlaps AA phases; diagnostic only; completed_queries={coverage_pairs}; '
                  f'started_queries={coverage_starts}; unfinished_queries={unfinished}; '
                  'rounded markers; partial totals omit unfinished queries'))
    else:
        interval('aa_coverage_counting', None, None, notes='Reliable coverage timing unavailable')
    return rows


def record_native_aa(run_dir, sample_id='', version='', role='current_execution', producer_id=''):
    """Emit AA subphase diagnostics or explicitly noncurrent producer references."""
    session = _SESSION.get()
    if session is None:
        return
    root = Path(run_dir)
    native = sorted((root / 'aa').glob('*.log'))
    log = native[0] if len(native) == 1 else root / 'aa.log'
    try:
        items = native_aa_intervals(log) if log.is_file() else [dict(
            stage='aa_native_phases', duration_seconds=None, status='unavailable',
            approximate=False, notes='Native AA timing log absent; no phase durations inferred')]
    except OSError as exc:
        items = [dict(stage='aa_native_phases', duration_seconds=None, status='unavailable',
            approximate=False, notes='Native AA timing log unreadable: ' + type(exc).__name__)]
    parent = _PARENT.get() or {}
    for item in items:
        notes = item.get('notes', '')
        if item.get('marker_start_line'):
            notes += f'; native log lines {item["marker_start_line"]}-{item["marker_end_line"]}'
        accounting = 'producer_reference' if role == 'producer_reference' else 'overlapping_diagnostic' if item['stage'] == 'aa_coverage_counting' else 'native_diagnostic'
        if role == 'producer_reference':
            notes = 'Reused producer timing; excluded from current execution cost; ' + notes
        session.emit(dict(event_id=uuid.uuid4().hex, parent_event_id=parent.get('event_id', ''),
            sample_id=sample_id, version=version, measurement_kind='native_log', accounting=accounting,
            cache_status='hit' if role == 'producer_reference' else 'miss', producer_id=producer_id,
            started_utc='', ended_utc='', log_path=str(log), method=item.get('method', 'native log'),
            stage=item['stage'], duration_seconds=item['duration_seconds'], status=item['status'],
            approximate=item['approximate'], notes=notes.strip('; ')))
