"""Hand off native raw coverage statistics within one verified Suite job."""
import json
import math
from pathlib import Path

from .common import digest
from .versions import load_lock


# These locked releases share the same 15-field format. The 1.5 readers below
# require the stored support to equal the AA option, whose native default is 2.
EQUAL_PAIR_SUPPORT_AA = {'1.5.r2', '1.5.r3', '1.5.r4', '1.5.r6'}
LEGACY_AA = {'1.3.r8', '1.4.r2'} | EQUAL_PAIR_SUPPORT_AA


def validated_values(path, bam):
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 16384:
        raise ValueError('missing or unexpected private coverage.stats file')
    lines = path.read_text().splitlines()
    if len(lines) != 1:
        raise ValueError('expected exactly one fresh native measurement')
    # The producer writes a tab before the first value, even for BAM paths with spaces.
    columns = lines[0].split('\t')
    if len(columns) != 16 or columns[0] != str(bam):
        raise ValueError('coverage.stats BAM path or native field count mismatch')
    values = [float(x) for x in columns[1:]]
    if (not all(math.isfinite(x) for x in values) or min(values) < 0
            or min(values[i] for i in (0, 1, 3, 4, 6, 7, 10, 12)) <= 0
            or values[11] < 2 or not values[11].is_integer()
            or values[9] > values[7] or values[7] > values[10]
            or values[13] != 3 or values[14] != bam.stat().st_size):
        raise ValueError('coverage.stats has invalid values or different BAM size/native parameters')
    return columns[1:]


def seed_stage_measurement(attempt, rows, sample, runtime, verified, config):
    """Consume only the just-completed upstream call, never search earlier attempts."""
    aa = verified.get('aa', {})
    if aa.get('version') not in LEGACY_AA or not verified.get('ready'):
        return None
    builtin = next((x['aa'] for x in load_lock()['versions'] if x['id'] == runtime['id']), {})
    if (aa.get('manifest_sha256') != builtin.get('manifest_sha256') or aa.get('commit') != builtin.get('commit')
            or aa.get('version') != builtin.get('version')
            or not aa.get('source_root') or Path(aa['source_root'])/builtin['entrypoint'] != Path(runtime['aa_script']).resolve()
            or aa.get('python') != str(Path(runtime['aa_python']).absolute())):
        return None
    if not rows or any(r.get('status') not in ('retained', 'not_retained')
                       or r.get('sample_id') != sample['sample_id'] or r.get('version') != runtime['id']
                       for r in rows):
        return None
    expected = {seed['seed_id'] for seed in sample['seeds']}
    if len(rows) != len(expected) or {r.get('seed_id') for r in rows} != expected:
        return None
    private_ref = attempt / 'private_AA_DATA_REPO'
    stats = private_ref / 'coverage.stats'
    bam = Path(sample['alignment']).resolve()
    helper = Path(runtime['aa_script']).resolve().with_name('amplified_intervals.py')
    build = sample['reference_build']
    try:
        identity = next((item for item in sample['_input_identities'] if item['path'] == str(bam)), {})
        stat = bam.stat()
        if stat.st_size != identity.get('bytes') or stat.st_mtime_ns != identity.get('mtime_ns'):
            raise ValueError('BAM metadata changed since bounded input validation')
        provenance = attempt / 'upstream_provenance.json'
        if any(r.get('provenance') != str(provenance) for r in rows):
            raise ValueError('upstream records reference a different producer')
        recorded = json.loads(provenance.read_text())
        args = recorded['amplified_intervals']['args']
        expected_no_cstats = config.get('parameter_policy') != 'suite_defaults'
        if (recorded['version'] != runtime['id'] or recorded['io_policy'] != config['_io_policy_path']
                or recorded['amplified_intervals']['sha256'] != digest(helper)
                or Path(recorded['amplified_intervals']['path']).resolve() != helper
                or ('--no_cstats' in args) != expected_no_cstats or args[args.index('--ref') + 1] != build
                or Path(args[args.index('--bam') + 1]).resolve() != bam
                or (private_ref/build).resolve(strict=True) != (Path(config['reference']['aa_data_repo'])/build).resolve(strict=True)):
            raise ValueError('upstream native measurement identity mismatch')
        with (attempt / 'seed_processing.log').open() as log:
            command = json.loads(log.readline())['argv']
        if (command[0] != runtime['aa_python'] or command[command.index('--script') + 1] != str(helper)
                or command[command.index('--policy') + 1] != config['_io_policy_path']
                or command[command.index('--') + 1:] != args):
            raise ValueError('seed filter command differs from recorded producer')
        fields = validated_values(stats, bam)
        return dict(mode='same_version_seed_stage', path=str(stats), sha256=digest(stats), fields=fields,
                    producer=str(provenance), producer_sha256=digest(provenance),
                    helper_sha256=digest(helper), aa_manifest_sha256=aa['manifest_sha256'],
                    reference_build=build)
    except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
        return dict(mode='native_measurement', reason=str(exc))


def plan_background_measurement(measurement, aa_version, bam, job, *, allow_alias):
    """Decide whether to stage a verified candidate for the locked AA reader."""
    if not measurement:
        return dict(mode='native_measurement', reason='No eligible same-version seed-stage measurement')
    if measurement['mode'] != 'same_version_seed_stage':
        return measurement

    pair_support = int(round(float(measurement['fields'][11])))
    bam_has_whitespace = any(c.isspace() for c in str(Path(bam).resolve()))
    if aa_version in EQUAL_PAIR_SUPPORT_AA and pair_support != 2:
        reason = f'AA {aa_version} requires cached pair support to equal pair_support_min=2; measured {pair_support}'
    elif bam_has_whitespace and not allow_alias:
        reason = 'Native AA cannot match whitespace in the original BAM path in coverage.stats'
    elif bam_has_whitespace and any(c.isspace() for c in str(Path(job).absolute())):
        reason = 'BAM and working paths both contain whitespace; native AA cannot match a private cache alias'
    else:
        return measurement
    candidate = {k: measurement[k] for k in ('path', 'sha256', 'producer') if k in measurement}
    return dict(mode='native_measurement', reason=reason,
                rejected_candidate=dict(candidate, pair_support=pair_support))


def native_background_recomputed(stdout, native_log):
    marker = 'cstats not found, generating coverage statistics'
    if marker in stdout:
        return True
    # Native AA writes this at DEBUG level, which stdout may exclude.
    try:
        with native_log.open(errors='replace') as log:
            return any(marker in line for line in log)
    except OSError:
        return False


def stage_measurement(measurement, attempt, bam, *, allow_alias=True):
    """Stage AA's own file format, preserving all numeric tokens verbatim."""
    original = Path(bam).resolve()
    fields = validated_values(Path(measurement['path']), original)
    if fields != measurement['fields'] or digest(measurement['path']) != measurement['sha256']:
        raise ValueError('Native AA measurement changed after the seed stage')
    # Historical AA splits its native cache on whitespace and calls getsize
    # on the BAM filename. This transient path points to the original BAM.
    native_name = original
    if allow_alias and any(c.isspace() for c in str(original)):
        native_name = attempt / 'coverage_input.bam'
        if any(c.isspace() for c in str(native_name)):
            raise ValueError('Historical AA needs a whitespace-free working directory')
        native_name.symlink_to(original)
    (attempt / 'reference' / 'coverage.stats').write_text(
        '\t'.join([str(native_name), *fields]) + '\n')
    return str(native_name)
