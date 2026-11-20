"""Native Suite WGS CNVkit discovery with provenance-checked reuse."""

import csv
import gzip
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import tempfile
from functools import lru_cache

from .common import atomic_json, digest, now
from .concordance_cache import stage_lock, signature
from .concordance_cnv import _environment, _probe_cnvkit, _sampled_identity, load_builtin_map, LocalFailure
from .concordance_timing import measure
from .runner import execute


def _stamp(path):
    st = Path(path).stat()
    return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns


@lru_cache(maxsize=256)
def _digest_stamped(path, stamp):
    return digest(path)


def _hash(path):
    return _digest_stamped(str(Path(path).resolve()), _stamp(path))


def _alignment_id(path):
    return dict(_sampled_identity(path), stat=_stamp(path))


def _options(words, command, allowed):
    if words.count(command) != 1:
        raise ValueError(f'Expected a single native {command} command')
    args, positional, opts = iter(words[words.index(command) + 1 :]), [], {}
    for token in args:
        if token.startswith('-'):
            if token not in allowed or token in opts:
                raise ValueError(f'Unexpected or repeated {command} option: {token}')
            opts[token] = next(args, '')
        else:
            positional.append(token)
    return positional, opts


def _same_file(a, b):
    return bool(a and b and Path(a).is_file() and Path(b).is_file() and os.path.samefile(a, b))


@lru_cache(maxsize=32)
def _archive_edges(cns, cns_sha, cnr, cnr_sha):
    # A cleanly truncated file can end on a valid CNS row. Prefer re-segmenting
    # the complete CNR to counting a missing chromosome/tail as no discovery.
    observed = {}
    with open(cns) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            c, start, end = row['chromosome'], int(row['start']), int(row['end'])
            prior = observed.get(c, (start, end))
            observed[c] = min(start, prior[0]), max(end, prior[1])
    opener = gzip.open if Path(cnr).suffix == '.gz' else open
    with opener(cnr, 'rt') as stream:
        rows = csv.DictReader(stream, delimiter='\t')
        for row in rows:
            if float(row.get('weight', 1)) > 0 and math.isfinite(float(row['log2'])):
                edges = observed.get(row['chromosome'])
                if not edges or int(row['start']) < edges[0] or int(row['end']) > edges[1]:
                    return False
    return True


def _historical(root, sample, config, version, spec):
    """Accept only explicitly recorded native WGS and CBS output in its run tree."""
    root = Path(root).resolve(strict=True)
    paths = list(root.glob('*_run_metadata.json'))
    if len(paths) != 1:
        raise ValueError('Expected one native run metadata JSON')
    path = paths[0]
    if path.stat().st_size > 4 * 1024**2 or path.is_symlink():
        raise ValueError('Unexpected native run metadata type/size')
    meta = json.loads(path.read_text())
    if not isinstance(meta, dict):
        raise ValueError('Native metadata must be an object')
    if meta.get('cnvkit_version') != version:
        raise ValueError('Recorded CNVkit version differs')
    source = spec is sample.get('source_result')
    expected_suite = sample.get('source_version') if source else spec.get('version')
    if (
        not meta.get('AmpliconSuite-pipeline_version')
        or expected_suite
        and str(meta['AmpliconSuite-pipeline_version']).lstrip('v') != expected_suite.removeprefix('suite-')
    ):
        raise ValueError('Recorded Suite version absent or inconsistent')
    build = config['reference']['build']
    if (
        meta.get('ref_genome') != build
        or source
        and sample.get('source_cnvkit_version')
        and sample['source_cnvkit_version'] != version
        or source
        and sample.get('_source_cnvkit_issue')
    ):
        raise ValueError('Build or source CNVkit identity inconsistent')
    words = shlex.split(meta['cnvkit_cmd'])
    if words.count(';') != 1:
        raise ValueError('Unrecognized batch/segment command history')
    split = words.index(';')
    bam, batch = _options(
        words[:split], 'batch', {'-m', '-r', '-d', '-p', '--normal', '--fasta', '--rscript-path'}
    )
    ratios, seg = _options(words[split + 1 :], 'segment', {'-p', '-m', '-o', '--rscript-path'})
    if (
        len(bam) != 1
        or not _same_file(bam[0], sample['alignment'])
        or batch.get('-m') != 'wgs'
        or len(ratios) != 1
        or seg.get('-m') != 'cbs'
    ):
        raise ValueError('Native tumor WGS/CBS command does not match')
    directory = Path(batch['-d'])
    name = Path(bam[0]).stem
    if Path(seg['-o']) != directory / (name + '.cns') or Path(ratios[0]) != directory / (name + '.cnr'):
        raise ValueError('Recorded raw CNS/CNR does not belong to this WGS batch')
    # Archive relocation is allowed as a unit, never by bare file basename.
    if directory.parent.name != root.name or not directory.name.endswith('_cnvkit_output'):
        raise ValueError('Native batch output is outside supplied Suite run')
    outdir = root / directory.name
    true_normal = sample.get('normal_alignment') if build != 'GRCh38_viral' else None
    fst = Path(config['reference']['fasta'])
    if true_normal:
        if (
            not _same_file(batch.get('--normal'), true_normal)
            or '-r' in batch
            or Path(batch.get('--fasta', '')) != directory / (fst.stem + '_reduced' + fst.suffix)
        ):
            raise ValueError('Paired normal or native reduced FASTA differs')
    elif (
        '--normal' in batch
        or '--fasta' in batch
        or not _same_file(
            batch.get('-r'),
            Path(config['reference']['aa_data_repo']) / build / (build + '_cnvkit_filtered_ref.cnn'),
        )
    ):
        raise ValueError('Suite tumor-only reference differs')
    cns, cnr = outdir / (name + '.cns'), outdir / (name + '.cnr')
    if source and sample.get('baseline_cns') and Path(sample['baseline_cns']).resolve() != cns.resolve():
        raise ValueError('Declared baseline_cns disagrees with native command')
    if not cnr.is_file():
        cnr = Path(str(cnr) + '.gz')
    reason = None
    if cns.is_file():
        expected = spec.get('cns_sha256')
        if expected and _hash(cns) != expected:
            reason = 'Archived CNS checksum differs from supplied cns_sha256'
        if not reason:
            try:
                _check_cns(cns, str(config['reference']['fasta']) + '.fai')
            except (OSError, KeyError, TypeError, ValueError) as exc:
                reason = 'Invalid archived CNS: ' + str(exc)
        if not reason and cnr.is_file() and not _archive_edges(str(cns), _hash(cns), str(cnr), _hash(cnr)):
            reason = 'Archived CNS does not span all positive-weight CNR chromosome extents'
    item = cns if cns.is_file() and not reason else cnr
    if not item.is_file() or not item.resolve().is_relative_to(root) or not item.stat().st_size:
        raise ValueError('Archived raw CNS/CNR missing')
    return dict(
        kind='cns' if item == cns else 'cnr',
        file=str(item),
        sha256=_hash(item),
        metadata=str(path),
        metadata_sha256=_hash(path),
        command=meta['cnvkit_cmd'],
        source_suite=meta['AmpliconSuite-pipeline_version'],
        cns_rejection=reason,
        evidence='native metadata and same current BAM paths; historical dependency and resource checksums not recorded',
    )


def _historical_inputs(sample, config, version):
    records, rejected, seen = [], [], set()
    for spec in [
        sample.get('source_result', {}),
        *sample.get('comparison_results', []),
        *sample.get('global_cnv_results', []),
    ]:
        root = spec.get('result_dir') or spec.get('run_dir')
        if not root or root in seen:
            continue
        seen.add(root)
        try:
            records.append(_historical(root, sample, config, version, spec))
        except (OSError, EOFError, KeyError, TypeError, ValueError) as exc:
            rejected.append(dict(root=str(root), reason=str(exc)))
    for kind in ('cns', 'cnr'):
        found = [r for r in records if r['kind'] == kind]
        if len({r['sha256'] for r in found}) > 1:
            rejected.append(dict(reason='Conflicting same-version archived ' + kind.upper()))
            return None, rejected
        if found:
            return found[0], rejected
    return None, rejected


def _check_cns(cns, fai):
    return _validated_segments(str(cns), _hash(cns), str(fai), _hash(fai))


@lru_cache(maxsize=64)
def _validated_segments(cns, cns_sha, fai, fai_sha):
    lengths = {f[0]: int(f[1]) for f in (x.split() for x in Path(fai).read_text().splitlines())}
    last, count = {}, 0
    with Path(cns).open() as f:
        rows = csv.DictReader(f, delimiter='\t')
        if (rows.fieldnames or [])[:5] != ['chromosome', 'start', 'end', 'gene', 'log2']:
            raise ValueError('Missing native CNS columns')
        for r in rows:
            c, s, e = r['chromosome'], int(r['start']), int(r['end'])
            if not (last.get(c, 0) <= s < e <= lengths.get(c, 0) and math.isfinite(float(r['log2']))):
                raise ValueError('Invalid, overlapping or unsorted native CNS segments')
            last[c], count = e, count + 1
    if not count:
        raise ValueError('Raw global CNS is empty')
    return count


def _reduced(fasta, sizes, out):
    # Same serialized bases as all locked Suite reduce_fasta.py scripts:
    # lexicographically sorted selected chromosome names, no wrapping.
    import pysam

    identity = dict(
        fasta=_hash(fasta), fai=_hash(str(fasta) + '.fai'), sizes=_hash(sizes), format='suite-reduced-v1'
    )
    key = signature(identity)
    root = out / 'reference' / key
    root.mkdir(parents=True, exist_ok=True)
    path = root / (fasta.stem + '_reduced' + fasta.suffix)
    receipt = root / 'complete.json'
    with stage_lock(out, 'reference', key):
        if receipt.exists():
            info = json.loads(receipt.read_text())
            if (
                path.is_file()
                and _hash(path) == info.get('sha256')
                and Path(str(path) + '.fai').is_file()
                and _hash(str(path) + '.fai') == info.get('index_sha256')
            ):
                return path
        names = sorted(
            {line.split()[0].lstrip('>') for line in Path(sizes).read_text().splitlines() if line.strip()}
        )
        if not names:
            raise ValueError('Empty official noAlt contig list')
        tmp = path.with_suffix('.partial')
        with pysam.FastaFile(str(fasta)) as ref, tmp.open('wb') as f:
            for name in names:
                f.write(('>' + name + '\n').encode())
                for i in range(0, ref.get_reference_length(name), 4 << 20):
                    f.write(ref.fetch(name, i, i + (4 << 20)).encode('ascii'))
                f.write(b'\n')
        os.replace(tmp, path)
        pysam.faidx(str(path))
        atomic_json(
            receipt, dict(identity=identity, sha256=_hash(path), index_sha256=_hash(str(path) + '.fai'))
        )
    return path


def _native_cmd(cmd, stage, attempt, env, timeout):
    log = attempt / (stage + '.log')
    with measure(stage, measurement_kind='subprocess', log_path=str(log)):
        rc = execute(list(map(str, cmd)), attempt, env, log, timeout=timeout)
        if rc:
            raise LocalFailure('failed', stage, f'Native CNVkit exited {rc}; see {log}')


def _private_bam(sample, alignment, index, attempt, stem):
    path = attempt / (stem + '.bam')
    path.symlink_to(Path(sample[alignment]).resolve())
    source_index = Path(sample[index])
    # Some CNVkit versions regenerate an index which is older than the BAM.
    # Copy the small index into the run directory so a native write can never
    # reach the user's original. Support native BAI inputs; CSI-only inputs
    # are indexed by CNVkit in this private directory if needed.
    suffix = '.csi' if source_index.suffix.lower() == '.csi' else '.bai'
    dest = Path(str(path) + suffix)
    shutil.copy2(source_index, dest)
    return str(path)


def _environment_key(env):
    chosen = {
        k: v
        for k, v in env.items()
        if k.startswith(
            ('LD_', 'DYLD_', 'R_', 'OMP_', 'OPENBLAS_', 'MKL_', 'PYTHON', 'HTS_', 'NUMEXPR_', 'LC_')
        )
        or k in ('LANG', 'TZ')
    }
    return signature(chosen)


def _depth_summary(attempt):
    counts = {}
    path = attempt / 'bedcov_events.jsonl'
    if path.is_file():
        for line in path.read_text().splitlines():
            event = json.loads(line)
            row = counts.setdefault(event.get('alignment', 'unspecified'), {'hits': 0, 'misses': 0})
            row['hits' if event['hit'] else 'misses'] += 1
    return counts


def _checkpoint_cnr(directory, key):
    try:
        metadata = json.loads((directory / 'batch_complete.json').read_text())
        path = Path(metadata['cnr'])
        if (
            metadata['key'] == key
            and path.resolve().is_relative_to(directory.resolve())
            and not path.is_symlink()
            and path.is_file()
            and _hash(path) == metadata['sha256']
        ):
            return path
    except (KeyError, OSError, ValueError, TypeError):
        pass
    return None


def _verified_receipt(directory, key):
    if directory.is_symlink():
        return None
    path = directory / 'complete.json'
    if not path.is_file() or path.is_symlink():
        return None
    try:
        record = json.loads(path.read_text())
        cns = Path(record['cns'])
        if (
            signature(record['recipe']) == key
            and cns.resolve().is_relative_to(directory.resolve())
            and cns.is_file()
            and not cns.is_symlink()
            and _hash(cns) == record['sha256']
        ):
            return record
    except (OSError, KeyError, ValueError, TypeError):
        pass
    return None


def _probe_stamp(identity, python, script):
    files = {Path(python), Path(script), Path(identity['Rscript']['path'])}
    for mod in identity['modules'].values():
        root = Path(mod['root'])
        files.add(root)
        files.update(root.rglob('*.py'))
        files.update(root.parent.glob('*.dist-info/RECORD'))
        files.update(root.parent.glob('*.dist-info/METADATA'))
    return signature({str(p): _stamp(p) for p in sorted(files)})


def _verified_probe(runtime, env, config, sample, out):
    python, script = runtime['cnvkit_python'], runtime['cnvkit_script']
    key = signature([python, script, runtime['cnvkit_version'], _environment_key(env), env['PATH']])
    probes = config.setdefault('_global_cnv_probes', {})
    if key in probes:
        identity, stamp, probe = probes[key]
        if _probe_stamp(identity, python, script) == stamp:
            return identity, probe
    probe = Path(tempfile.mkdtemp(prefix='probe-', dir=out))
    policy = probe / 'probe_policy.json'
    atomic_json(policy, dict(alignments={sample['alignment']: sample['index']}))
    identity = _probe_cnvkit(runtime, policy, probe, env)
    if identity['Rscript'].get('DNAcopy_version', {}).get('returncode') != 0:
        raise LocalFailure('unsupported', 'version_check', 'CBS requires working Rscript and DNAcopy')
    probes[key] = identity, _probe_stamp(identity, python, script), probe
    return identity, probe


def _fresh(sample, runtime, config, out, normal, processes, timeout, inputs, archived, rejected):
    # Probe in the native interpreter even when an archive only supplies CNR.
    if not runtime.get('cnvkit_python') or not runtime.get('cnvkit_script'):
        details = '; '.join(r['reason'] for r in rejected)
        raise LocalFailure('missing', 'version_check',
                           f"CNVkit {runtime['cnvkit_version']} required for computation; no verified complete archived CNS. {details}")
    python, script = runtime['cnvkit_python'], runtime['cnvkit_script']
    cnv_env = dict(
        runtime,
        env=runtime.get('cnvkit_env', runtime.get('env', {})),
        extra_path=runtime.get('cnvkit_extra_path', runtime.get('extra_path', [])),
    )
    env = _environment(cnv_env, python)
    out.mkdir(parents=True, exist_ok=True)
    identity, probe = _verified_probe(runtime, env, config, sample, out)
    version = identity['version']
    ref = config['reference']
    fasta, build = Path(ref['fasta']), ref['build']
    by_build = Path(ref['aa_data_repo']) / build
    recipe = dict(
        schema=3,
        runner_sha256=_hash(__file__),
        raw_depth_helper_sha256=_hash(Path(__file__).with_name('concordance_global_depth.py')),
        scope='global_wgs',
        cnvkit_identity=identity,
        sample=sample['sample_id'],
        environment=_environment_key(env),
        build=build,
        inputs=inputs,
        processes=processes,
        method='native_batch_then_explicit_cbs',
        reference_mode='matched_normal' if normal else 'suite_reference',
    )
    reference = by_build / (build + '_cnvkit_filtered_ref.cnn')
    sizes = by_build / (build + '_noAlt.fa.fai')
    if normal:
        recipe['reduced_reference_inputs'] = {
            'fasta_sha256': _hash(fasta),
            'noalt_sha256': _hash(sizes),
            'format': 'suite-reduced-v1',
        }
    else:
        recipe['suite_reference_sha256'] = _hash(reference)
    recipe['reference_fai_sha256'] = _hash(str(fasta) + '.fai')
    if archived:
        recipe['archived_cnr'] = archived
    evidence = 'historical_metadata_only' if archived else 'native_verified'
    key = signature(recipe)
    directory = out / 'cnvkit' / key
    directory.mkdir(parents=True, exist_ok=True)
    receipt = directory / 'complete.json'
    with stage_lock(out, 'cnvkit', key):
        current = _verified_receipt(directory, key)
        if current:
            return dict(
                cns=current['cns'],
                cnvkit_version=version,
                reuse='same_recipe_cache',
                provenance=str(receipt),
                key=key,
                evidence_level=evidence,
            )
        prior_root = config.get('reuse_global_cnv_from')
        if prior_root:
            previous = Path(prior_root).resolve() / 'global_cnv_cache' / 'cnvkit' / key
            previous_record = _verified_receipt(previous, key)
            if previous_record:
                attempt = Path(tempfile.mkdtemp(prefix='import-', dir=directory))
                cns = attempt / 'tumor.cns'
                shutil.copyfile(previous_record['cns'], cns)
                if _hash(cns) != previous_record['sha256']:
                    raise ValueError('Prior global CNS changed during verified import')
                previous_receipt = previous / 'complete.json'
                record = dict(
                    recipe=recipe,
                    reuse='prior_same_recipe_cache',
                    evidence_level=evidence,
                    cns=str(cns),
                    sha256=_hash(cns),
                    directory=str(attempt),
                    completed=now(),
                    segments=_check_cns(cns, str(fasta) + '.fai'),
                    prior_receipt=str(previous_receipt),
                    prior_receipt_sha256=_hash(previous_receipt),
                    archive_rejections=rejected,
                )
                atomic_json(receipt, record)
                return dict(
                    cns=str(cns),
                    cnvkit_version=version,
                    reuse=record['reuse'],
                    provenance=str(receipt),
                    key=key,
                    evidence_level=evidence,
                )
        failures = config.setdefault('_global_cnv_failed_recipes', {})
        if key in failures:
            raise LocalFailure('failed', 'global_cnv', failures[key])
        attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=directory))
        record = dict(
            recipe=recipe,
            created=now(),
            archive_rejections=rejected,
            runtime_identity=str(probe / 'cnvkit_identity.json'),
        )
        try:
            cnr, cns = attempt / 'tumor.cnr', attempt / 'tumor.cns'
            env.setdefault('NUMEXPR_MAX_THREADS', str(processes))
            env.setdefault('MPLCONFIGDIR', str(attempt / 'mplcache'))
            cli = [python, script]
            previous_cnr = None if archived else _checkpoint_cnr(directory, key)
            if previous_cnr:
                cnr = previous_cnr
                record['batch_checkpoint'] = str(directory / 'batch_complete.json')
                reuse = 'same_recipe_cnr_checkpoint'
            elif archived:
                source = Path(archived['file'])
                opener = gzip.open if source.suffix == '.gz' else open
                with opener(source, 'rb') as frm, cnr.open('wb') as to:
                    shutil.copyfileobj(frm, to)
                reuse = 'archived_same_version_cnr'
            else:
                tumor = _private_bam(sample, 'alignment', 'index', attempt, 'tumor')
                alignments = {tumor: [inputs['alignment'], inputs['index']]}
                args = ['batch', tumor, '-m', 'wgs', '-p', str(processes), '-d', str(attempt)]
                if normal:
                    reference = _reduced(fasta, sizes, out)
                    baseline = _private_bam(sample, 'normal_alignment', 'normal_index', attempt, 'normal')
                    alignments[baseline] = [inputs['normal_alignment'], inputs['normal_index']]
                    args += ['--fasta', str(reference), '--normal', baseline]
                else:
                    args += ['-r', str(reference)]
                spec = dict(
                    cache=str(out / 'raw_bedcov'),
                    events=str(attempt / 'bedcov_events.jsonl'),
                    environment=_environment_key(env),
                    alignments=alignments,
                    reference=str(reference.resolve()),
                    reference_sha256=_hash(reference),
                )
                atomic_json(attempt / 'bedcov_spec.json', spec)
                wrapper = Path(__file__).with_name('concordance_global_depth.py')
                _native_cmd(
                    [python, wrapper, script, *args],
                    'global_cnv_batch',
                    attempt,
                    dict(env, AS_GLOBAL_DEPTH_SPEC=str(attempt / 'bedcov_spec.json')),
                    timeout,
                )
                if any(_alignment_id(sample[k]) != v for k, v in inputs.items()):
                    raise ValueError('BAM/index changed during global batch')
                atomic_json(
                    directory / 'batch_complete.json',
                    dict(
                        key=key,
                        cnr=str(cnr),
                        sha256=_hash(cnr),
                        completed=now(),
                        bedcov_counts=_depth_summary(attempt),
                    ),
                )
                reuse = 'native_global_batch'
            # The Suite deliberately segments the CNR again after batch. Use
            # precisely that output, never the .call.cns or the batch-only CNS.
            cns.unlink(missing_ok=True)
            _native_cmd(
                [*cli, 'segment', cnr, '-p', str(processes), '-m', 'cbs', '-o', cns],
                'global_cnv_segment',
                attempt,
                env,
                timeout,
            )
            count = _check_cns(cns, str(fasta) + '.fai')
            if any(_alignment_id(sample[k]) != v for k, v in inputs.items()):
                raise ValueError('BAM/index changed during global CNVkit')
            record.update(
                cns=str(cns),
                sha256=_hash(cns),
                segments=count,
                reuse=reuse,
                directory=str(attempt),
                completed=now(),
                evidence_level=evidence,
                bedcov_counts=_depth_summary(attempt),
            )
            atomic_json(receipt, record)
        except (OSError, ValueError, KeyError, TimeoutError, LocalFailure) as exc:
            failures[key] = f'Global recipe failed: {exc}; see {attempt}'
            atomic_json(attempt / 'failed.json', dict(record, reason=str(exc)))
            raise
    return dict(
        cns=str(cns),
        cnvkit_version=version,
        reuse=reuse,
        provenance=str(receipt),
        key=key,
        evidence_level=evidence,
    )


def global_cns(sample, runtime, entry, config, out):
    """Resolve one complete uncalled CNS. No Suite filter or AA is invoked here."""
    try:
        version = load_builtin_map()['versions'][entry['id']]['cnvkit_version']
        if (
            runtime['id'] != entry['id']
            or runtime.get('cnvkit_version', version) != version
            or entry['suite']['commit'] != load_builtin_map()['versions'][entry['id']]['suite_commit']
        ):
            raise LocalFailure(
                'unsupported', 'version_check', 'Runtime differs from pinned Suite CNVkit recipe'
            )
        out = Path(out).resolve()
        build = config['reference']['build']
        if sample.get('reference_build') != build:
            raise ValueError('Sample reference build differs from global CNVkit reference')
        normal = bool(sample.get('normal_alignment')) and build != 'GRCh38_viral'
        keys = ['alignment', 'index'] + (['normal_alignment', 'normal_index'] if normal else [])
        inputs = {k: _alignment_id(sample[k]) for k in keys}
        archived, rejected = _historical_inputs(sample, config, version)
        if archived and archived['kind'] == 'cns':
            try:
                count = _check_cns(archived['file'], str(config['reference']['fasta']) + '.fai')
            except (ValueError, OSError, KeyError) as exc:
                rejected.append(dict(root=archived['file'], reason=str(exc)))
                archived = None
            else:
                key = signature(dict(archive=archived, inputs=inputs, version=version, scope='global_wgs'))
                receipt = out / 'archives' / key / 'complete.json'
                with stage_lock(out, 'archives', key):
                    atomic_json(
                        receipt,
                        dict(
                            archived,
                            inputs=inputs,
                            created=now(),
                            segments=count,
                            archive_rejections=rejected,
                            cnvkit_version=version,
                            evidence_level='historical_metadata_only',
                        ),
                    )
                return dict(
                    cns=archived['file'],
                    cnvkit_version=version,
                    reuse='archived_same_version_cns',
                    provenance=str(receipt),
                    key=key,
                    evidence_level='historical_metadata_only',
                )
        settings = config.get('global_cnv', {})
        processes = settings.get('processes', 1)
        timeout = settings.get('timeout_seconds', 86400)
        if type(processes) is not int or processes < 1 or type(timeout) is not int or timeout < 1:
            raise ValueError('global_cnv.processes and timeout_seconds must be positive integers')
        return _fresh(sample, runtime, config, out, normal, processes, timeout, inputs, archived, rejected)
    except LocalFailure:
        raise
    except (KeyError, OSError, EOFError, TypeError, ValueError, TimeoutError) as exc:
        raise LocalFailure('failed', 'global_cnv', str(exc)) from exc
