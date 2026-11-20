"""Pair native whole-genome seed calls on common loci; keep AA inputs native."""
from bisect import bisect_right
from collections import defaultdict

from .runner import empty_record
from .statistics import merge_intervals

GLOBAL_SCOPE = 'whole_genome_native_seed_discovery_and_sampled_locus_aa_ac'
GLOBAL_METRIC = 'global_native_seed_max_directional_loss'
LOCUS_RULE = 'connected_components_of_positive_length_native_seed_overlap'


def _interval(seed):
    chrom, start, end = (seed[key] for key in ('chrom', 'start', 'end'))
    if (not isinstance(chrom, str) or not chrom or type(start) is not int or type(end) is not int
            or not 0 <= start < end):
        raise ValueError(f'Invalid native seed interval: {seed}')
    return chrom, start, end


def global_loci(seed_sets):
    """Strict overlap only. Adjacent half-open intervals are separate comparison units."""
    intervals = sorted({_interval(s) for call in seed_sets.values() if call['status'] == 'ok'
                        for s in call['seeds']})
    components = []
    for chrom, start, end in intervals:
        if components and chrom == components[-1][0] and start < components[-1][2]:
            components[-1][2] = max(components[-1][2], end)
        else:
            components.append([chrom, start, end])
    return [dict(chrom=c, start=a, end=b, seed_id=f'{c}:{a}-{b}') for c, a, b in components]


def seeds_by_locus(loci, seeds):
    """Assign each native seed, without changing its boundaries, to its containing locus."""
    starts = [(l['chrom'], l['start']) for l in loci]
    members = defaultdict(list)
    for seed in seeds:
        c, a, b = _interval(seed)
        idx = bisect_right(starts, (c, a)) - 1
        if idx < 0 or loci[idx]['chrom'] != c or b > loci[idx]['end']:
            raise ValueError('Native seed does not belong to a declared global locus')
        members[loci[idx]['seed_id']].append(seed)
    return members


def global_report_units(sample, joints, versions):
    """AA output is projected onto every sampled locus, including loci absent in its own seed BED."""
    from .concordance import project_record

    calls = sample['global_seed_sets']
    loci = global_loci(calls)
    ids = {l['seed_id'] for l in loci}
    if set(sample['upstream_seed_ids']) != ids or not set(sample['selected_seed_ids']) <= ids:
        raise ValueError('Global report locus universe differs from the declared sampling panel')
    selected = set(sample['selected_seed_ids'])
    sid, source = sample['sample_id'], sample['source_version']
    upstream, downstream = [], []
    for vid in versions:
        call = calls.get(vid, {'status': 'missing', 'stage': 'seed_processing', 'reason': 'No global seed record'})
        ok = call['status'] == 'ok'
        stage = call.get('stage', 'seed_processing')
        reason = call.get('reason', '')
        members = seeds_by_locus(loci, call['seeds']) if ok else {}
        joint = joints.get(vid)
        selected_native_ids = [s['seed_id'] for key in selected for s in members.get(key, [])]
        executed = sample.get('global_executed_seed_ids')
        if ok and executed is not None and sorted(executed.get(vid, ())) != sorted(selected_native_ids):
            raise ValueError(f'Native seed selection disagrees with the sampling record for {vid}')
        has_native_selection = bool(selected_native_ids)
        for locus in loci or [None]:
            key = locus['seed_id'] if locus else '__no_seeds__'
            seeds = members.get(key, [])
            intervals = merge_intervals([_interval(s) for s in seeds])
            status = 'retained' if seeds else 'not_retained' if locus else 'no_seeds'
            context = dict(seed=locus, seed_id=key, sample_id=sid, source_version=source, version=vid,
                           global_scope=True, native_seed_count=len(seeds) if ok else None,
                           native_seeds=seeds if ok else None)
            upstream.append(dict(context, cnvkit_version=call.get('cnvkit_version'),
                native_recipe_match=call.get('native_recipe_match'), stage=stage,
                reason=reason, matching_rule=LOCUS_RULE,
                reference_mode=call.get('reference_mode'), status=status if ok else call['status'],
                retained_intervals=intervals if ok else None,
                retention_fraction=(sum(b-a for _, a, b in intervals)/(locus['end']-locus['start'])
                                    if ok and locus else None),
                calibration_role='source_control' if vid == source else 'candidate_observation',
                calibration_status=sample.get('global_source_seed_calibration', 'unverified') if vid == source else 'unverified',
                comparison_eligible=bool(ok and locus)))
            if locus and key not in selected:
                continue
            if not locus or not ok:
                down = empty_record(sid, key, vid, 'no_seeds' if ok else call['status'], stage,
                                    'Completed global discovery found no native seeds' if ok else reason,
                                    comparison_eligible=False, detected=None)
            elif joint:
                down = project_record(joint, locus)
            elif has_native_selection:
                down = empty_record(sid, key, vid, 'missing', 'aa', 'Missing selected-seed reconstruction',
                                    comparison_eligible=False, detected=None)
            else:
                down = empty_record(sid, key, vid, 'not_detected', 'seed_processing',
                                    'This version contributed no native seeds to any sampled locus',
                                    comparison_eligible=True, detected=False, supported_labels=call.get('supported_labels'),
                                    aa_executed=False, ac_executed=False, aa_reconstruction_performed=False,
                                    provenance_kind='completed_empty_native_selection')
            down.update(context, analysis_kind='aa_ac')
            downstream.append(down)
    return upstream, downstream
