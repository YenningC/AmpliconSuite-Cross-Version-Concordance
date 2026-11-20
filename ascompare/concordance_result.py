"""Reader-facing, sample-level comparisons. Formal decisions stay at cohort level."""
from pathlib import Path

from .concordance_global_report import GLOBAL_METRIC
from .reporting import write_tsv


def write_result(out, downstream, upstream, mode, reference, synthetic, matched_samples=(), *, global_mode=False):
    lookup = {(r['sample_id'], r['target_version'], r['comparator_version']): r for r in upstream}
    matched_by_suite = {(r['sample_id'], version): r for r in matched_samples for version in r['target_suite_versions']}
    rows = []
    for item in downstream:
        key = (item['sample_id'], item['target_version'], item['comparator_version'])
        local = lookup.get(key, {}) if mode == 'with_bam' and item['analysis_kind'] == 'aa_ac' else {}
        row = {k: item[k] for k in (
            'sample_id', 'analysis_kind', 'reference_version', 'source_version', 'target_version',
            'comparator_version', 'comparator_roles', 'identity_anchor', 'n_expected_units',
            'n_evaluable_units', 'n_signal_units', 'target_statuses', 'comparator_statuses',
            'raw_aa_interval_jaccard', 'raw_aa_interval_jaccard_n_units',
            'ac_feature_interval_jaccard', 'ac_feature_interval_jaccard_n_units',
            'ac_class_footprint_jaccard', 'full_label_agreement', 'full_label_agreement_count',
            'full_label_compared_count', 'common_label_agreement', 'common_label_agreement_count',
            'common_label_compared_count', 'common_label_scopes', 'label_comparison_reasons', 'raw_aa_comparison_reasons')}
        for k in ('upstream_interval_jaccard', 'upstream_intersection_bp', 'upstream_union_bp',
                  'upstream_retained_count', 'upstream_evaluable_count', 'upstream_expected_seeds',
                  'upstream_retained_ratio', 'upstream_panel_retained_ratio', 'upstream_unevaluable_count',
                  'upstream_count_status', 'upstream_paired_seeds',
                  'upstream_target_cnvkit_versions', 'upstream_comparator_cnvkit_versions'):
            row[k] = local.get(k)
        if global_mode:
            for k in (GLOBAL_METRIC, 'locus_discordance', 'upstream_similarity_scope', 'upstream_pair_positive_loci', 'upstream_extra_loci',
                      'upstream_missed_loci', 'target_native_seed_count', 'comparator_native_seed_count',
                      'comparator_native_seeds_retained_count', 'comparator_native_seeds_retained_ratio',
                      'target_native_seeds_retained_count', 'target_native_seeds_retained_ratio', 'seed_match_fraction'):
                row[k] = local.get(k)
        row['upstream_similarity_status'] = local.get('upstream_similarity_status', 'not_run' if mode == 'with_bam' and item['analysis_kind'] == 'aa_ac' else 'not_applicable')
        matched = matched_by_suite.get((item['sample_id'], item['target_version']), {}) if mode == 'with_bam' and item['analysis_kind'] == 'aa_ac' else {}
        for field in ('source_cnvkit_version', 'target_cnvkit_version', 'interval_jaccard', 'retention_delta',
                      'paired_seed_count', 'expected_seed_count', 'extra_count', 'missed_count', 'baseline_original_calibration_status',
                      'target_retained_count', 'target_evaluable_count', 'baseline_retained_count', 'baseline_evaluable_count'):
            row['local_cnv_' + field] = matched.get(field)
        row['local_cnv_status'] = matched.get('status', 'not_applicable' if not local or mode != 'with_bam' else 'unavailable')
        if global_mode:
            row['local_cnv_status'] = 'not_applicable_global_workflow'
        row['synthetic'] = synthetic
        rows.append(row)
    out = Path(out)
    write_tsv(out / 'result.tsv', rows)
    def percent(value):
        return '—' if value is None else f'{value:.4%}'
    def counted(value, num, den):
        return '—' if value is None else f'{num}/{den} ({percent(value)})'
    lines = ['# 逐样本跨版本一致性结果', '', f'固定面板参考 **{reference}**。模式 `{mode}`。', '']
    if synthetic:
        lines += ['**人工测试数据。**', '']
    if global_mode:
        lines += ['上游对照两版在全基因组产生的最终种子区段。一个有重叠的连续位点可能覆盖多条原生种子。'
                  '位点保留数与原生种子条数分别列出。下游对抽中的共同位点逐一比较，各版本仅用自己的原生种子启动 AA。', '']
    elif mode == 'with_bam':
        lines += ['空种子 BED 记为“无候选，未评估”，不计入检出和一致性的分母。局部复核只覆盖提供的候选区域及其分析窗口。', '',
                  '## 以输入样本原用 CNVkit 为标准的局部对照', '',
                  '标准版与比较版使用相同的肿瘤和参考策略、分析窗口和来源 AA 筛种器。Jaccard 衡量原候选内最终种子区段的重叠。'
                  '保留数按各版单独计数；新增与遗漏只统计成功配对的原候选，方向为比较版相对标准版。', '',
                  '| 样本 | 标准 CNVkit | 比较 CNVkit | 标准保留/可评 | 比较版保留/可评 | 配对/预期 | 区段 Jaccard | 新增/遗漏 | 来源局部复现原结果 | 比较状态 |',
                  '|---|---|---|---:|---:|---:|---:|---:|---|---|']
        status_text = {'ok': '可评价', 'both_empty': '双方区段均空', 'partial': '仅部分可配对',
                       'not_recorded_in_run': '旧运行未测', 'source_baseline_unavailable': '标准不可用',
                       'target_unavailable': '比较版不可用', 'source_cnvkit_version_unknown': '原版本未知',
                       'no_candidates': '无候选', 'candidate_universe_undeclared': '缺少候选清单',
                       'pair_context_mismatch': '计算条件不一致'}
        def counts(a, b):
            return '—' if a is None or b is None else f'{a}/{b}'
        for r in matched_samples:
            label = status_text.get(r['status'], r['status'])
            if r['identity_anchor']:
                label += '（自身）'
            lines.append('| ' + ' | '.join([str(r['sample_id']), r['source_cnvkit_version'] or '未知',
                r['target_cnvkit_version'] or '未知', counts(r['baseline_retained_count'], r['baseline_evaluable_count']),
                counts(r['target_retained_count'], r['target_evaluable_count']), counts(r['paired_seed_count'], r['expected_seed_count']),
                percent(r['interval_jaccard']), counts(r['extra_count'], r['missed_count']),
                {'passed': '通过', 'failed': '失败', 'unverified': '未核实'}[r['baseline_original_calibration_status']], label]) + ' |')
        lines += ['', '[逐种子对照](local_cnv_comparator_seeds.tsv) · [每例保留比例和差值](local_cnv_comparator_samples.tsv) · '
                  '[整组均值与描述性置信区间](local_cnv_comparator_cohort.tsv)', '', '## Suite 配方与下游结果', '']
    lines += ['全局上游的交并比按完整种子覆盖碱基计算。位点保留比的分母是面板所有成功版本发现的位点并集。'
              '原生种子保留数按目标种子并集覆盖对照种子的比例达到 seed_match_fraction 计数。新增、遗漏的位点均相对于表中对照。两版都没有种子时，分母为零的指标留空。' if global_mode else '表中每一行比较同一样本的一对版本。上游区段相似度是原候选范围内两版最终种子区段的交并比（Jaccard）。保留数按目标版本计算，分母为可评价原候选数，另列原候选总数及全候选保留比。全候选保留比要求所有预期候选均可评价。双方均无区段信号时相似度留空。',
              '', '下游分别比较 AA 重建区段、AC 分类区段及标签集合。标签在至少一侧有相应信号的可评价单位中计数；完整标签要求类别支持一致，共同标签仅检查双方支持的类别。两个 Jaccard 均只比较位置。表内百分数是单样本观察值。', '',
              '[精确数值和状态](result.tsv) · [下游逐单位证据](downstream_similarity_units.tsv) · [整组统计和兼容判定](report.md)', '']
    if mode == 'with_bam':
        lines += ['[上游全部版本对及交并长度](upstream_similarity_pairs.tsv) · [上游整组描述性区间](upstream_similarity_cohort.tsv)', '']
    def cell(value):
        return str(value).replace('|', r'\|').replace('\n', ' ')
    groups = {}
    for r in rows:
        groups.setdefault((r['sample_id'], r['analysis_kind'], r['comparator_version']), []).append(r)
    for (sample, kind, comp), group in groups.items():
        role = '、'.join('固定参考' if x == 'fixed_reference' else '输入来源版本' for x in group[0]['comparator_roles'])
        lines += [f'## {cell(sample)} · {kind} · 对照 {cell(comp)}（{role}）', '',
                  '| 目标 Suite | 配方种子区段 | 种子保留数/可评数 | 候选总数 | 全候选保留比 | AA 区段 | AC 区段 | 完整标签一致率 | 共同标签一致率 |',
                  '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        if global_mode:
            lines[-2] = lines[-2].replace('种子保留数/可评数', '有种子位点/可评位点').replace('候选总数', '并集位点数').replace('全候选保留比', '位点保留比') + ' 原生种子数（目标/对照） | 对照原生种子保留数/总数 | 新增/遗漏位点 | 双向较大未保留率 |'
            lines[-1] += '---:|---:|---:|---:|'
        for r in group:
            retained = counted(r['upstream_retained_ratio'], r['upstream_retained_count'], r['upstream_evaluable_count'])
            n = r['upstream_expected_seeds']
            cells = [cell(r['target_version']) + ('（自身）' if r['identity_anchor'] else ''),
                     percent(r['upstream_interval_jaccard']), retained, '—' if n is None else str(n),
                     percent(r['upstream_panel_retained_ratio']), percent(r['raw_aa_interval_jaccard']), percent(r['ac_feature_interval_jaccard']),
                     counted(r['full_label_agreement'], r['full_label_agreement_count'], r['full_label_compared_count']),
                     counted(r['common_label_agreement'], r['common_label_agreement_count'], r['common_label_compared_count'])]
            if global_mode:
                for a, b in [('target_native_seed_count', 'comparator_native_seed_count'),
                             ('comparator_native_seeds_retained_count', 'comparator_native_seed_count'),
                             ('upstream_extra_loci', 'upstream_missed_loci')]:
                    cells.append(' / '.join('—' if r[k] is None else str(r[k]) for k in (a, b)))
                cells.append(percent(r[GLOBAL_METRIC]))
            lines.append('| ' + ' | '.join(cells) + ' |')
        lines.append('')
    lines += ['“—”表示该指标缺少可评价证据、类别范围不符、双方均空或模式不适用。具体原因见精确数值表及逐单位证据。自身比较只用于核对。', '']
    (out / 'result.md').write_text('\n'.join(lines), encoding='utf-8')
    return rows
