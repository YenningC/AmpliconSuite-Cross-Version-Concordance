# 全局 CNVkit 与各版本原生种子

从 0.5.0 开发版开始，新 `with_bam` 清单使用 `"cnv_scope": "global"`。推荐把范围写入清单。既有清单带有 `local_cnv` 且未填写 `cnv_scope` 时仍执行原局部设计。升级该清单时删除 `local_cnv`、添加 `"cnv_scope": "global"` 并选择新输出目录。原自定义匹配比例可写入顶层 `seed_match_fraction`（默认 1.0）。配对抽样上限、随机种子、AA 读取预算、统计容忍度和固定参考 Suite 1.6.2 沿用清单设置。

## 一例样本怎样运行

1. 提供排序肿瘤 BAM、现成索引、匹配的 FASTA/FAI 及 AA_DATA_REPO。配对设计另提供正常 BAM 和索引。无正常对照时，使用该 build 的 Suite 原生 CNVkit 参考资源。原始 Suite 输出根目录填入 `samples[].source_result.result_dir`，程序查找运行记录、CNVkit 文件和原生种子。原种子 BED 对全局发现可选；原结果缺少可核对证据时，来源复现显示未核实。
2. 每个 Suite 按锁定配方取得未经 call 的全局 CNS。配方使用 0.9.10 或 0.9.13。程序读取历史记录中的实际 CNVkit 版本。来源为其他版本时，保留该来源身份，分别取得目标版本结果。不同样本的肿瘤深度和匹配正常对照由输入身份区分。常规 WGS 有匹配正常 BAM 时用配对策略；GRCh38_viral 按 Suite 的肿瘤参考策略处理。
3. 各 Suite 调用各自锁定的 CNS 转换、染色体臂预筛和 AA 的 amplified_intervals 筛种。所有版本完成发现或留下失败记录后，把成功版本的种子按正长度重叠连成比较位点；端点仅相接的区间保持独立。每例按原 `sampling.random_seed` 对位点进行一次配对抽样，默认最多 10 个，`sampling.max_seeds: "all"` 可全选。每版 AA 仅接收该版在抽中位点的原生种子，保留坐标、顺序和 BED 字段。同一位点可以有多个原生种子，实际 AA 输入行数可大于抽中的位点数。AA 仍能向相连区域扩展。
4. 各版 AC 在相应的 AA 图上分类。种子额外检出、遗漏和上下游差异分别报告。上游按全体已观察位点统计，下游按抽中位点统计。全局成功但零种子、部分版本失败、缺原始文件、AA 或 AC 失败分别记录。种子本身只是启动重建的候选。

## 文件复用

| 已有材料 | 处理 |
|---|---|
| 该版 Suite 原始结果根目录，含原生运行记录与完整原始 `.cns` | 检查版本、全局 WGS/CBS 命令、肿瘤/正常文件、参考策略、输出命名与数据内容。核对通过直接使用；记录历史证据等级。|
| 同版本原生全局运行仅保留 `.cnr` 或 `.cnr.gz` | 执行同版本的 CBS 分段。压缩输入仅在本次工作目录内展开。|
| 同一次分析有多个 Suite 需要同一 CNV 配方，或同目录续跑 | 核对输入、版本、参考、参数、环境和完成文件摘要。匹配即共用全局 CNS。|
| 不同 CNVkit 配方发出完全相同的原生覆盖度查询 | 在 BAM/索引、底层计数程序、BED 和计数选项一致时，共用原始计数输出；每版自行处理覆盖度、构建参考、校正和分段。|
| 只有其他版本的 `.cnr`、`.cns`、来源不明的 `.cnn`、局部结果或 `.call.cns` | 保留为历史资料；需要目标版本全局结果时，从可核实的步骤计算。源码、输入或环境不明时不能直接证明两个版本的全流程可复用。|

同一样本的其他历史 Suite 结果可以列在 `samples[].global_cnv_results`，每项填写 `{"result_dir": "/data/another_suite_result"}`。另有 `comparison_results[].result_dir` 时也会检查。可在归档项填写已有的原始 `cns_sha256`。若历史 CNS 与配套 CNR 的染色体覆盖端点相悖，程序优先从同版 CNR 重新分段。若已有 Concordance 全局运行，顶层 `reuse_global_cnv_from` 可指向其分析输出根目录，从只读目录复用经校验的同配方 CNS。有历史 CNS 而缺少运行记录，或者记录中的 BAM 已无法对应当前数据时，程序保留拒绝原因。历史记录常缺生成时的依赖和参考资源摘要，观察值与正式兼容结论的证据等级分开标记。归档需保留原生目录布局与能识别的 BAM/正常对照路径。

## 执行与验收

填写本机对应的组件运行环境后，在本地运行以下命令。

```sh
asconcordance check --manifest concordance.json --out runs/global
asconcordance run --manifest concordance.json --out runs/global
asconcordance report --from-run runs/global --out reports/global
```

`global_cnv.processes` 默认 1，可单独设置 CNVkit 工作进程数。`global_cnv.timeout_seconds` 默认每个全局命令 86400 秒。AA/AC 的限时与读取预算仍由原字段控制。直接重新运行相同命令可恢复完成的有效阶段。全局 CNV 全流程的实际运行时间以本机日志为准。

| 输出 | 怎样读 |
|---|---|
| `global_upstream/<sample>/discovery.json` | 各 Suite 的 CNVkit 版本、CNS 路径及复用来源、筛种是否完成、原生种子和失败原因 |
| `sampling/<sample>.json` | 全部比较位点、抽中的位点及各版实际交给 AA 的种子 ID |
| `source_bridge/<sample>.seeds.json` | 来源真实 CNVkit 与面板配方、原种子数和重现状态 |
| `stage_timing.tsv`、`global_cnv_cache/` | 区分全局计算、筛种、AA、AC、历史文件或缓存命中；原生中间文件和命令日志保存在本地 |
| `result.md`、`result.tsv`、`sample_reports.md`、`cohort_compatibility_metrics.tsv` | 各样本观察差异、逐样本图以及按患者聚合的总体比较 |

[人工 WGS 验证结果与对照](validation/global_wgs.synthetic.json)

上游正式指标分别计算两个方向的原生种子未保留比例，取两者中较大者。只有一版有种子时记为 1，双方都没有时不可估计。两个版本均无种子时不把空集计作 100% 一致。区段碱基 Jaccard 能显示同一位点内边界的变化。种子数保留比按 `seed_match_fraction` 判断对照原生种子被目标种子并集覆盖的比例。兼容容限 `statistics.upstream_loss_margin` 用于该版本对较差方向的原生种子未保留比例；升级研究时应确认这个指标对应的可接受差异。AA/AC 的下游指标、患者聚合及多重比较沿用既有方法。由于传入 AA 的种子允许随版本变化，下游差异同时包含上游选种与 AA、AC 的影响。

[CNVkit 官方流程](https://cnvkit.readthedocs.io/en/stable/pipeline.html) · [固定参考 Suite 源码](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/9b6cf5696ddc/AmpliconSuite-pipeline.py)
