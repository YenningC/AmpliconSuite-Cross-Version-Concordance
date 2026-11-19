# Batch samples and cohort reports

版本 0.4.0 增加 `import-samples`，把一张本地 TSV 样本表与共享 JSON 配置合成一个清单。后续仍使用同一个 `check` / `run` 入口，按样本首次出现的顺序逐例分析，再汇总整组。样本编号、输出文件夹名由使用者指定，无需按固定格式重新命名样本。

## 1. 准备样本表

选择 `examples/samples.with_bam.tsv` 或 `examples/samples.results_only.tsv`，将人工占位路径替换为真实本地路径。文件是 UTF-8 的制表符分隔文本，支持表格软件导出的 UTF-8 BOM。示例中的 SYNTHETIC 路径不存在，不能直接执行科学分析。

| 列 | 要求与含义 |
|---|---|
| `sample_id` | 每行必填。同一份样本跨版本使用同一个编号 |
| `source_version` | 每行必填。原输入来源，例如 `suite-1.6.1` |
| `reference_build` | 每行必填，例如 `hg19`，必须与实际结果及配置匹配 |
| `result_version` | 该行现有结果的版本；留空即 `source_version` |
| `result_dir` / `run_dir` / `aa_dir` + `ac_dir` | 每行明确一份现有结果目录。Suite 顶层目录用 `result_dir`；含 `aa`、`ac` 子目录的规范目录用 `run_dir`；也可分别指定两个目录。存在歧义时报错 |
| `patient_id` | 强烈建议填写。同一患者的多个样本使用相同 ID；使用时所有样本都必须有值 |
| `alignment`, `index`, `seeds_bed` | `with_bam` 必需。排序 BAM、现成索引和源结果行的最终候选种子 BED |
| `normal_alignment`, `normal_index` | 成对填写。`matched_normal` 需要同一患者的正常对照 |
| `aa_version`, `ac_version` | 该行结果原始组件版本声明；不会被当作历史环境已核验的证据 |
| `source_cnvkit_version` | 原样本实际使用的精确 CNVkit 版本，如 `0.9.14`；用于该样本的局部标准。来源结果有单份记录 CNVkit 执行命令与版本的 `*_run_metadata.json` 时可留空；声明和原记录冲突时标为不可用 |
| `baseline_cns`, `source_seed_bed`, `baseline_provenance` | 原候选再保留的历史校准证据；缺失时保留未核实状态，仍可描述已实测的局部版本差异 |
| `stages_json` | 可选的规范化运行阶段证据文件 |

每行是一份现有版本结果。`results_only` 中，同一个 `sample_id` 可以出现多行，分别填写不同 `result_version`，归并到该样本的 `comparison_results`。每例必须恰好有一行来源版本结果。没有 BAM 时不能生成未提供的 AA 重建。`with_bam` 通常每例一行来源结果，程序将对内置目标版本运行局部上下游流程。

同一样本各行的 `sample_id`、`source_version`、`reference_build` 每行重复填写，其他样本级字段可以只填一次；重复填写的非空值必须相同。重复版本、未知列、缺列、字段冲突会被拒绝。样本和患者编号使用字母、数字、下划线、短横线或点，避免空格及路径分隔符。

## 2. 准备共享配置

选择 `examples/batch.with_bam.json` 或 `examples/batch.results_only.json`。配置不带 `samples`，或设置为空列表；样本全部来自 TSV。参考、内置目标版本、实际运行环境、计算预算和统计参数沿用单样本 JSON 配置。

`with_bam` 示例的 `runtimes` 为空，使用前须配置实际精确版本环境。参考与正常参考策略必须适合整个批次。共享 `background_bed` 应事先验证适合这些样本及测序条件；程序逐例去掉与候选及延伸窗口重叠的部分，因此每例实际有效背景可能不同。不适用同一背景或实验条件的样本应分批运行。来源版本不同的样本自动分层，不将跨来源结果混在一个统计结论内。

共享配置可加入顶层 `cnvkit_runtimes`，逐版填写 `version`、`python` 和 `script`。为整批所有样本实际使用的 CNVkit 准备精确环境。可用 `local_cnv.comparison_versions` 追加不在所选 Suite 配方中的目标。各例从自身来源 CNVkit 起步，并对所有目标固定该例来源 Suite 的 AA 筛种器。[完整示例](local_cnv_calibration.md)

学习阶段可用 `statistics.primary_pairs: []` 只报告观察差异与区间。验证阶段先锁定版本对、容忍阈值和分析计划，再导入独立验证数据。默认不提供保证通过的容忍阈值；详见 [队列设计](cohort_study.md)。

## 3. 导入、检查、运行

```sh
asconcordance import-samples --samples samples.tsv --config batch.json --out batch.manifest.json
asconcordance check --manifest batch.manifest.json --out runs/batch01
asconcordance run --manifest batch.manifest.json --out runs/batch01
```

TSV 内的相对文件路径以 TSV 所在目录为基准；共享配置中的路径以 JSON 所在目录为基准。生成的清单使用绝对路径，并记录两个输入文件的 SHA-256 和样本顺序。导入不会打开 BAM 或启动 AA、AC、CNVkit；后续检查和运行才执行相应核验。导入拒绝覆盖已有清单。

进度显示 `[当前样本/总样本]`。`sample_progress.json` 和 `sample_progress.tsv` 保存顺序、最近状态、时间及已产生结果的状态计数。`finished` 表示该例处理结束，具体版本是否成功仍须看状态计数与异常记录。单例输入或已捕获的组件失败不阻止其他可处理样本；进程中断或未处理异常保留已完成缓存。重复同一个 `run` 命令继续，不能同时向同一输出目录运行两个任务。

## 4. 查看结果

| 输出 | 用途 |
|---|---|
| `report.md` | 整组首页，配对检出变化 Δ、区间、计数与结论 |
| `sample_reports.md`, `sample_plot_index.tsv` | 全部样本图链接，每例每种分析类型一张 |
| `sample_ladders/*.{png,pdf,svg}` | 各版本检出、标签、区段与类别区段相似度；标出额外/遗漏/分类改变和有效分母。`with_bam` 增加已知候选保留面板 |
| `sample_similarity.tsv`, `upstream_samples.tsv` | 每例观察数值和状态 |
| `local_cnv_comparator_{seeds,samples,cohort}.tsv` | 以每例原 CNVkit 的局部结果为标准，报告配对保留数、额外/遗漏、保留率差、区段 Jaccard 和整组描述性置信区间 |
| `local_cnv_comparator_cohort.{png,pdf,svg}` | 有 BAM 时的来源 CNVkit 分层队列对照图 |
| `similarity_ladder.{png,pdf,svg}` | 来源分层的整组相似性及保守区间 |
| `cohort_compatibility.{png,pdf,svg}` | 整组四项下游失配及适用的上游候选丢失、区间、阈值、独立单位数、有效/应评价数和联合判定 |
| `cohort_compatibility_metrics.tsv` | 整组统计兼容图的精确绘图数据，包含区间类型、校正后 alpha、阈值、可评价状态和判定理由 |
| `compatibility_table.*`, `compatibility_combined.tsv` | 按来源与目标版本分层的三态兼容汇总 |
| `downstream/` | 固定参考及其他版本之间的配对统计明细 |

逐样本图只给观察相似度，不能从一个人的多个种子推断总体等效。整组统计先在每例内平均种子，再在同患者内平均样本，按独立患者等权；未填患者编号时假定样本独立。不同来源 Suite 和 AA＋AC / 仅 AC 分层。

整组兼容图的主要比较使用与原判定完全相同、经过误差预算和多重比较校正的单侧界，实心点和实线表示这组判定界。其下界与上界分别为单侧界，不是普通双侧置信区间。探索性比较显示描述性双侧区间和空心点，不转成正式兼容结论。图中界值线表示分析前定义的最大允许失配。单项失配在界内也仍需满足完整性、类别支持、上游校准和来源桥接等联合条件。

没有信号、失败、缺文件和类别不支持保留为明确状态，不补为零失配。区间宽或阈值未设时输出证据不足。相似度排名不是真实性或完整 Suite 可互换概率。固定参考始终在图标题中标明，目前为 `suite-1.6.2`。

## English summary

Use `import-samples --samples samples.tsv --config batch.json --out batch.manifest.json`, then the ordinary `check` and `run` commands. Each TSV row describes one existing version result. Required columns are `sample_id`, `source_version`, `reference_build`, and an explicit result directory specification. `result_version` defaults to the source version. Repeated sample IDs merge distinct existing versions; each sample needs exactly one source result. First appearance defines serial execution order. Use `patient_id` consistently for linked samples. With BAM, also provide alignment, index and source-seed BED; matched-normal analysis needs normal alignment and index.

Paths in the TSV and shared JSON are resolved against their respective parent directories. Templates contain invented placeholder paths. Configure actual runtimes before with-BAM execution. Import reads no BAM and starts no algorithms. Shared background and reference settings must suit the cohort; otherwise use separate batches.

For each with-BAM sample, record `source_cnvkit_version` or supply the original Suite run metadata with its CNVkit command and version. An optional shared `cnvkit_runtimes` registry maps exact versions to executable paths. `local_cnv.comparison_versions` can add further targets. Each patient’s fixed-source-filter CNVkit comparison is summarized in `local_cnv_comparator_samples.tsv`.

Open `sample_reports.md` for every sample and `report.md` for cohort results. `cohort_compatibility_metrics.tsv` is the exact figure data. Primary endpoint limits are the existing multiplicity-adjusted one-sided decision limits; exploratory results use descriptive intervals only. Patients, when declared, are the independent units. Missing margins, inadequate samples and failed scientific prerequisites remain insufficient evidence. No patient data belongs in the repository.
