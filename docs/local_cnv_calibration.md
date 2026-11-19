# 局部 CNVkit 版本对照与原候选再保留校准

局部分析衡量两个量。其一是其他 CNVkit 与该样本来源版本在相同局部条件下的结果相似性；其二是来源版本在有限窗口重新计算后，能再保留多少原始候选。默认窗口为候选两侧各 100 kb，另加用户提供的可信背景 BED；所有窗口合并后的长度默认最多 25 Mb。所有 BAM 均要求已有索引；不会扩大为全基因组扫描。

`local_cnv.reference_mode` 默认 `matched_normal`，要求配对正常 BAM。只有肿瘤 BAM 时，可以**明确**配置 `flat`：程序使用原生 `cnvkit reference -t <有限 targets> -f <FASTA>` 构建平坦参考，不读取正常样本，也不静默切换参考策略。两种策略仍都要求可靠背景 BED。`flat` 不是正常样本，不能声称完全复现原来的正常参考/归一化；未经同策略校准时保持 `unverified`。策略写入每条结果、provenance 和缓存身份，切换策略会重新计算。

CNVkit 原生 coverage、reference、fix、CBS segment 计算保留；WGS 风格禁用 edge correction。coverage 中原先从 BAM 开头取样的排序/读长辅助步骤，改为排序头检查和有限 BED 内最多 1,000 条 reads 的读长取样，明确记录 `bounded_indexed_auxiliary_probes`。不修改安装的上游源码。CBS 跨未测窗口间隙的输出仅裁回实际测量窗口，避免把未测区段当作 CNV 证据。所有局部 CN 区段，包括中性区段，都传入该版 AA 的 `amplified_intervals.py`；该步骤保留其必要的有限、索引式背景估计，并受读取预算限制。

局部中位数归一化、背景和分段上下文与全基因组不同。**候选保留率可以描述，但未经校准不能当作整条 Suite 流程兼容证据。** `min_seed_overlap` 默认 1.0（完全覆盖），是可配置的基准匹配阈值，不是 CNVkit/AA 的官方阈值。它计算合并后的保留区段与原种子的交集长度/原种子长度，重叠不会重复计数。

`builtin_cnvkit_map.json` 来自锁定 Suite 源码中的容器构建配方。1.5.2 及以前的七个版本对应 CNVkit 0.9.10，1.5.3 和三个 1.6.x 版本对应 0.9.13。程序核查实际 CNVkit 版本和导入源码。0.9.14 等受控版本记录为 `controlled_cnvkit`，报告保留其观察结果和实验范围。

## 0.4.2 的同版本局部标准

对于每个样本，程序首先识别原 CNVkit 版本。原生 Suite 结果可通过 `source_result.result_dir` 或 `run_dir` 下单份 `*_run_metadata.json` 提供。程序在该文件记录有效 `cnvkit_cmd` 时读取精确的 `cnvkit_version`，并核对其中已记录的 Suite 版本及 build。也可以用样本级 `source_cnvkit_version` 声明原版本。两者冲突、多份元数据有歧义或版本未知时，保留缺失原因。`source_cnvkit_evidence` 记录版本依据。Suite 的构建配方可用于规划目标面板；历史输入使用的实际 CNVkit 由原始记录或样本声明指定。

例如样本来自 Suite 1.6.1，原运行使用 CNVkit 0.9.14。局部标准为 0.9.14。若所选 Suite 涵盖 0.9.10、0.9.13，则这两版自动成为目标。若面板仅含配方 0.9.13，仍可在 `comparison_versions` 显式加入 0.9.10。版本不同的样本各用自己的标准。下列字段合并到已有的单样本或批量共享 JSON，路径替换为本机已安装的精确版本。示例中的 0.9.13 由既有 Suite 运行环境提供。

```json
{
  "cnvkit_runtimes": [
    {"version": "0.9.14", "python": "/env/cnv014/bin/python", "script": "/env/cnv014/bin/cnvkit.py"},
    {"version": "0.9.10", "python": "/env/cnv010/bin/python", "script": "/env/cnv010/bin/cnvkit.py"}
  ],
  "local_cnv": {
    "reference_mode": "matched_normal",
    "background_bed": "/resources/verified_background.bed",
    "comparison_versions": ["0.9.10"]
  }
}
```

- **登记环境。** `cnvkit_runtimes` 是按精确版本索引的环境列表，每版一条，可提供 `extra_path` 和字符串映射 `env`。程序优先采用该版的显式登记，也可复用已配置的相应 Suite 环境。存在多个不同可选环境时，显式登记可以消除歧义。实际运行还要核对探针版本和组件身份。
- **选择目标。** `local_cnv.comparison_versions` 是可选的额外精确版本号列表，与原标准和所选 Suite 对应的 CNVkit 去重合并；仅登记环境本身不会将它加入目标。某版对应多套 Suite 时，`target_suite_versions` 列出这些 Suite，保留一份局部对照。单独指定的版本可没有对应的 Suite 标签。
- **控制计算条件。** 用同一样本的原始候选、邻近窗口、背景、肿瘤/正常输入、参考策略和筛种条件。每版均应用来源 Suite 的同一 AA 筛种器。输入、窗口与筛种器记录形成局部设计标识。存在已核实且条件相符的来源原生作业时可复用。没有精确环境时保留缺失状态。
- **阅读结果。** `local_cnv_comparator_seeds.tsv` 可以核查每个候选双方保留状态、覆盖比例、实测坐标及原因。`local_cnv_comparator_samples.tsv` 的 `baseline_retained_count/baseline_evaluable_count` 和 `target_retained_count/target_evaluable_count` 分别描述相对原候选的保留数。对双方可配对候选计算 `extra_count` 和 `missed_count`，方向是目标相对标准。完整面板的 `retention_delta` 为（额外数 − 遗漏数）/ 可配对数。`interval_jaccard` 为原候选内两个最终区段并集的交集长度/并集长度，多个候选重叠的碱基只计一次。双方都空时 Jaccard 留空；只有部分候选能配对时单列 `observed_*`。自身版本比较标为 `identity_anchor`。
- **整组结果。** `local_cnv_comparator_cohort.tsv` 按来源 Suite、来源 CNVkit、目标 CNVkit、build 和正常参考策略分层，给出逐样本/逐患者等权均值、独立单位数及描述性双侧区间。单例观察和这一组探索性指标用于解释版本差异。正式 Suite 兼容判定仍使用预设端点、重复性和原始来源校准门槛。

`results_only` 显示 `not_applicable_results_only`。旧运行缺少这套固定筛种器与局部标准的记录，重绘时显示 `not_recorded_in_run`。获取该终点需新增一次受限 `with_bam` 运行。

<a id="synthetic-control-fixture"></a>

## 生成人工配对 BAM 并接入自己的环境

从源码根目录使用已安装 `pysam` 的 Python，提供一个尚不存在的目录。脚本随项目按 MIT 许可发布。它只写测试输入，耗时主要是生成两份 BAM 及索引。[生成脚本](../scripts/prepare_synthetic_cnv_control.py) · [许可](../LICENSE)

```sh
python scripts/prepare_synthetic_cnv_control.py --out ./SYNTHETIC_CNV_control
```

输出目录中的 `fixture.json` 是带有 `SYNTHETIC_ONLY` 标签的文件清单。所有路径都相对于该清单所在目录。两条 15 Mb 的人工染色体采用重复 ACGT 序列，`hg19` 仅作为 AA 识别的目录标签。肿瘤和配对正常 BAM 分别含 154,314 和 149,980 条人工比对记录，肿瘤在 chr1 的一个局部区段增加了覆盖。清单还列出候选 BED（chr1:5,100,000–5,210,000，起点从 0 计、终点不含）、两个背景区段、FASTA 和 AA 资源目录。

要把这些输入接入现有的 `with_bam` 测试配置，请用清单中的文件替换下表字段。把路径转换为绝对路径，或改写成相对于**运行清单**的路径。保留原有的 Suite 运行环境和与该人工数据对应的 `source_result` 配置；生成器只提供上游测试输入，完整 CLI 另需来源 AA/AC 结果。给运行清单加上顶层 `"synthetic": true`。

| `fixture.json` 字段 | `with_bam` 运行清单字段 |
|---|---|
| `tumor_bam`、`tumor_index` | 样本的 `alignment`、`index` |
| `normal_bam`、`normal_index`、`seeds_bed` | 样本的 `normal_alignment`、`normal_index`、`seeds_bed` |
| `build_label`、`fasta`、`aa_data_repo` | `reference.build`、`reference.fasta`、`reference.aa_data_repo`；样本的 `reference_build` 同为 `hg19` |
| `background_bed` | `local_cnv.background_bed` |

本次人工对照将 `suite-1.6.1` 设为来源 Suite、`0.9.14` 声明为人工样本的来源 CNVkit；Suite 1.6.1 和固定参考 Suite 1.6.2 的原生配方均为 0.9.13。按上文环境登记示例提供 0.9.14、0.9.13、0.9.10 的实际解释器、入口及各自的 R/DNAcopy 依赖，并将 0.9.10 加入 `comparison_versions`。局部条件为 `reference_mode=matched_normal`、`flank_bp=100000`、`bin_size=1000`、`max_window_bp=2000000`、`min_seed_overlap=1.0`、`cngain=4.5`、`cnsize_min=50000`。该单候选实验适合检查实际组件调用、区段比较和来源版本对齐；评估临床样本条件仍需相应样本及预先规定的判定指标。

## 原始候选何时能通过校准

需要保留原始运行时的输入和输出身份记录，以及原始 CNS、最终候选 BED。仅声明 Suite/CNVkit 版本、设置 `calibrated=true`，或者事后给现有文件计算指纹，都不会构成原运行的可追溯证据。来源的 `reference_mode` 也要一致；来源未记录时仅允许按带有正常 BAM 身份证据的 `matched_normal` 处理，不能据此校准 `flat`。

样本配置可提供以下字段（均可使用本地文件）：

```json
{
  "source_version": "suite-1.3.4",
  "source_cnvkit_version": "0.9.10",
  "baseline_cns": "/data/original_run/sample.cns",
  "source_seed_bed": "/data/original_run/sample_AA_CNV_SEEDS.bed",
  "baseline_provenance": "/data/original_run/upstream_provenance.json"
}
```

`baseline_provenance` 要求以下结构。下面是字段示意，**不是可以直接填入当前文件指纹以追认历史运行的模板**：

```json
{
  "provenance_kind": "original_upstream_run",
  "cnvkit_version": "0.9.10",
  "reference_build": "hg19",
  "reference_mode": "matched_normal",
  "inputs": {
    "alignment": {
      "path": "/data/tumor.bam",
      "bytes": 123456789,
      "mtime_ns": 1234567890000000000,
      "method": "sampled_first_middle_last_64KiB",
      "blocks": [
        {"offset": 0, "sha256": "original-first-block-digest"},
        {"offset": 61695626, "sha256": "original-middle-block-digest"},
        {"offset": 123391253, "sha256": "original-last-block-digest"}
      ]
    },
    "normal_alignment": {"same_fields": "as alignment, recorded during original run"}
  },
  "outputs": {
    "baseline_cns_sha256": "original-output-digest",
    "source_seed_bed_sha256": "original-output-digest"
  }
}
```

未来来源运行可在开始时调用 `ascompare.concordance_cnv._sampled_identity(path)` 保存这种有限输入指纹，并在实际完成后记录两个小型结果的 SHA-256。历史记录采用其他身份格式时，需要根据**已留存的原记录**转换或单独核验，当前程序不会在缺少证据时自动认定其相同。抽样指纹本身也不是完整 BAM 校验值。

校准执行的检查为：来源 CNVkit 版本与本次实际版本相同；build 和参考策略相同；原始肿瘤和正常 BAM 抽样身份与当前输入一致（`flat` 只检查肿瘤）；CNS 和候选 BED 的完整 SHA-256 与原运行记录相同；原候选及 CNS 对请求种子的覆盖均达到 `min_seed_overlap`。满足来源条件后，本次局部流程仍保留候选记为 `calibration_status=passed`；不能再保留记为 `failed`；来源缺失或不匹配记为 `unverified`。

独立校准只发生在运行版本与样本来源 Suite 相同的对照重跑，标为 `calibration_role=source_control`。`passed` 的范围是 `selected_source_candidate_retention_only`，表示原来源版本对所选候选的再保留。目标版本行标记为 `candidate_observation`，其保留或漏种均作为结果记录。

条件性联合判定使用独立来源校准，要求双方有完整可评价候选、共同的匹配设置和参考策略，且使用各自锁定的 CNVkit 配方。跨版本的结果仍适用于本次已知候选和局部设置。原数据与同版本 AA＋AC 重跑另做来源复现检查。

来源对照保存 `calibration_evidence.json`，包括判据、来源文件哈希、逐种子原始候选覆盖和本次保留率。缺少可追溯历史证据时仍生成描述性差异和图表。

## 合成集成测试

`tests/test_concordance_cnv.py` 包含默认运行的独立测试，以及可选真实 CNVkit 测试。后者自行生成微型参考、两个微型 BAM 与索引，只读取这些人工数据：

```sh
CNVKIT_SMOKE_PYTHON=/absolute/cnvkit-env/bin/python \
CNVKIT_SMOKE_VERSION=0.9.13 \
CNVKIT_SMOKE_REFERENCE_MODE=matched_normal \
python -m unittest discover -s tests -p test_concordance_cnv.py -k test_actual_cnvkit -v
```

从仓库根目录运行，需先安装比较程序及其测试依赖，并准备已安装 CNVkit、R 和 DNAcopy 的组件环境。测试验证有限 BED 的肿瘤/正常 coverage、reference、fix、CBS 和 I/O 预算审计；它不等于真实样本性能或完整 AA 种子筛选已验证。

若解释器来自借用基础科学依赖的新 venv，可用 `CNVKIT_SMOKE_EXTRA_PATH=/absolute/base-env/bin` 指定 Rscript 等外部工具路径。`CNVKIT_SMOKE_REFERENCE_MODE=flat` 测试只有肿瘤数据的分支。`CNVKIT_SMOKE_ARTIFACTS=/absolute/new/directory` 保留全部人工输入、输出和审计，路径必须不存在。

CNVkit 源码身份核验包含真实 `cnvlib/skgenome` 导入位置、安装发行版 RECORD 对所有带哈希 Python 文件的核对、源码清单指纹、完整 `pip freeze`、R/DNAcopy 版本。修改过源码、错误导入路径或无可核验的 RECORD 会明确拒绝。缓存同时包含 I/O 策略、guard 源码、实际依赖、局部参数、参考策略和输入抽样身份；收紧预算不会悄悄复用旧缓存。

发布包提供便携安装/测试脚本。基础 Python 应已有 CNVkit0.9.14 所需的科学依赖、R/DNAcopy 和本程序测试依赖；安装脚本只在指定新目录建立隔离 venv，不改基础环境，已有环境或安装证据会拒绝覆盖。`--no-deps` 是明确复用基础依赖，不能据此认为所有历史分段算法的依赖已安装。

```sh
python scripts/prepare_cnvkit_smoke_envs.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.10 --legacy-cbs-dependencies
python scripts/prepare_cnvkit_smoke_envs.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.13
python scripts/run_cnvkit_smoke.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.10
python scripts/run_cnvkit_smoke.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.13 --reference-mode flat
```

安装报告保存公共包来源和SHA-256；测试脚本通过 `unittest discover` 只选择人工集成测试，并将人工输入、所有结果、guard审计和日志复制到新的 attempt 目录。失败尝试也保留。脚本不启动真实样本任务。

原生 `bedcov` 的 BED 查询长度计入预算，但其 C 层实际 reads 数不能从返回值直接观测。审计中的 Python fetch reads 计数不能冒充对原生 `bedcov` 全部 reads 的计数；有限 BED 总长度和查询预算仍严格执行。

## 已完成的人工组件验证

2026-11-17 的实际执行结果如下；全部是测试脚本自行创建的人工输入，不是患者结论。

| CNVkit | 参考策略 | 实际覆盖/参考/fix/CBS | guard审计 |
|---|---|---|---|
| 0.9.10 | matched_normal | 通过 | 7个受限子进程（含2个身份探针），blocked=0；肿瘤/正常查询合计分别41/59 kb |
| 0.9.13 | matched_normal | 通过 | 7个受限子进程，blocked=0；41/59 kb |
| 0.9.13 | flat | 通过 | 6个受限子进程，blocked=0；肿瘤41 kb，没有正常样本coverage |
| 0.9.14 | matched_normal | 通过 | 7个受限子进程，blocked=0；41/59 kb |

0.9.10 和0.9.13已通过新增的发行版RECORD和真实导入位置核验。0.9.10保留原源码，在新建隔离环境固定 NumPy1.26.4、pandas2.2.3、SciPy1.14.1、pytz2025.2；pomegranate1.1.2仅满足原模块的无条件导入。本基准固定使用CBS，**此环境未验证也不支持旧pomegranate HMM接口**。其他所需依赖来自明确记录的基础环境，R4.5.3、DNAcopy1.86.0；pip来源、包哈希、完整运行依赖和源码指纹均留档。该依赖组合解决了0.9.10在缺少pomegranate或现代NumPy/pandas/SciPy组合下的已观察兼容问题。

这些环境验证了指定历史CNVkit代码的CBS局部流程，**不声称重建了历史容器的全部依赖发行版，也不能把依赖同时不同的差异全归因于CNVkit源码版本**。AA原生种子筛选目前在本模块经过命令/状态/解析测试，尚未用真实患者完成这条局部联合流程；实际候选保留结论仍需真实条件验证。
