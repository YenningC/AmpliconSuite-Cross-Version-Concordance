# 局部 CNV 候选保留的来源校准

局部基准只回答：已有候选区域，在有限窗口内重新估算 CNV 并经过对应 AA 种子筛选后，仍有多少区段被保留。默认窗口为候选两侧各 100 kb，另加用户提供的可信背景 BED；所有窗口合并后的长度默认最多 25 Mb。所有 BAM 均要求已有索引；不会扩大为全基因组扫描。

`local_cnv.reference_mode` 默认 `matched_normal`，要求配对正常 BAM。只有肿瘤 BAM 时，可以**明确**配置 `flat`：程序使用原生 `cnvkit reference -t <有限 targets> -f <FASTA>` 构建平坦参考，不读取正常样本，也不静默切换参考策略。两种策略仍都要求可靠背景 BED。`flat` 不是正常样本，不能声称完全复现原来的正常参考/归一化；未经同策略校准时保持 `unverified`。策略写入每条结果、provenance 和缓存身份，切换策略会重新计算。

CNVkit 原生 coverage、reference、fix、CBS segment 计算保留；WGS 风格禁用 edge correction。coverage 中原先从 BAM 开头取样的排序/读长辅助步骤，改为排序头检查和有限 BED 内最多 1,000 条 reads 的读长取样，明确记录 `bounded_indexed_auxiliary_probes`。不修改安装的上游源码。CBS 跨未测窗口间隙的输出仅裁回实际测量窗口，避免把未测区段当作 CNV 证据。所有局部 CN 区段，包括中性区段，都传入该版 AA 的 `amplified_intervals.py`；该步骤保留其必要的有限、索引式背景估计，并受读取预算限制。

局部中位数归一化、背景和分段上下文与全基因组不同。**候选保留率可以描述，但未经校准不能当作整条 Suite 流程兼容证据。** `min_seed_overlap` 默认 1.0（完全覆盖），是可配置的基准匹配阈值，不是 CNVkit/AA 的官方阈值。它计算合并后的保留区段与原种子的交集长度/原种子长度，重叠不会重复计数。

`builtin_cnvkit_map.json` 来自锁定 Suite 源码中的容器构建配方。1.5.2 及以前的七个版本对应 CNVkit 0.9.10，1.5.3 和三个 1.6.x 版本对应 0.9.13。程序核查实际 CNVkit 版本和导入源码。0.9.14 等受控版本记录为 `controlled_cnvkit`，报告保留其观察结果和实验范围。

## 何时能校准

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

未来来源运行可在开始时调用 `ascompare.benchmark_cnv._sampled_identity(path)` 保存这种有限输入指纹，并在实际完成后记录两个小型结果的 SHA-256。历史记录采用其他身份格式时，需要根据**已留存的原记录**转换或单独核验，当前程序不会在缺少证据时自动认定其相同。抽样指纹本身也不是完整 BAM 校验值。

校准执行的检查为：来源 CNVkit 版本与本次实际版本相同；build 和参考策略相同；原始肿瘤和正常 BAM 抽样身份与当前输入一致（`flat` 只检查肿瘤）；CNS 和候选 BED 的完整 SHA-256 与原运行记录相同；原候选及 CNS 对请求种子的覆盖均达到 `min_seed_overlap`。满足来源条件后，本次局部流程仍保留候选记为 `calibration_status=passed`；不能再保留记为 `failed`；来源缺失或不匹配记为 `unverified`。

独立校准只发生在运行版本与样本来源 Suite 相同的对照重跑，标为 `calibration_role=source_control`。`passed` 的范围是 `selected_source_candidate_retention_only`，表示原来源版本对所选候选的再保留。目标版本行标记为 `candidate_observation`，其保留或漏种均作为结果记录。

条件性联合判定使用独立来源校准，要求双方有完整可评价候选、共同的匹配设置和参考策略，且使用各自锁定的 CNVkit 配方。跨版本的结果仍适用于本次已知候选和局部设置。原数据与同版本 AA＋AC 重跑另做来源复现检查。

来源对照保存 `calibration_evidence.json`，包括判据、来源文件哈希、逐种子原始候选覆盖和本次保留率。缺少可追溯历史证据时仍生成描述性差异和图表。

## 合成集成测试

`tests/test_benchmark_cnv.py` 包含默认运行的独立测试，以及可选真实 CNVkit 测试。后者自行生成微型参考、两个微型 BAM 与索引，只读取这些人工数据：

```sh
CNVKIT_SMOKE_PYTHON=/absolute/cnvkit-env/bin/python \
CNVKIT_SMOKE_VERSION=0.9.13 \
CNVKIT_SMOKE_REFERENCE_MODE=matched_normal \
python -m unittest discover -s tests -p test_benchmark_cnv.py -k test_actual_cnvkit -v
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
