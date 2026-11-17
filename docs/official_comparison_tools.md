# 官方比较工具复用评估

**当前范围：**本仓库保留以下独立适配器。`asbench` 不自动执行 `compare_completed_records` 或把官方距离合入兼容判定；主报告使用可审计的检出、多标签及区段指标。需要路径/断点诊断时，按本页独立命令对已完成结果运行；不读取BAM或启动全基因组操作。

核查日期：2026-11-17。评估的官方工具为 **EC-comparator** 与 AmpliconClassifier 的 **`results_comparator.py`**。本项目为独立工具，不代表官方对本基准或统计结论的认可。

| 工具 | 固定源码 | 适合本项目的用途 | 不能替代的部分 |
|---|---|---|---|
| [EC-comparator](https://github.com/AmpliconSuite/EC-comparator) | v0.1.3；`13cb961dc10aca54e44b71da38d64d817f91c5ac` | AA 路径/环的覆盖、拷贝数、片段、环数量和断点匹配诊断 | AA/AC 多版本启动、运行失败状态、样本配对统计及兼容判断 |
| [AC results_comparator.py](https://github.com/AmpliconSuite/AmpliconClassifier/blob/2ec0c146f78e74b4857860fc30b5e8a1de7920b7/scripts/results_comparator.py) | 当前项目固定参考的 AC2.0.1；`2ec0c146f78e74b4857860fc30b5e8a1de7920b7` | 官方结果表中 feature 的增检、丢失、分类变化、拆分、合并及 Jaccard | 原始输入检查、跨版本缺失类别、多种状态区分、置信区间和等效/一致性检验 |

实现：提供可选官方诊断适配器 `ascompare.official_compare`，调用原始官方程序，保留全部原始输出；主程序继续负责版本固定、配对状态、样本聚类统计与允许差异范围。没有重新实现官方断点匹配、CN 距离、feature 拆分/合并算法。

## 官方距离范围与兼容阈值核查

固定 v0.1.3 的 [README “Distances explained”](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/README.md#distances-explained) 将六种归一化距离描述为0到1，越小越相似，并将默认总成本描述为0到5。源码 [distances.py](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/eccomparator/metrics/distances.py) 另输出未经归一化的 `cn_hamming_dist`，其单位是差异覆盖碱基数，不能当作0到1距离。`total_cost` 是 [compare.py](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/eccomparator/metrics/compare.py) 按 [configs.py](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/eccomparator/configs.py) 权重求和，本适配器保持六项开启和权重总和5；禁用指标后总成本范围会改变。

**在本次核查的固定版本 README、CLI、configs、距离实现以及 AC results_comparator.py 中，没有找到可用于“跨软件版本兼容通过”的官方分数界限。** 这是对上述已检查来源的结论，不能排除其他或未来文档提出标准。断点1000/3000距离、1000bp feature匹配以及 AC 自身分类阈值都是算法参数，不是兼容性容许差异，不能挪用。

因此逐样本版本图可直接展示官方距离（越低越近）及主程序区段Jaccard（越高越近），不任意换算“真实性概率”。兼容表仍使用预先配置的指标容许差异与样本聚类统计；未预设容许差异或独立样本太少时保持证据不足。`official_metrics.tsv` 导出每个 sample/seed/version 对的原生标量、单位、方向、文档范围和状态；缺失/不支持保留空值，越界值单列标记，不强行裁剪或填0。

## results_comparator.py 的实际语义

输入是 `*_result_table.tsv`，不是 `*_amplicon_classification_profiles.tsv`。其必要列为 `Sample name`、`Feature ID`、`Classification`、`Location`。每次必须比较相同样本和种子；两个结果表中的样本身份需一致，不能把不同患者的文件直接拼成同一样本。

官方程序以至少 **1000 bp** 区段交集连接两个版本的 features，按二部图连通分量报告 `GAIN`、`LOSS`、`CLASS_CHANGE`、`SPLIT`、`MERGE`、`COMPLEX_MERGE_SPLIT` 或 `STABLE`。1000 bp 在该版本脚本内固定，没有对应命令行参数。适配器固定 run1 为参考、run2 为候选，因此 GAIN 指候选额外 feature。

`STABLE` 只表示一对一且分类相同，不能解释为结构完全相同或统计兼容。无坐标、畸形坐标或部分坐标错误的 feature 在原工具中会被全部或部分忽略，适配器提前拒绝这类错误；原工具本身无法区分没有结果和真正未检出。因此只向它传入已确认成功且可评价的任务。原工具的 per-feature Jaccard 与主程序的每个种子全部分类区段并集 Jaccard 是不同层级，不能直接混成一列。

还有多标签注意事项：官方匹配只看区段重叠，不要求机制类别相同。如果两个版本都有彼此重叠的 ecDNA 与 BFB/FAN features，即使完全不变，也可能形成多对多连通分量并被命名 `COMPLEX_MERGE_SPLIT`。这是官方的重叠关系标签，不能据此认定发生了真实结构拆分/合并；应同时查看主程序多标签差异、区段变化和原始 AA 结构。

合成 native-format 测试已执行原脚本：一个 ecDNA feature 拆为 ecDNA 与 BFB 两个 feature，正确输出 `SPLIT`；另一个不变线性 feature 输出 `STABLE`，两者 Jaccard 均为1。对应的可复用人工输入保留在 `examples/official_synthetic/`；不能把这些结果当成患者验证。

## EC-comparator 的输入转换与局限

官方输入是带表头的 BED-like 路径表，至少含 `#chr,start,end,strand,circ_id,estimated_cn`，v0.1.3 还识别 `iscyclic`。顺序与方向携带结构信息；不能使用合并后的 BED 区段替代路径。

源码检查发现官方辅助脚本 `cycles_to_bed_modified.py` 会丢弃 AA 的0号源/汇节点，未输出 `iscyclic`；EC-comparator 缺少这一列时默认所有路径为环。辅助脚本还直接复制 AA 的闭区间末端。因此本项目仅编写必要的**格式转换**：保留每次片段出现、方向、顺序和原始 Copy_count；AA末端+1转为半开区间；两端0号节点的路径明确写 `iscyclic=False`；不认识的内部0号节点结构拒绝转换。多个 amplicon 的 cycle ID 加文件序号，防止同号混合。

转换后的整套 AA 路径包含未被 AC 判为 ecDNA 的结构，报告只能称“重建路径/环的结构诊断”。相同 seed 的所有 amplicon 一起比较，保留相连区域，不按 amplicon 编号强行一一匹配。不能把该工具内部 `true/reconstructed` 变量名当作已知生物学真值；本项目只使用“固定参考/候选”标识。

多 amplicon 可以安全表示为**独立路径/环的集合**：官方 [breakpoints.py](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/eccomparator/metrics/breakpoints.py) 按 `circ_id` 分组产生断点，适配器的文件序号前缀保证不同 amplicon 之间不会新建连边。官方 [features.py](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/eccomparator/metrics/features.py) 聚合覆盖与CN、[distances.py](https://github.com/AmpliconSuite/EC-comparator/blob/13cb961dc10aca54e44b71da38d64d817f91c5ac/eccomparator/metrics/distances.py) 按基因组bin计算片段/路径数量差异。因此部分重叠或重复的 amplicon 分解会累加贡献；它比较的是软件重建集合，不能声称拼成一个实体完整分子，也不能把累计路径CN当作独立测得的基因组真值。

该工具的 `main()` 捕获异常后可能仍退出0。适配器同时要求产生有效且数值有限的 `metrics.json`，不能仅看退出码。固定默认断点候选距离1000、代价距离3000、min-cn0，原参数与 metrics 原样保存。这些诊断阈值并非本项目的统计兼容容许误差。

## 运行

源码和环境均在本地，所有患者文件留在本地。两个诊断都要求新的空输出目录；保存源码提交、文件哈希、原始输入副本或路径转换、实际 Python/依赖版本、参数、日志和官方原始输出。此步骤不重跑 AA，不做全基因组注释。

先按主 README 进入本程序目录，替换下列本地路径：

```sh
python -m ascompare.official_compare results \
  --first /path/reference_result_table.tsv \
  --second /path/candidate_result_table.tsv \
  --source /path/locked/AmpliconClassifier-2ec0c146f78e \
  --python /path/python-with-pandas \
  --out /path/new-official-feature-comparison \
  --reference-version suite-1.6.2 --candidate-version suite-1.6.1

python -m ascompare.official_compare ec \
  --first /path/reference_amplicon1_cycles.txt /path/reference_amplicon2_cycles.txt \
  --second /path/candidate_amplicon1_cycles.txt \
  --source /path/EC-comparator-13cb961dc10a \
  --python /path/ec-env/bin/python --extra-path /path/bedtools/bin \
  --out /path/new-official-structure-comparison \
  --reference-version suite-1.6.2 --candidate-version suite-1.6.1
```

EC-comparator 建议使用独立 Python3.10 环境，安装固定源的 `requirements.txt`（numpy1.26.4、pandas2.2.2、scipy1.13.1、pybedtools0.12.0、pyranges0.0.129 等）并提供 bedtools。即使关闭绘图/报告，官方代码启动时仍导入相关库，因此需要完整依赖。适配器检查官方精确依赖要求并记录全部实际依赖；缺失或不符输出 `unsupported`，不会替换算法或假装完成。

```sh
git clone --depth 1 --branch v0.1.3 https://github.com/AmpliconSuite/EC-comparator.git /path/EC-comparator-13cb961dc10a
python3.10 -m venv /path/ec-env
/path/ec-env/bin/python -m pip install -r /path/EC-comparator-13cb961dc10a/requirements.txt
```

开发验证曾在macOS ARM64上安装并运行固定官方 EC 代码的合成测试。原生扩展的编译器、SDK和bedtools需与部署平台匹配；不要复制其他机器的二进制环境。

实际官方 EC 合成测试中，两个相同的双片段环得到全部距离和 total_cost 为0，两个断点匹配，退出码0并通过 `metrics.json` 校验。读取器的 `total_cost_description` 权重字典处理纳入回归。**这只验证官方程序与适配器在该人工条件下能运行，不验证真实样本准确性或版本兼容性。** 运行本仓库常规测试不会自动下载或启动上述外部软件。

官方源码清单锁在 `ascompare/official_tools.lock.json`；AC 工具使用项目原 `versions.lock.json` 的参考版本源码验证。版本身份来自固定提交与源码指纹核验，不根据镜像名推定。

格式/保护测试覆盖末端坐标、重复遍历、方向、开路径、跨amplicon同号、缺失片段、源/汇异常与缺坐标拒绝。真实病例的兼容结论始终来自主程序的有信号区域、整样本统计和预先设定误差界限；官方诊断分值不能转译为真实性概率。

## 与自动运行结果衔接

程序接口 `compare_completed_records(records, out, toolconfig, reference_version)` 只比较同一 `sample_id + seed_id` 的已成功任务，生成 `out/official_comparisons/index.tsv`、`index.json` 和每对版本的独立原始诊断目录。`toolconfig` 格式如下：

```json
{
  "results": {"source": "/path/locked/AmpliconClassifier-2ec0c146f78e", "python": "/path/python-with-pandas"},
  "ec": {
    "source": "/path/EC-comparator-13cb961dc10a",
    "python": "/path/ec-env/bin/python",
    "extra_path": ["/path/bedtools/bin"],
    "allow_multiple_amplicons": true
  }
}
```

未配置的工具不运行。历史版本没有 native result table 时明确记 `missing`；主流程失败记 `not_evaluable`，不生成假阴性表。对于使用 `__whole_sample__` 标识的整样本已完成结果（诊断目录映射为 `whole_sample`），EC 默认允许比较多个 amplicon 的独立结构集合；对于固定种子结果，需设置 `allow_multiple_amplicons=true` 才进行多 amplicon 集合比较，否则记 `unsupported`。也可以显式设为false关闭集合比较。不会按输出文件顺序连接两个 amplicon，也不会把集合分值说成一一对应结构相同。

适配器回归覆盖上述保护、官方 metrics schema、标量指标单位/方向/状态及绘图失败隔离，见 `tests/test_official_compare.py`。显式调用自动配对 helper 时，会在其输出的 `official_comparisons/` 下生成 `official_metrics.tsv`；这不是 `asbench` 默认输出的一部分。

## 逐样本官方距离天梯图

`ascompare.official_plots.plot_official_metrics(out)` 读取上述标量表，按 sample和seed分别绘图：六项原生归一化距离置于六个独立坐标轴（0–1），官方总成本另绘0–5图。均为越低越接近固定参考，保留实际参考版本名；单样本只画观测值，不加样本总体置信区间。不可评价状态用灰带和文字表示，不填0；不同种子不合并成一个样本结果。

自动配对 helper 在保存标量表后自动生成 PNG、PDF、SVG 及 `official_plots.tsv`，无需单独启动绘图。`plots_status.json` 记录绘图是否成功；绘图异常不会改写官方原始比较状态或统计兼容结论，刷新前清空图表索引以防误用旧图。参考版本无重建结构时使用灰色状态，不画0距离锚点。在源记录声明不完整人类参考时图中保留范围提示；人工测试显示明确测试标识。绘图回归见 `tests/test_official_plots.py`，包含整样本标识匹配、缺失结构和人工标识保护。
