# 来源复现与局部 CNV 校准

源结果由历史运行产生。“来源版本重跑”指 Concordance 在本次冻结输入上用同一 Suite 版本再次运行。两者需要独立核对。下面的对照可以利用旧结果中的小文件，前两步不需要 BAM。

## 1．定位来源 AA＋AC 的差异

1. 在 `source_bridge/<sample>.json` 找到未复现的种子、增减的标签和区段。查看原样本的 `*_run_metadata.json`、原生日志以及当前 `downstream_jobs/.../stages.json`。固定 build、种子、AA 的 EXPLORE/FULL、下采样、solver，以及 AC 最小环长度、最小流量和参考数据版本。AC 的 `min_size` 和上游种子长度 `cnsize_min` 是不同的参数。
2. 按种子和重建区段对应原始与重跑的 amplicon。分别核对 AA graph、cycles、AC profile、feature BED。graph/cycles 文本字节有差异时检查具体字段；文本相同、区段相同、标签相同是不同的核对结果。保留整套相连区段。
3. 用已归档的两套 AA 输出各做同版本 AC 重分类。在每套图内一次只改变一个有记录差异的参数，对照原 feature BED 和标签。参考资源、输入图和其他参数固定。可以沿用 results_only 的 AC-only 模板，每组条件使用独立目录。原生默认值要从对应版本源码或原运行配置核实。若来源确实用了 5,000 bp，下次统一比较时可明确填写 `"ac_parameters": {"min_size": 5000, "min_flow": 1.0}`，并对每版使用同一预设。历史 50,000 bp 面板的结论继续对应它自己的设置。
4. 将选定设置写入新清单，再通过来源桥接。报告中的 AC 固定图对照只能解释分类步骤。它不补齐一个失败的 AA 重建，也不能替代跨患者验证。

## 2．把 CNV 差异定位到具体步骤

1. 保留**运行当时**的原始输入身份、CNVkit/AA 版本、build、参考策略、targets、分箱、校正、分段和筛种参数，以及原 CNS 和最终 BED 校验值。配置字段与既有格式见 [局部校准输入规范](local_cnv_calibration.md)。仅找到旧命令可以还原参数；缺少的原运行身份仍标为 `unverified`。
2. 优先读取已有的原始 CNR（分箱）和 CNS（分段），分别和局部流程的 CNR、CNS 在候选区域的共同有效碱基上配对。记录每个种子的可比长度、log2 偏移和跨过同一增益阈值的范围。CNR 已出现整体偏移时先核对正常参考、中心化、GC/重复序列校正和分箱。CNR 相近而 CNS 不同，再查 CBS 窗口、邻近边界与未测间隙。CNS 相近而最终 BED 不同，再查原生 AA 筛种条件、匹配阈值和覆盖比例。此顺序用于定位疑点。
3. 在新目录运行 0.4.2 的 `with_bam`。程序按每个样本的实际来源 CNVkit 建立局部标准，对每个目标 CNVkit 固定同一来源 AA 筛种器。若原用 0.9.14，打开 `local_cnv_comparator_samples.tsv` 查看 0.9.10、0.9.13 等目标相对局部 0.9.14 的区段交并比、额外保留、遗漏和保留率差。用标准行的保留个数比单独检查原候选的再保留；`baseline_original_calibration_status` 同时考虑历史来源证据。原生 Suite 配方组合继续独立运行。来源标准未执行或局部设计不一致时，保留其缺失原因。
4. 在训练患者上选择局部窗口/背景和候选匹配规则，冻结后用独立患者或预留患者验证。比较所有预定候选并保留失败。只针对原候选的校准支持“候选再保留”这个终点。缺少历史证据时可以建立前瞻性来源基线，保留历史样本的未核实状态。

CNVkit 各阶段的输入输出见[官方流程说明](https://cnvkit.readthedocs.io/en/stable/pipeline.html)。AA 输出与版本说明见[上游项目](https://github.com/AmpliconSuite/AmpliconArchitect)。

### 用保存的结果做分箱与分段对照

在源码目录运行下面的命令，以本地路径替换示例。三个 build 声明须一致，CNR/CNS 可为 gzip 文件。CNS 同时受双方实测 CNR 范围约束。每例输出到一个新文件。

```sh
python3 scripts/compare_cnv_profiles.py --seeds-bed input/source_seeds.bed \
  --source-cnr input/source.cnr.gz --source-cns input/source.cns \
  --local-cnr local/tumor.cnr --local-cns local/tumor.cns \
  --source-build hg19 --local-build hg19 --seed-build hg19 \
  --gain-cn 4.5 --out control/profile.tsv
```

### 历史 AA 触及唯一访问区段上限

符合条件的旧版种子处理已生成该版原生 `coverage.stats`。程序检查同次运行的样本、BAM、版本、参考与原生测量，复用到该版 AA 的私有参考目录。数值来自种子步骤，重建继续用自身的下采样与向外扩展。检查该次重建目录的 `aa_background.json` 中的模式和 `reuse_confirmed`。缺少有效证据时使用原有受限测量。若仍超出 100 Mb 守护上限，查看 `bam_access.*.json` 的 `first_rejected_unique_request`，记录触发区间与当时的累计唯一碱基数；该次重建仍记为失败。不同种子子集属于新的分析范围，需要对所有版本一致冻结并单独报告。

当历史 AA 以读段名称的 Python 哈希作下采样时，应记录实际哈希种子。0.4.1 对受限 AA 进程使用固定的数字 `PYTHONHASHSEED`，未指定时取 0，并写入该次背景记录与缓存身份。配置为随机值时须先改为可复现的数字。历史运行未记录其有效值时，原始下采样无法仅凭参数声明复现。
