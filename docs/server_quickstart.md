# 服务器快速开始：AmpliconSuite Cross-Version Concordance

1. 在仓库根目录按README创建并激活Python环境、安装程序，再执行asconcordance demo检查导入和报告。
2. 准备已有结果、来源Suite版本和参考版本声明；不要上传患者文件。
3. 选择results_only或with_bam示例，填写实际本地路径。
4. 需要调用算法时准备锁定源码和对应运行环境，先check再run。

## 无BAM，比较已有输出

```sh
asconcordance check --manifest concordance.results_only.json --out /output/existing_results
asconcordance run --manifest concordance.results_only.json --out /output/existing_results
```

同一样本的其他版本输出放在comparison_results列表中。若只有一个版本输出，其他AA版本明确未验证；要比较AC，使用concordance.ac_only.json配置参考资源和AC解释器/源码。

## 有BAM，局部上下游

使用concordance.with_bam.json。每例提供已有种子BED、原AA/AC目录、来源Suite版本、BAM和现有索引。上游需明确背景BED。默认matched_normal参考模式还需配对正常BAM；只有肿瘤BAM时可明确设 `local_cnv.reference_mode: "flat"` 并去掉normal字段。不会自动由matched_normal切换flat。不得用种子BED冒充正常背景，flat结果仍须校准，不能宣称复现了原正常参考。

内置版本不包含二进制环境。可用锁定源准备结果与运行模板生成多版本配置：

```sh
asconcordance prepare-panel --out /work/locked_sources
python scripts/build_concordance_manifest.py \
  --input /config/base_manifest.json \
  --prepared /work/locked_sources/prepared_panel.json \
  --runtime-templates /work/runtime_templates.json \
  --cnvkit-runtimes /config/cnvkit_runtimes.json \
  --out /config/concordance.json
asconcordance check --manifest /config/concordance.json --out /output/local_concordance
asconcordance run --manifest /config/concordance.json --out /output/local_concordance
```

生成器仅填来源明确的路径，不能自动证明环境已安装正确。runtime_templates使用 `scripts/create_runtime_envs.py` 返回的AA/AC模板，精确CNVkit路径由独立清单提供，禁止0.9.14静默替代0.9.10/0.9.13。输入和输出manifest需同目录，以保留相对路径含义。

## 参数和结论

`sampling.max_seeds`控制每例下游种子数，默认10；`random_seed`固定抽样。全部输入候选仍用于上游。`local_cnv`控制两侧扩展、bin、总窗口预算和种子覆盖匹配；这些是可配置试验条件，不是官方兼容标准。

`statistics.margins`四项：detection_discordance、label_discordance、interval_distance、class_footprint_distance；`upstream_loss_margin`是允许的样本平均候选丢失比例。它们默认未设。只有在分析前依据用途确定，才填入0到1之间的允许差异；不提供任意数字来保证通过。

默认min_samples=20是保守实施门槛，不是“20例必然有足够检验效能”。程序仍要求保守置信界位于允许范围内。小样本可用于描述初步差异和检查工程正确性，通常不足以证明总体等效。详细解释见statistical_methods.md。

先看report.md中的执行覆盖率和限制，再看天梯图。compatibility_combined.tsv包含“证据不足”，不得强制改成二分类。来源身份、失败、类别不支持、局部校准和来源复现检查都有独立字段。

## 续跑和迁移

重复相同run命令可复用验证过的CNV、AA和AC重分类结果。缺件或哈希变化会重试；输入或锁变更需新输出目录。修改允许差异属于新的分析计划，应保留旧报告。

```sh
asconcordance report --out /output/local_concordance
```

此命令仅根据已有记录重画图，不读取BAM、不重新运行算法。真实样本应使用 `asconcordance` 入口；旧 `ascompare run/check` 已禁用。原始数据、真实结果、环境配置和许可证应存放在仓库外，不能随Git提交上传。
