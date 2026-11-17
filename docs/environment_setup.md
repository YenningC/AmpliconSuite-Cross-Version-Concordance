# 0.2 环境与依赖

比较/报告环境需Python≥3.10及requirements.txt中的NumPy、SciPy、Matplotlib、pysam。AA/AC/CNVkit各自使用显式解释器，不能由PATH猜版本。安装源码包或wheel后使用asbench；ascompare仅用于公共源码准备等辅助功能，不用于真实样本启动。

```sh
python3 -m venv .venv
.venv/bin/python -m pip install .
.venv/bin/asbench demo --out /absolute/output/SYNTHETIC_demo
. .venv/bin/activate
```

后续命令在已激活的比较程序环境、仓库根目录运行。示例中的 `/absolute`、`/srv` 路径需要替换为自己的目录。

无BAM、已有多套AA/AC输出时不需要安装算法运行环境。重跑AC需对应版本源码、AA参考资源及该AC依赖；不需要AA解释器和BAM。

## 获取确切官方源码

```sh
ascompare prepare --version suite-1.6.2 --out /absolute/work/locked_sources
python scripts/prepare_panel.py --out /absolute/work/locked_sources
```

第一条只获取一个版本；第二条获取锁内面板，生成prepared_panel.json。只下载公共软件，不读取患者数据。锁内源码归档、commit及逐文件清单均核对，不用镜像标签推测版本。算法不在安装包内，以免未经核验的软件替代原版。

## AA与AC环境

可复用既有已核实的环境，或根据官方对应版本依赖建立隔离环境。推荐服务器使用Conda/Mamba提供科学库、pysam、CBC、R和DNAcopy；历史AA还需Mosek及有效许可，新AA支持Clarabel。不得为绕过历史许可而改成新AA。

已有 `scripts/create_runtime_envs.py` 可从包含科学依赖的Conda基础Python创建reference-env和ac200-env。该脚本返回的模板含CNVkit0.9.14，仅将其AA/AC环境部分用于本面板；局部CNV另设0.9.10与0.9.13解释器。它不是复刻所有历史依赖的保证。

参考Suite1.6.2要求AC2.0.1/BFBArchitect1.0.2；Suite1.6.0/1.6.1要求AC2.0.0/BFBArchitect1.0.1。现代与旧版BFB环境不可混为相同依赖。PuLP/CBC的处理器架构、Mosek许可、BFB实际回退后端均需检查。程序保存依赖清单、源码核验和探针；启动成功不等于真实样本验证通过。

本程序直接编写AC输入表，不调用需要Bash关联数组的make_input.sh。实际AA/AC文件数量、完成日志和解析结果均会检查，不仅依赖退出码。

## 两个原生配方CNVkit环境

前七个内置Suite的容器配方使用CNVkit0.9.10，后三个用0.9.13。Conda/独立安装最低约束不足以证明实际版本，必须核验。

建立两个隔离环境，分别安装确切的 `CNVkit==0.9.10` 和 `CNVkit==0.9.13`，配齐各版本依赖与R/DNAcopy。不能删除导入、伪造模块或静默改成最新版。发布脚本支持复用明确指定的基础环境，且拒绝覆盖已有环境或安装证据：

```sh
python scripts/prepare_cnvkit_smoke_envs.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.10 --legacy-cbs-dependencies
python scripts/prepare_cnvkit_smoke_envs.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.13
```

基础Python应已有本程序测试依赖及CNVkit所需的科学库；外部工具R/DNAcopy应可从基础环境的 `bin` 找到。脚本通过 `--system-site-packages` 复用这些依赖，再以 `--no-deps` 安装确切CNVkit。0.9.10的已验证CBS依赖组合为 **NumPy1.26.4、pandas2.2.3、SciPy1.14.1、pytz2025.2、pomegranate1.1.2**。后者只满足旧模块的无条件导入，**此环境不支持或验证旧pomegranate HMM接口**；本程序固定使用CBS。安装报告保存公共包下载来源和SHA-256，不能把 `--no-deps` 理解为所有历史算法依赖均已完整安装。

这些是指定历史CNVkit源码的隔离运行环境，**不是历史整套容器的复刻**。其他科学库和R/DNAcopy来自明确记录的基础环境；不同依赖同时变化时，观察到的差异不能全归因于CNVkit源码版本。实际版本、真实导入位置、发行版RECORD、源码清单和完整pip freeze会保存。

配置示例在examples/cnvkit_runtimes.example.json，明确version/python/script三个字段。程序用该解释器查询版本、计算cnvlib代码清单和pip freeze；不匹配立即失败。若故意用0.9.14作受控诊断，会保留结果并标 `native_recipe_match=false`，不能形成历史原生Suite的联合接受结论。

局部流程需要有限背景BED，默认matched_normal还需配对正常BAM。显式flat模式可仅用肿瘤BAM，由该CNVkit在有限bins上构建平坦参考；不会自动切换，也不等同原始正常参考。仅调用coverage、reference、fix和CBS segment；不调用batch/autobin/access，不重建BAM索引。CNVkit预检和读长日志辅助函数的索引适配会明确记录；不能称完全未适配的容器发行运行。

## 已完成的人工CNV验证与复现

2026-11-17在Apple Silicon、Python3.12开发验证环境下，使用测试脚本自行生成的微型FASTA、BAM及索引，实际完成以下验证：

| CNVkit | 参考策略 | 有限coverage→reference→fix→CBS | guard审计 |
|---|---|---|---|
| 0.9.10 | matched_normal | 通过 | 7个受限子进程（含2个身份探针），blocked=0；肿瘤/正常累计请求41/59 kb |
| 0.9.13 | matched_normal | 通过 | 7个受限子进程，blocked=0；41/59 kb |
| 0.9.13 | flat | 通过 | 6个受限子进程，blocked=0；肿瘤41 kb，无正常样本coverage |
| 0.9.14 | matched_normal | 通过 | 7个受限子进程，blocked=0；41/59 kb |

0.9.10使用上述兼容依赖组合；0.9.13/0.9.14使用记录的现代科学依赖，外部R4.5.3、DNAcopy1.86.0。验证期间发现的pomegranate缺失及NumPy/pandas/SciPy兼容问题通过显式依赖配置解决，没有修改CNVkit算法源码。详细校准判据和验证边界见 [局部CNV校准](local_cnv_calibration.md)。

```sh
python scripts/run_cnvkit_smoke.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.10
python scripts/run_cnvkit_smoke.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.13
python scripts/run_cnvkit_smoke.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.13 --reference-mode flat
# 基础环境实际安装CNVkit0.9.14时，验证受控诊断版本：
python scripts/run_cnvkit_smoke.py \
  --base-python /srv/base-env/bin/python --out /srv/work/cnvkit-smoke \
  --version 0.9.14
```

脚本使用 `unittest discover` 只选择人工集成测试，保留每次attempt的命令、日志、人工输入/输出和读取审计；已有结果不覆盖。若CNVkit解释器位于其他位置，可显式指定 `--cnvkit-python`。它不启动患者任务。

以上不是Linux服务器全流程验证，也不构成真实样本保留率或兼容结论；局部CNV与真实AA原生种子筛选、AA/AC重建的联合患者流程尚未由这些人工测试验证。部署Linux后应先运行相同人工测试，再根据真实样本条件验证。C层bedcov的实际reads数不可直接观测；审计中的Python fetch计数不能冒充全部reads计数，有限BED和请求长度预算仍执行。

## 资源与路径

提供已下载、固定且与比对参考一致的AA_DATA_REPO；新程序不替用户在线更新参考，不上传输出。FASTA和BAM的索引必须已经存在。BAM头中的参考长度只提供有限一致性证据；不扫描BAM或整份FASTA来声称核实每一条序列。CRAM转换不在本版范围。

环境可在新目录安装，原数据和原研究环境不可覆盖。运行日志只记录许可证路径或服务器配置，不收集许可证内容。本仓库只提供当前有限区域流程的运行说明。
