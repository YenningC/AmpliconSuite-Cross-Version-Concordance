# 上传到 GitHub

项目名称为 **AmpliconSuite Cross-Version Concordance**。GitHub 仓库名为 `AmpliconSuite-Cross-Version-Concordance`。简介可填写：

> Bounded, reproducible cross-version comparisons of AmpliconSuite outputs, with per-sample similarity ladders and cautious compatibility assessment.

## 网页上传

1. 解压 GitHub 源码包，进入包含 `README.md`、`pyproject.toml`、`ascompare/` 的那一层目录。
2. 在你的 GitHub 账户创建仓库，自行选择公开或私有。已有README和忽略文件，不需要再让GitHub生成同名文件。
3. 上传这个目录里的文件和子目录，使 `README.md` 位于仓库根目录；不要只把ZIP作为唯一仓库文件上传。
4. 确认隐藏的 `.gitignore` 和 `.github/workflows/tests.yml` 也已上传。图示位于 `docs/assets/`，需要一同上传。
5. 检查首页中英文链接、人工示例图片和安装步骤。包含的Actions只跑人工测试；首次GitHub执行结果需以上传后的实际日志为准。

macOS默认可能隐藏以点开头的文件，可以按 `Command + Shift + .` 显示。若网页不便保留目录和隐藏文件，可用GitHub Desktop提交解压目录的内容。

## 使用 Git 上传

在解压后的仓库根目录执行；先在GitHub创建空仓库，再将最后一条地址替换为你的真实仓库地址。占位地址不是已创建的仓库。

```sh
git init -b main
git add .
git diff --cached --stat
git commit -m "Add AmpliconSuite Cross-Version Concordance"
git remote add origin https://github.com/YOUR_ACCOUNT/AmpliconSuite-Cross-Version-Concordance.git
git push -u origin main
```

源码包只包含程序、说明和人工示例。不要把原研究目录、真实配置、患者结果、BAM、参考库、环境或Mosek许可证一起上传。`.gitignore`是辅助措施，不能自动识别所有敏感文件。

## 许可证与发布附件

项目原创代码和文档采用 [MIT 许可证](LICENSE)。重新分发时保留版权和许可证声明。

可以在仓库 Release 中另附 wheel、源码 ZIP/TAR.GZ 及 SHA-256 清单。0.3.0 新增报告和解析能力。请等候 Actions 测试完成，并保留 0.2.0 历史发布包以便复现旧分析。
