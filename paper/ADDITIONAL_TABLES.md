# 附加模拟：从合并结果重建论文表格

本包包含已生成的两张 LaTeX 表、两页 PDF 预览、完整精度 CSV，以及生成它们的 base R 脚本。表格使用 N=1,000、B=300、锁定的 delta=0.1 配置。制表不运行估计器或 bootstrap，不需要 cluster、CVXR 或 MrDAG。

## 1. 放置文件

将 `additional_tables_delta01_b300.zip` 解压到本机 `MRrr-rebuild` 仓库根目录，使脚本位于 `paper/scripts/30_make_additional_simulation_tables.R`。将从 cluster 下载的 `additional_results_delta01_b300.tar.gz` 也放在仓库根目录。

在 Windows Git Bash 执行：

```bash
cd "/d/Users/YuexiangPeng/Documents/UW/Research/Ye Ting/MRrr-rebuild"
sha256sum additional_results_delta01_b300.tar.gz
```

此结果包的 SHA-256 应为：

```text
681840cee148b2fdb5044d6157c0771e3816cb466a2c4350fd5ed34f0e15b004
```

下面提取合并结果和恢复记录。原结果包还包含代码快照；这条命令只提取所列结果目录，当前仓库代码由 Git 管理。

```bash
tar -xzf additional_results_delta01_b300.tar.gz \
  paper/output/additional_simulations/cluster_full_b300_delta01/point_merged_recovered \
  paper/output/additional_simulations/cluster_full_b300_delta01/bootstrap_merged \
  paper/output/additional_simulations/cluster_full_b300_delta01/sparse_cholesky_recovery
```

如果这些结果目录已完整存在，可跳过提取。

## 2. 一条命令重建表格

保持在仓库根目录执行：

```bash
Rscript --vanilla paper/scripts/30_make_additional_simulation_tables.R --overwrite=true
```

`--overwrite=true` 允许重新生成已有的表格文件。合并结果只读。输出在 `paper/output/additional_simulations/manuscript_tables/`。

在 cluster 使用同一脚本时，先进入 `/home/peng.1276/MRrr-additional-cluster`；已有合并结果可直接作为输入。

可显式指定路径：

```bash
Rscript --vanilla paper/scripts/30_make_additional_simulation_tables.R \
  --input-dir=paper/output/additional_simulations/cluster_full_b300_delta01 \
  --output-dir=paper/output/additional_simulations/manuscript_tables \
  --overwrite=true
```

路径包含空格时，引用整个参数，例如 `"--input-dir=D:/Research/MR results/run"`。

脚本读取：

| 输入（相对于 input-dir） | 用途 |
|---|---|
| `point_merged_recovered/additional_point_results.rds` | 已恢复的点估计、完整真值、逐重复偏差及恢复记录 |
| `bootstrap_merged/additional_bootstrap_results.rds` | Bootstrap SE、置信区间、覆盖指标及配置 |
| `bootstrap_merged/additional_bootstrap_summary.csv` | 独立核对现有 SE/CP 汇总 |

脚本拒绝不完整、非有限、配置未锁定或 N/B 不符的输入。其默认输入是恢复后的点估计文件。

## 3. 输出及 LaTeX 使用

| 输出 | 内容 |
|---|---|
| `table_rank_misspecification.tex` | 真秩 2，工作秩 1/2/3；四个 setting，共 24 行 |
| `table_approximate_low_rank.tex` | 真奇异值 (1,1,0.1)，降秩方法工作秩 2；共 28 行 |
| `additional_tables_summary.csv` | 52 行未主动舍入的数值汇总，方便再排版 |
| `additional_tables_entrywise.csv` | 52×27=1,404 行逐系数统计 |
| `additional_tables_preview.tex` | 可独立编译的两表文档 |
| `additional_tables_preview.pdf` | 本交付包中已编译好的两页预览 |
| `table_validation.txt` | 本次制表检查结果 |
| `table_provenance.rds` | 输入 MD5、脚本 MD5、R 环境、原结果元数据及恢复记录 |

将两份 `table_*.tex` 上传至 Overleaf，在论文导言区加入：

```latex
\usepackage{booktabs}
```

在正文所需位置加入：

```latex
\input{table_rank_misspecification.tex}
\input{table_approximate_low_rank.tex}
```

表号由论文自动编号；预览中的 Table 1/2 不是固定稿件表号。独立预览使用 letterpaper、10pt 和 0.7in 页边距。不同期刊模板可能需要调整表格字号或宽度。

R 脚本生成 LaTeX 和 CSV。若修改输入后还需更新 PDF，在有 LaTeX 的环境执行：

```bash
cd paper/output/additional_simulations/manuscript_tables
pdflatex -interaction=nonstopmode -halt-on-error additional_tables_preview.tex
```

也可把预览 `.tex` 与两份表格 `.tex` 一起上传至 Overleaf 编译。

## 4. 统计定义

每个 setting 和方法先对 C 的每个元素 j 计算 1,000 次重复的统计量，再在全部 27 个元素之间取中位数、25% 和 75% 分位数（R `quantile(type=7)`）。每个表格单元格第一行为中位数，第二行为 (Q1, Q3)。

| 指标 | 每个元素的定义 | 显示精度 |
|---|---|---|
| Bias | `abs(mean(C_hat_j - C_true_j))` | 3 位小数 |
| SD | 1,000 个点估计的样本标准差，分母 N−1 | 3 位小数 |
| SE | 1,000 次重复各自 bootstrap SE 的均值 | 3 位小数 |
| CP | 95% percentile bootstrap 区间覆盖真值的比例 ×100 | 1 位小数，百分数 |

Bias 是平均误差的绝对值，不是平均绝对误差。Approximate-low-rank 的 Bias 和 CP 使用包含第三奇异分量的完整真值。Sparse MR-rr 没有 bootstrap 推断，SE/CP 在 LaTeX 中显示为破折号，CSV 中留空。

## 5. 本次数据与验证记录

- 输入压缩包与 cluster 给出的 SHA-256 一致。
- 240 个 point chunk 已通过原合并校验，bootstrap 合并覆盖重复 1–1,000，B=300。
- 本脚本逐项核对点估计与已存偏差、完整真值、重复编号、点估计/bootstrap 数据种子、standard/MrDAG 的共同重采样种子。
- 对所有 36 行有 bootstrap 的方法，重新计算的 SE/CP 与合并 RDS 和原 CSV 一致；覆盖指标也与保存区间是否包含真值逐项一致。
- 点估计使用 recovery v2 的合并结果，保留 setting 1、replicate 840、Sparse 工作秩 1/2/3 的三个恢复拟合。恢复审计随表格 provenance 保存；其记录表明此前成功的估计保持不变。
- 本次运行使用 WebR 中的 R 4.6.0 读取 cluster R 4.5.1 生成的真实 RDS；制表脚本只依赖 base/recommended R 功能。PDF 经 LaTeX 编译并逐页检查。

成功时应出现 `ADDITIONAL MANUSCRIPT TABLES: PASS`、`Input files unchanged: PASS`，并给出 24/28/1,404 的行数。原始 Slurm 的失败任务记录仍是历史记录；本次制表使用已通过恢复合并校验的结果。

## 6. 把生成方法记录到 Git

在仓库根目录查看修改后，提交新增脚本和本文档：

```bash
git status --short
git add paper/scripts/30_make_additional_simulation_tables.R paper/ADDITIONAL_TABLES.md
git diff --cached --check
git commit -m "Add reproducible tables for robustness simulations"
git push
```

结果和表格位于已忽略的 `paper/output/`。不要用 `git add -f` 强制提交 RDS、压缩包或整批结果。要修改表格排版，修改生成脚本中的 `additional_tables_latex()` 后重新生成；保留从 RDS 到 LaTeX 的可重复流程。
