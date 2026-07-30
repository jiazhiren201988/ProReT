# ProReT 中文使用说明

ProReT（Program-space Reversal of Transcription）将疾病转录签名与 LINCS
药物扰动签名投影到同一套用户上传的 gene-by-program 基上，并依据
program 空间中的反向一致性对药物排序。

本 R 包只包含模型核心流程，不包含论文专用的 baseline、SOTA、RepoDB、
置换检验、bootstrap 或外部验证脚本。

## 用户需要准备

1. GEO Series Matrix：`*_series_matrix.txt` 或 `.txt.gz`；
2. gene-by-program 文本矩阵，第一列为基因，其余列为 program；
3. 明确的病例/对照分组规则、配对信息和必要协变量。

推荐先运行：

```r
geo <- read_geo_series("GSEXXXXX_series_matrix.txt.gz")
colnames(geo$phenotype)
head(geo$phenotype)
```

确认真实表型字段后，再指定 `group_column` 和互不重叠的正则表达式。

## 安装与 LINCS

依赖安装和完整示例见英文版 `README.md`。LINCS Level 5 GCTX 约 5.6 GB，
不会塞进 R 包；首次运行 `download_lincs()` 自动下载到用户缓存，后续复用。
如果服务器已有 LINCS 文件，可把三份清单文件放入同一目录，并通过
`lincs_cache` 指向该目录。

## 核心原则

- 疾病端使用 limma moderated t，而不是只保留显著 DEG；
- 疾病和药物必须使用完全相同的 basis 与 core genes；
- signed `gene_spectra_score` 的输出称为 signed program-template
  alignment，不称为 cNMF activity/usage；
- 默认使用 Gram 矩阵正则化后的 program 几何和 negative cosine；
- 同一药物×细胞的 LINCS 条件取中位数，再跨至少 3 个细胞系取均值；
- 疾病 core-gene 覆盖率低于 80% 时停止；
- 所有输入与关键设置写入 `run_manifest.rds`。

## 结果解释

药物排名表示“转录 program 逆转候选”，不等同于临床疗效预测。论文中的
gold standard 富集、SOTA 对比与独立验证应由独立复现 workflow 完成。
