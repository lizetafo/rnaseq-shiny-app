# Validation

**Question:** does the app give exactly the same results as running DESeq2 directly?

For each example dataset, [`compare_to_deseq2.R`](compare_to_deseq2.R) runs:
- **(a)** the app's own functions (`read_file`, `qc_check`, `prepare_meta`, `run_deseq`, `extract_results`), loaded from `app.R`, and
- **(b)** a plain DESeq2 analysis written independently: `DESeqDataSetFromMatrix` → filter genes with < 10 total counts → `DESeq` → `results(contrast, alpha = 0.05)`.

It then compares every result column (`baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`) and the normalized counts.

Re-run it in RStudio from the repository folder:
```r
source("validation/compare_to_deseq2.R")
```

## Results (R 4.3, DESeq2 1.42)

| Dataset | Comparison | Genes tested | Same genes | Largest difference: results / normalized counts | Significant genes (padj < 0.05): app / DESeq2 |
|---|---|---|---|---|---|
| prostate_hypoxia | Hypoxia vs Normoxia | 25,062 | ✅ | 0 / 0 | 302 / 302 |
| tumor_vs_normal | cancer vs normal | 28,225 | ✅ | 0 / 0 | 6,834 / 6,834 |
| cd4_tcell_activation | Activated vs Resting | 25,831 | ✅ | 0 / 0 | 17,622 / 17,622 |
| mouse_knockout | KO vs WT | 16,565 | ✅ | 0 / 0 | 5,733 / 5,733 |

The app's results are **identical** to plain DESeq2 for all four datasets. The app adds checks, annotation and plots without changing the statistics.

## Full app run
The app was also run in a browser on the three human datasets: upload → automatic design detection → Run → all plots → download of all three tables. The downloaded CSV files matched plain DESeq2 to within 10⁻¹⁴, which is rounding from saving the CSV. Gene symbols were assigned to 75–76% of Ensembl IDs and 96% of Entrez IDs. Ensembl IDs without a symbol are mostly non-coding or novel genes.
