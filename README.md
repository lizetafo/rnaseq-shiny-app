# Bulk RNA-seq Differential Expression Shiny App

An interactive **R Shiny** app that runs a complete bulk RNA-seq differential expression analysis with **DESeq2**, from a raw count matrix to annotated results and publication-ready plots, **without writing any code**.

Upload a count matrix and a sample sheet, check the experimental design the app detects, and click **Run**. The app runs DESeq2, adds gene names, and produces volcano plots, heatmaps, bar plots and result tables you can download.

### ▶️ [Try the app online](https://01a1c6c1-cbf8-a124-ddcd-1f24a83e1a02.share.connect.posit.cloud/)
An upgraded online version, with quality-control plots, an MA plot and a gene explorer, is in [RNAseq-web](https://github.com/lizetafo/RNAseq-web).

![Volcano plot](docs/screenshots/2_volcano_plot.png)

---

## ✨ Features
- **Flexible input:** count matrices and sample sheets as CSV, TSV or TXT, up to 1 GB
- **Built-in checks:** missing or mismatched samples, sample order, duplicated samples and gene IDs, missing counts
- **Guided design setup:** detects the grouping factor and design formula, asks you to choose the control group, then fills in the contrast, all editable
- **Multi-factor designs**, for example `~Activation + Infection`
- **Gene annotation** for human, mouse, rat, fly, zebrafish and yeast, or any other Bioconductor organism package, from Ensembl, Entrez or symbol IDs. Each organism database is downloaded only the first time it is used.
- **Plots:** volcano plot with labeled top genes, clustered heatmap (normalized, VST or rlog), bar plot, all downloadable
- **Tables:** annotated results, raw DESeq2 results and normalized counts, all downloadable
- **Validated:** results are identical to running DESeq2 directly ([validation](validation/VALIDATION.md))

## 🚀 Quick start
1. Install the packages once, in RStudio with this folder as the working directory:
   ```r
   source("setup/install_packages.R")
   ```
2. **Session → Restart R**, open `app.R`, and click **Run App**.
3. Upload `example_data/prostate_hypoxia/counts.csv` and `metadata.csv`, then click **Upload**.
4. Choose **Normoxia** as the **reference level (control group)**. The contrast fills in as `Condition, Hypoxia, Normoxia`. Click **Run**.

Step-by-step instructions with screenshots: **[docs/user_guide.md](docs/user_guide.md)**

## 🧪 Example datasets
| Dataset | Organism | Comparison |
|---|---|---|
| `prostate_hypoxia` | Human | Prostate cancer cell lines, hypoxia vs normoxia |
| `tumor_vs_normal` | Human | Colorectal tumor vs normal colon (37 samples) |

Both datasets are public GEO data, from Guo et al., *Nature Communications* 2019 and Kim et al., *Molecular Oncology* 2014 ([GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305) and [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760)). Sources, settings for each dataset and the required file format: [example_data/README.md](example_data/README.md)

## ✅ Validation
For both datasets, the app's results were compared with DESeq2 run directly. Every gene, fold change, p-value and normalized count was **identical**, with the same significant genes (302 and 6,834). Details: [validation/VALIDATION.md](validation/VALIDATION.md).

## 📁 Repository structure
```
├── app.R                         # the Shiny app
├── setup/install_packages.R      # one-time package installation
├── docs/
│   ├── user_guide.md             # step-by-step instructions
│   ├── outputs.md                # what each table column and plot means
│   ├── methods.md                # analysis steps and where the app fits in an RNA-seq workflow
│   └── screenshots/
├── example_data/                 # 2 example datasets (counts + metadata)
└── validation/
    ├── VALIDATION.md
    └── compare_to_deseq2.R       # re-runs the comparison
```

## 🛠️ Built with
R · Shiny · DESeq2 · AnnotationDbi · ComplexHeatmap · ggplot2 · ggrepel · dplyr · DT

## 👩🏾‍🔬 Author
**Elizabeth Etafo**, Ph.D. candidate, Pharmaceutical Sciences, University of Nebraska Medical Center
[elizabethetafo.com](https://elizabethetafo.com) · [LinkedIn](https://www.linkedin.com/in/eetafo/)
