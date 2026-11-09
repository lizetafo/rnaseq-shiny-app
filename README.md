# Bulk RNA-seq Differential Expression Shiny App

An interactive **R Shiny** app that runs a complete bulk RNA-seq differential expression analysis with **DESeq2**, from a raw count matrix to annotated results and publication-ready plots, **without writing any code**.

Upload a count matrix and a sample sheet, check the experimental design the app detects, and click **Run**. The app runs DESeq2, adds gene names, and produces volcano plots, heatmaps, bar plots and result tables you can download.

![Volcano plot](docs/screenshots/2_volcano_plot.png)

---

## ✨ Features
- **Flexible input:** count matrices and sample sheets as CSV, TSV or TXT, up to 1 GB
- **Built-in checks:** missing or mismatched samples, sample order, duplicated samples and gene IDs, missing counts
- **Automatic design detection:** suggests the grouping factor, reference group, design formula and contrast, all editable
- **Multi-factor designs**, for example `~Activation + Infection`
- **Gene annotation** for human, mouse, rat, fly, zebrafish and yeast, or any other Bioconductor organism package, from Ensembl, Entrez or symbol IDs
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
4. Set the **reference level** to `Normoxia` and the **contrast** to `Condition, Hypoxia, Normoxia`, then click **Run**.

Step-by-step instructions with screenshots: **[docs/user_guide.md](docs/user_guide.md)**

## 🧪 Example datasets
| Dataset | Organism | Comparison |
|---|---|---|
| `prostate_hypoxia` | Human | Prostate cancer cell lines, hypoxia vs normoxia |
| `tumor_vs_normal` | Human | Tumor vs normal tissue (37 samples) |
| `cd4_tcell_activation` | Human | CD4⁺ T cells, activated vs resting (± latent HIV) |
| `mouse_knockout` | Mouse | Embryonic fibroblasts, knockout vs wild type |

Settings for each dataset and the required file format: [example_data/README.md](example_data/README.md)

## ✅ Validation
For all four datasets, the app's results were compared with DESeq2 run directly. Every gene, fold change, p-value and normalized count was **identical**, with the same significant genes (302, 6,834, 17,622 and 5,733). Details: [validation/VALIDATION.md](validation/VALIDATION.md).

## 📁 Repository structure
```
├── app.R                         # the Shiny app
├── setup/install_packages.R      # one-time package installation
├── docs/
│   ├── user_guide.md             # step-by-step instructions
│   ├── outputs.md                # what each table column and plot means
│   ├── methods.md                # analysis steps and where the app fits in an RNA-seq workflow
│   └── screenshots/
├── example_data/                 # 4 example datasets (counts + metadata)
└── validation/
    ├── VALIDATION.md
    └── compare_to_deseq2.R       # re-runs the comparison
```

## 🛠️ Built with
R · Shiny · DESeq2 · AnnotationDbi · ComplexHeatmap · ggplot2 · ggrepel · dplyr · DT

## 👩🏾‍🔬 Author
**Elizabeth Etafo**, Ph.D. candidate, Pharmaceutical Sciences, University of Nebraska Medical Center
[elizabethetafo.com](https://elizabethetafo.com) · [LinkedIn](https://www.linkedin.com/in/eetafo/)
