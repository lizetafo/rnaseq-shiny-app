# Example data

Four ready-to-use datasets. Each folder has a raw count matrix (`counts.csv`, genes × samples) and a sample sheet (`metadata.csv`).

| Folder | Organism | Gene IDs | Samples | Comparison | App settings |
|---|---|---|---|---|---|
| `prostate_hypoxia/` | Human | Ensembl (`ENSG…`) | 8 (4 + 4) | PC3 and LNCaP prostate cancer cells, **Hypoxia vs Normoxia** | Human · ENSEMBL · design `~Condition` · reference `Normoxia` · contrast `Condition, Hypoxia, Normoxia` |
| `tumor_vs_normal/` | Human | Ensembl | 37 (18 cancer + 19 normal) | **cancer vs normal** tissue | Human · ENSEMBL · `~type` · reference `normal` · contrast `type, cancer, normal` |
| `cd4_tcell_activation/` | Human | Entrez (numbers) | 16 (8 + 8) | CD4⁺ T cells, **Activated vs Resting** (also uninfected vs latent HIV infection) | Human · ENTREZID · `~Activation` · reference `Resting` · contrast `Activation, Activated, Resting` |
| `mouse_knockout/` | Mouse | Ensembl (`ENSMUSG…`) | 6 (3 + 3) | Mouse embryonic fibroblasts, **KO vs WT** | Mouse · ENSEMBL · `~genotype` · reference `WT` · contrast `genotype, KO, WT` |

`cd4_tcell_activation/metadata.csv` also has an `Infection` column, so you can try a two-factor design such as `~Activation + Infection`.

**Note:** in `cd4_tcell_activation/metadata.csv`, four samples were originally labeled `Activated ` (with a trailing space). R treats that as a separate group, so the label was corrected to `Activated`. Check your own sample sheets for stray spaces in group names.

## Input format

**Counts:** first column = gene IDs; one column per sample; raw integer counts (for example from featureCounts or HTSeq). Do not use normalized values such as TPM or FPKM.
```
GeneID,Sample1,Sample2,Sample3
ENSG00000000003,348,609,693
```

**Metadata:** first column = sample names, matching the count column names **exactly and in the same order**; then one column per experimental factor.
```
ID,Condition
Sample1,Normoxia
Sample2,Hypoxia
```
CSV, TSV and TXT files are all accepted.
