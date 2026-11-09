# ================================================================
# One-time installation of every package the app uses.
# Run once in the RStudio Console, from the repository folder:
#   source("setup/install_packages.R")
# The app also installs anything missing on first launch, but doing it
# here first avoids a long wait inside the app.
# ================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c("shiny", "dplyr", "ggplot2", "ggrepel", "forcats",
                   "RColorBrewer", "circlize", "DT"))

# DESeq2, annotation and heatmap packages, plus the organism databases.
# Remove organisms you will never use to save download time.
BiocManager::install(c("DESeq2", "AnnotationDbi", "ComplexHeatmap",
                       "org.Hs.eg.db",   # human
                       "org.Mm.eg.db",   # mouse
                       "org.Rn.eg.db",   # rat
                       "org.Dm.eg.db",   # fruit fly
                       "org.Dr.eg.db",   # zebrafish
                       "org.Sc.sgd.db"), # yeast
                     ask = FALSE, update = FALSE)

cat("\nDone. Restart R (Session > Restart R), open app.R and click Run App.\n")
