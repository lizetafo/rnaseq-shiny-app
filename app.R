##################################################################################
############################ FINAL PROJECT #######################################
############################ ELIZABETH ETAFO #####################################
############################### 12/19/25 #########################################
##################################################################################

library(shiny)   #load shiny package

#####################################################################
######################Increasing Shiny Upload Limit##################
######################################################################

#default r shiny upload limit is 5mb. Some sequencing count files could be above this limit, 
#hence it is essential to increase the upload size limit. The below code increases 
#the upload limit to 1gb.

options(shiny.maxRequestSize=1000*1024^2)

##########################################################################
#####################Installing necessary packages#########################
##########################################################################

loadlib <- function(extra_Bioc=NULL, other_extra=NULL){  #incase you need to load extra packages
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager") 
  }
  Bioc_<- c("DESeq2", "AnnotationDbi", 
            "org.Hs.eg.db","org.Mm.eg.db",
            "org.Rn.eg.db", "org.Dm.eg.db",  
            "org.Dr.eg.db","org.Sc.sgd.db", "ComplexHeatmap")
  Bioc_pack<- unique(c(extra_Bioc, Bioc_)) #concatenates extra biocmanager packages to default pacakages
  Bioc_missing<- Bioc_pack[!vapply(Bioc_pack, requireNamespace, 
                                   FUN.VALUE = logical(1), quietly=TRUE)]
  if(length(Bioc_missing)>0){
    BiocManager::install(Bioc_missing, ask = FALSE, update=FALSE) 
  }
  
  other_<- c("dplyr", "ggplot2", "ggrepel", "forcats", "grid", "RColorBrewer", "circlize", "DT")
  other_pack<- unique(c(other_extra, other_))#concatenates extra cran packages to default cran pacakages
  other_missing<- other_pack[!vapply(other_pack, requireNamespace, 
                                     FUN.VALUE = logical(1), quietly=TRUE)]
  if(length(other_missing)>0){
    install.packages(other_missing, ask = FALSE, update=FALSE) 
  }
  
  suppressPackageStartupMessages({
    lapply(c(Bioc_pack, other_pack), library, character.only= TRUE)
  })
  
  message("All required packages have been successfully installed and loaded")
}

########################################################################
#####################Function to read count matrix######################
#######################################################################
read_file<- function(cnts, met, na_to_zero = FALSE){
  #reading count matrix file
  #example use: read_file("counts.csv", "metadata.csv")
  #some feature counts can be saved as either .csv, .txt
  if(grepl("\\.csv$", cnts, ignore.case = TRUE)) {   #if cnts cends with .csv
    counts<- read.csv(cnts, check.names = FALSE, as.is = TRUE)  
  } else {    #else read .txt or .tsv 
    counts<- read.delim(cnts, check.names = FALSE, as.is = TRUE)
  }
  rownames(counts)<- counts[[1]]  #GeneIds are typically found on the first column
  counts<- counts[, -1, drop=FALSE] #removing the column of gene ID
  counts<- data.matrix(counts) #making the columns numeric
  storage.mode(counts)<- "integer"   #ensuring the values are integers
  
  #handles NAs
  if(anyNA(counts)){
    if(na_to_zero){
      counts[is.na(counts)]<- 0
    }else{
      stop("NA values in counts, set na_to_zero = TRUE")
    }
    
  }
  
  #handles duplicate sample names/ gene IDs
  if (anyDuplicated(colnames(counts))){
    stop("Duplicated sample found in count data")
  }
  if (anyDuplicated(rownames(counts))){
    dup<- sum(duplicated(rownames(counts)))
    message("Removing", dup,"duplicated genes IDs")
    counts<- counts[!duplicated(rownames(counts)),]
  }
  
  #reading metadata
  #similar to feature counts, metadata may be .csv, .txt
  if(grepl("\\.csv$", met, ignore.case = TRUE)){
    meta<- read.csv(met, check.names = FALSE, as.is=TRUE)
  } else{
    meta<- read.delim(met, check.names = FALSE, as.is=TRUE)
  }
  meta[[1]]<- as.character(meta[[1]])
  if (anyDuplicated(meta[[1]])){
    meta[[1]]<- make.unique(meta[[1]])
  }
  rownames(meta)<- trimws(meta[,1]) #trimws removes whitespace 
  meta<- meta[, -1, drop=FALSE]
  
  invisible(list(counts = counts, metadata = meta))
  
}

########################################################################
################Function to quality check data##########################
#######################################################################
qc_check<- function(counts, meta){
  if(!all(colnames(counts)%in% rownames(meta)))
    stop("Some samples in counts data missing from metadata")
  if(!all(colnames(counts) == rownames(meta)))
    stop("Mismatch in sample order between counts and metadata")
  if(ncol(counts) != nrow(meta))
    stop("Difference in number of samples in counts and metadata")
  message("Quality checks passed")
}

#######################################################################
#######################Function to prepare metadata####################
#######################################################################
prepare_meta<- function(meta, Exp_groups=NULL, ref_levels = list()){
  #meta: metadata dataframe
  #Exp_groups: character vector of columns to use as experimental factors, e.g. c("Activation", "Infection")
  #ref_levels: list of reference levels, e.g.list(Infection = "Uninfected", Activation = "Resting")
  
  #Some quality check of metadata 
  if(is.null(Exp_groups)){ #warning message if Exp_group is not provided
    stop("Provide experimental grouping column names (Exp_groups).")
  }
  missing_cols<- setdiff(Exp_groups, names(meta))
  if (length(missing_cols) > 0){  #warning if experiment groups are not found in metadata
    stop(sprintf("These Exp_groups are missing from metadata: %s", paste(missing_cols, collapse = ",")))
  }
  
  #convert all experiment groups columns to factors
  for(col in Exp_groups){
    meta[[col]]<- as.factor(meta[[col]])
  }
  
  #set reference levels
  
  for (nm in names(ref_levels)){
    if (nm %in% names(meta)){
      meta[[nm]]<- relevel(meta[[nm]], ref=ref_levels[[nm]])
    }else{
      warning(sprintf("Column '%s' not found in metadata, skipping reference level.", nm))
    }
  }
  invisible(meta)
}

#######################################################################
#########################Function to run DESeq2#######################
#######################################################################
run_deseq<- function(counts, meta, Exp_design=NULL, min_count=10){
  #counts and meta defined from previous function
  #Exp_design: to be used later in function (e.g. Exp_design= ~type)
  
  #Check if DESeq2 is installed and load the package
  if(!requireNamespace("DESeq2", quietly = TRUE)){
    if(!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager") 
    BiocManager::install("DESeq2", ask = FALSE, update=FALSE)
  }
  suppressPackageStartupMessages(library(DESeq2))
  
  #validate function
  if(is.null(Exp_design))
    stop("Provide a valid DESeq2 design formula (Exp_design)")
  if(!inherits(Exp_design, "formula"))
    stop("Exp_design must be a formula object, e.g. ~Condition or ~Activation + Infection")
  
  #Create DESeq2 dataset
  dds<- DESeqDataSetFromMatrix(countData = counts, 
                               colData = meta,
                               design = Exp_design )
  #Prefilter low count genes
  keep<- rowSums(counts(dds)) >= min_count
  dds <- dds[keep,]
  
  #Run DESeq normalization and model fitting
  deg <- DESeq(dds)
  
  message("DESeq2 analysis completed")
  deg
}

#########################################################################
##############Function to extract results and normalization##############
#########################################################################
extract_results<- function(deg, contrast=NULL, alpha=NULL, normalized= TRUE){
  #Example use: contrast= c("type", "normal", "cancer")----in metadata we can 2 controls in the "type" column, normal and cancer
  #some basic checks
  if(!inherits(deg, "DESeqDataSet"))
    stop("'deg' must be a DESeqDataSet (output of DESeq())")
  if (is.null(alpha)) alpha <- 0.05
  if (!is.numeric(alpha) || length(alpha) !=1 || alpha <=0 || alpha >= 1)
    stop("'alpha' must be a single number in (0, 1)")
  
  #result with and without contrast
  if(is.null(contrast)){
    res<- results(deg, alpha = alpha)
  } else if (is.character(contrast) && length(contrast) ==1){
    res <-results(deg, name = contrast, alpha = alpha)
  } else{
    res<- results(deg, contrast = contrast, alpha= alpha)
  }
  
  #Normalized counts matrix
  norm<- counts(deg, normalized= normalized)
  
  list(count_results= res, normalized_counts= norm)
}

########################################################################
#################Function to annotate gene names#######################
########################################################################
annotate_results<- function(res, orgdb, id_type=NULL, 
                            to_symbol="SYMBOL"){
  
  #basic checks
  if(missing(orgdb)) stop("Please pass an OrgDb object, e.g. org.Hs.eg.db or org.Mm.eg.db")
  if(!inherits(orgdb, "OrgDb")) stop("'orgdb' must be an OrgDb (e.g. org.Hs.eg.db) ")
  
  if(!requireNamespace("AnnotationDbi", quietly = TRUE)) 
    stop("Install package 'AnnotationDbi'. Install via BiocManager::install('AnnotationDbi').")
  
  resdf<- as.data.frame(res)
  genekeys<- rownames(resdf)
  if(is.null(genekeys)) stop("Row names of 'res' must contain gene IDs") 
  
  
  #guess ID type if not provided
  guess_idtype <- function(ids){
    if (all(grepl("^ENSG\\d+$", sub("\\.\\d+$", "", ids)))) return ("ENSEMBL")
    if (all(grepl("^ENSMUSG\\d+$", sub("\\.\\d+$", "", ids)))) return("ENSEMBL")
    if (all(grepl("^\\d+$", ids))) return("ENTREZID")
    if (all(grepl("^[A-Za-z0-9._-]+$", ids))) return("SYMBOL")
    return(NA_character_)
  }
  if(is.null(id_type)){
    id_type<- guess_idtype(genekeys)
  }
  
  #validate ID types
  ktypes <- AnnotationDbi::keytypes(orgdb)
  if (is.na(id_type) || !(id_type%in% ktypes)) {
    stop(sprintf("Could not determine valid id_type. Please specify one explicitly.
                  Available types: %s", paste(ktypes, collapse= ",")))
  }
  
  #perform gene mapping  
  resdf$GeneName<- AnnotationDbi::mapIds(orgdb,
                                         keys= genekeys,
                                         keytype=id_type,
                                         column= to_symbol,
                                         multiVals = "first")
  
  attr(resdf, "id_type_used")<- id_type
  resdf
}

#######################################################################
######################Function to get top genes#######################
#######################################################################
top_genes<- function(resdf, padj="padj", n=15){ #can change n to number of best genes you want to extract
  if(!requireNamespace("dplyr", quietly= TRUE)) #warning message if you dont have package
    stop("Install 'dplyr' package first")
  
  if(!padj %in% names(resdf))  #warning message is padj is not in your result dataframe
    stop("Column 'padj' not found in results dataframe")
  
  best_genes<- dplyr::arrange(resdf, !!as.name(padj))
  head(best_genes, n)
}

#######################################################################
##################Function to generate volcano plot####################
########################################################################
plot_volcano<- function(resdf, best_genes= NULL, out_png= NULL, 
                        log2FoldChange="log2FoldChange",
                        padj= "padj", 
                        lfc_thr= 1, padj_thr= 0.05,
                        label_col= "GeneName", 
                        main= "Differential Expression Volcano Plot") {
  
  #Some checks
  if(!all(c(log2FoldChange, padj) %in% names(resdf))) #checks if columns required for plotting are in data frame
    stop("Missing some columns for volcano plot")
  if(!requireNamespace("ggplot2", quietly = TRUE)) #checks if ggplot2 is installed
    stop("'ggplot2' package required for plotting volcano plot")
  
  #removes NAs
  keep<- !is.na(resdf[[padj]]) & !is.na(resdf[[log2FoldChange]]) 
  vol<- resdf[keep, , drop=FALSE]
  
  # .sig marks significat genes (pvalue <= threshold and log2foldchange >= threshold)
  vol$neglog10<- -log10(vol[[padj]])
  
  vol$status<- "NS"
  vol$status [vol[[padj]]<= padj_thr&vol[[log2FoldChange]]>= lfc_thr]<- "Up"
  vol$status [vol[[padj]]<= padj_thr&vol[[log2FoldChange]]<= -lfc_thr]<- "Down"
  vol$status<- factor(vol$status, levels=c("Up","Down", "NS"))
  
  #building the ggplot using the ggplot2 package
  vol_plot<- ggplot2::ggplot(
    vol,
    ggplot2::aes(x=.data[[log2FoldChange]], y=.data[["neglog10"]])) + 
    ggplot2::geom_point(ggplot2::aes(color=.data[["status"]]), size= 1.3, alpha=0.7) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype= "dashed") +
    ggplot2::geom_hline(yintercept = -log10(padj_thr), linetype ="dashed") +
    ggplot2::scale_color_manual(
      values=c(Up="#F72585", Down="#3B4CC0", NS="grey50"),
      name="Differential Expression"
    )+
    ggplot2::labs(x="log2 fold change", y="-log10(adjusted p-value)", title= main)+
    ggplot2::theme_minimal(base_size = 14)+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face="bold"),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.title = ggplot2::element_text(face="bold", size = 12),
                   legend.text = ggplot2::element_text(size=10)
    )
  
  if(!is.null(best_genes) && nrow(best_genes)>0){
    lab_rows<- intersect(row.names(best_genes), row.names(vol))
    if (length(lab_rows)) {
      label_df <- vol[lab_rows, , drop = FALSE]
      # choose label text: use label_col if present, else fallback to rownames
      label_df$..label.. <- if (label_col %in% names(label_df)) label_df[[label_col]] else rownames(label_df)
      
      #adding labels in lot for top genes
      vol_plot<- vol_plot +
        ggrepel::geom_text_repel(data= label_df, 
                                 ggplot2::aes(x=.data[[log2FoldChange]],
                                              y=.data[["neglog10"]],
                                              label= .data[[label_col]]), 
                                 size=3,
                                 max.overlaps = 20,
                                 box.padding = 0.4,
                                 point.padding = 0.3,
                                 segment.color="grey30",
                                 force = 2,
                                 min.segment.length = 0)  
      
    }
  }
  if (!is.null(out_png)){
    ggplot2::ggsave(filename = out_png, 
                    plot=vol_plot, width = 7,
                    height = 6, dpi = 150)
  } else {
    print(vol_plot)
  }
  invisible(vol_plot)
}
########################################################################
#######################Function to plot heat map########################
########################################################################
plot_Heatmap<- function(deg, resdf, out_png= NULL, n=30, main="Top Differentially Expressed Genes", mat=NULL){
  #check for required package
  if(!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("'ComplexHeatmap' package required")
  if(!requireNamespace("circlize", quietly = TRUE))
    stop("'circlize' package required")
  if(!requireNamespace("grid", quietly = TRUE))
    stop("'grid' package required")
  
  use_genes<- top_genes(resdf, n=n) #extracts topgenes from the top genes function
  
  if(is.null(mat)){
    mat<- counts(deg, normalized=TRUE) #extract normalized counts from DESeq object 
  }
  
  genes<- intersect(rownames(use_genes), row.names(mat)) #only keep genes presne ton both top_genes and expression matrix
  mat<- mat[genes, , drop=FALSE]
  mat.z<- t(scale(t(mat))) #zscore normalize each gene across samples
  mat.z[is.na(mat.z)] <- 0 #repace NA with zero
  labs<- use_genes[genes, "GeneName"]
  labs <- ifelse(is.na(labs) | labs == "", genes, labs) 
  rownames(mat.z)<- make.unique(labs)
  
  design_var<- attr(terms(design(deg)), "term.labels")[1] #extracts first design variable from DESeq2 design
  
  group_vec<- NULL
  top_annot<- NULL
  my_cols <- c(    "#4D90BD","#E06C9F", "#8BC77B","#CFA6E8","#F2C14E")
  if(!is.null(design_var)){
    cd<- as.data.frame(SummarizedExperiment::colData(deg))
    cd<- cd[colnames(mat.z), , drop=FALSE]  #match sample metadata to mat.z column
    if(design_var%in% names(cd)){
      group_vec<- factor(cd[[design_var]])
      levels(group_vec) <- tools::toTitleCase(levels(group_vec)) #get groouping factor for samples
      cols<- setNames(my_cols[seq_len(nlevels(group_vec))],
                      levels(group_vec))
      
      #building the heatmap annotaion bar showing sample groups
      top_annot<- ComplexHeatmap::HeatmapAnnotation(
        Group=group_vec,
        col = list(Group=cols),
        annotation_legend_param = list(
          title="Groups",
          title_gp=grid::gpar(fontsize=12, fontface="bold"),
          labels_gp=grid::gpar(fontsize=10)
        )
      )
    }
  }
  #generating colo scale for zscores
  heat_colors<- circlize::colorRamp2(
    c(-3,0,3), c("#3B4CC0", "white","#F72585")
  )
  #generating heatmap using the complexheatmap package
  Heat_plot<- ComplexHeatmap::Heatmap(mat.z, 
                                      name = "Z-score",
                                      heatmap_legend_param = list(
                                        title_gp=grid::gpar(fontsize=12, fontface="bold"),
                                        labels_gp=grid::gpar(fontsize=10)
                                      ),
                                      col=heat_colors,
                                      top_annotation = top_annot,
                                      cluster_rows=TRUE, 
                                      cluster_columns=TRUE,
                                      column_split = group_vec,
                                      show_row_names = TRUE,
                                      show_column_names = FALSE,
                                      column_title = main,
                                      column_title_gp = grid::gpar(fontsize=16, fontface="bold"),
                                      rect_gp = grid::gpar(col="grey70", lwd=0.3))
  if (!is.null(out_png)) {
    png(out_png, width = 1200, height = 1200, res = 150)
    on.exit(dev.off(), add = TRUE)
  }
  suppressMessages(ComplexHeatmap::draw(Heat_plot))
}
######################################################################
####################Function for Barplot##############################
#####################################################################
plot_bar<- function(resdf, n=30,log2FoldChange="log2FoldChange",
                    padj= "padj",label_col="GeneName", x_limit=NULL,
                    direction= c("up", "down", "all"),
                    bar_color="#F72585", color_by_direction=TRUE,
                    up_color="#F72585",
                    down_color="#3B4CC0",
                    neutral_color="grey70",
                    title="Differential Expression Bar Plot",
                    out_png=NULL){
  
  #Some checks
  if(!all(c(log2FoldChange, padj) %in% names(resdf))) #checks if columns required for plotting are in data frame
    stop("Missing some columns for bar plot")
  if(!requireNamespace("ggplot2", quietly = TRUE)) #checks if ggplot2 is installed
    stop("'ggplot2' package required")
  direction<- match.arg(direction)
  #removes NAs
  keep<- !is.na(resdf[[padj]]) & !is.na(resdf[[log2FoldChange]]) 
  bardf<- resdf[keep, , drop=FALSE]
  
  use_genes<- top_genes(bardf, n=n)
  
  #specifies direction of the bar plot
  if(direction == "up") use_genes<- dplyr::filter(use_genes, log2FoldChange>0)
  if(direction == "down") use_genes<- dplyr::filter(use_genes, log2FoldChange<0)
  use_genes<- dplyr::filter(use_genes, !is.na(.data[[label_col]]))
  
  #reorder factor levels by Log2FC so bars are correctly sorted
  use_genes[[label_col]] <- forcats::fct_reorder(
    use_genes[[label_col]], use_genes[[log2FoldChange]]
  )
  
  if(color_by_direction&&direction=="all"){
    use_genes<- dplyr::mutate(
      use_genes,
      Direction=dplyr::case_when(
        .data[[log2FoldChange]]>0 ~"Up",
        .data[[log2FoldChange]]<0 ~"Down",
        TRUE  ~"No change"
      )
    )
    #building bar plot with color matching the direction
    barpl<- ggplot2::ggplot(
      use_genes,
      ggplot2::aes(x=.data[[label_col]],
                   y=.data[[log2FoldChange]],
                   fill=.data$Direction)
    )+
      ggplot2::geom_col(width=0.7)+
      ggplot2::scale_fill_manual(
        values = c(Up=up_color, Down=down_color, `No change`=neutral_color),
        name="Expression Change"
      )
  }else{
    barpl<- ggplot2::ggplot(
      use_genes,
      ggplot2::aes(x=.data[[label_col]],y=.data[[log2FoldChange]])
    )+
      ggplot2::geom_col(width=0.7, fill=bar_color)
  }
  #generates barplot
  barpl<- barpl+ 
    ggplot2::coord_flip()+
    ggplot2::geom_hline(yintercept = 0, color="grey30")+
    ggplot2::labs(x=NULL, y="log2(fold change)", title = title)+
    ggplot2::theme_classic(base_size = 12)+
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust=0.5, face="bold"),
      axis.ticks.y = ggplot2::element_line(),
      axis.line = ggplot2::element_line(color="black"),
      axis.text.y = ggplot2::element_text(size=8, face = "bold"),
      axis.text.x = ggplot2::element_text(size=8, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      legend.text=ggplot2::element_text(size = 8),
      legend.title=ggplot2::element_text(size=9,face = "bold")
    )
  if (!is.null(x_limit)) {
    barpl <- barpl + ggplot2::scale_y_continuous(limits = x_limit)
  }
  if (!is.null(out_png)) {
    ggplot2::ggsave(filename = out_png, plot = barpl, width = 6, height = 4, dpi = 300)
  }
  print(barpl)
  invisible(barpl)
}

#####################################################################
###########################Defining UI##############################
######################################################################
ui<- fluidPage(
  #App title:
  titlePanel("Bulk RNAseq Differential Gene Expression Pipeline"),
  
  tags$head(tags$style(HTML("
                             .mySidebar{
                             font-size: 14px; #font size of the side bar text
                            }
                            "))
  ),
#######################################################################  
#############Sidebar Layouts with Inputs and Output Definitions:#########
##########################################################################
  sidebarLayout(
    sidebarPanel(
      class= "mySidebar",
      
      fileInput(inputId = "Count",  #using fileinput to ulpad count and meta file
                "Upload Counts File", 
                accept = c(".csv", ".txt", ".tsv")),
      helpText("Upload counts file in one of the following formats: CSV, TSV, or TXT format."), #added help text as an explanation to user on what to do
      fileInput(inputId = "Metadata", 
                "Upload Metadata File", 
                accept = c(".csv", ".txt", ".tsv")),
      helpText("Upload metadata file in one of the following formats: CSV, TSV, or TXT format."),
      
      checkboxInput(inputId = "NA_to_zero", #using checkbox input for users to select if they want to remove NA, default is set to TRUE(remove NA)
                    "Convert NA in counts to zero",
                    value = TRUE),
      
      actionButton("upload_data", "Upload", class="btn-info"), #using action button to upload  data
      helpText("Upload data and automatically detect experimental design"),
      
      tags$hr(),
      
      conditionalPanel( #this will show additonal sidebar panel inputs only after user clicks on upload 
        condition = "input.upload_data>0",
      
      #organism database selection
      selectInput(inputId = "Orgdb", 
                  "Organism Database", 
                  choices = c("Human (org.Hs.eg.db)"= "human", 
                              "Mouse (org.Mm.eg.db)"="mouse",
                              "Rat (org.Rn.eg.db)" ="rat",
                              "Fly (org.Dm.eg.db)"= "fly",
                              "Zebrafish (org.Dr.eg.db)"="zfish",
                              "Yeast (org.Sc.sgd.db)"= "yeast",
                              "Other (not listed)"= "other")
      ),
      #Custom Organism database, incase user has other's that are not listed above
      conditionalPanel(
        condition = "input.Orgdb=='other'",
        textInput("OtherOrg", "Enter organism name or annotation package:")
      ),
      helpText("Choose appropriate organism reference database for gene annotation."),
      
      #gene ID type selction
      selectInput(inputId = "idtype", 
                  "ID Type for Gene Annotation", 
                  choices = c("ENSEMBL" ="ensembl",
                              "ENTREZID"= "entrezid",
                              "SYMBOL"= "symbol",
                              "Other (not listed)"= "other")
      ),
      #similarly, custom ID type can be entered if not provided in the select options
      conditionalPanel(
        condition = "input.idtype=='other'",
        textInput("OtherIDtype", "Enter custom ID type (e.g., REFSEQ, UNIPROT):")
      ),
      helpText("Choose appropriate ID Type for gene names annotation (example: ENSEMBL)."),
      
      #text input to enter experiment design types and conditions
      textInput(inputId = "Exp_groups", "Experimental Grouping Factor(s)"),
      helpText("Enter experimental grouping factor(s) (comma-seperated) (example: Activation, Infection)"),
      
      textInput(inputId = "ref_col", "Reference Column"),
      helpText("Enter reference column name from metadata (example: Activation)"),
      textInput(inputId = "ref_levels", "Reference level"),
      helpText("Enter reference level from metadata (example: Resting)"),
      
      textInput(inputId = "Exp_design", "Experimental Design Formula"),
      helpText("Enter experiment design formula to be used for DESeq2 (example ~Activation)"),
      
      #text input to select for contrast
      textInput(inputId = "contrast", "Contrast (Factor, Case, Control)"),
      helpText("Enter contrast (from metadata)( example: Activation, Activated, Resting)"),
      
      #user can select thresholds
      numericInput("alpha", "Adjusted p-value Threshold (alpha)", value = 0.05, 
                   min = 0.0001, max= 0.5, step = 0.01),
      numericInput("lfc_thr", "log2FC Threshold", value=1, step = 0.1),
      numericInput("padj_thr_plot", "Adjusted p-value Threshold", value=0.05, min=1e-6, max=0.5),
      
      #Options for heatmap transformation types
      selectInput("heatmap_transform", "Transformation for Heatmap",
                  choices = c("Normalized counts"="norm",
                              "VST"="vst",
                              "rlog"="rlog"),
                  selected = "norm"),
      #User can determine the number of tope genes they want displayed in their volcano plot and heatmap
      numericInput("n", "Number of Top genes",
                   value = 30, min=10, max = 200, step=5),
      helpText("Enter number of top genes to use for plots"),
      
      #action button to run DESeq2 pipeline
      actionButton("run", "Run", class= "btn-primary"),
      helpText("Run RNA DESeq Pipeline")
      
    )
    ),
########################################################################################
###############################Main Panel##########################################
##################################################################################
    mainPanel(
      tabsetPanel(
        #Graphical plots tab
        tabPanel("Graphical Plots",
                 br(),
                 h4("Differential Gene Expression Plots"),
                 tabsetPanel(
                   #tab for volcano plot
                   tabPanel(
                     "Volcano Plot",
                     br(),
                     fluidRow(
                       column(
                         width = 8,
                         plotOutput("volcano_plot", height = "1000px"),
                         downloadButton("download_volcano", "Download Volcano Plot")
                       ),
                       column(
                         width = 4,
                         conditionalPanel(
                           condition="input.run>0",
                           h4("Top 15 genes (by padj)"),
                           DT::DTOutput("volcano_top15")
                         )
                       )
                     )
                   ),
                   #tab for heatmap
                   tabPanel("Heatmap Plot",
                            br(),
                            plotOutput("heatmap", height = "1000px"),
                            downloadButton("download_heatmap", "Download Heatmap")),
                   #tab for bar plot
                   tabPanel("Bar Plot",
                            br(),
                            plotOutput("bar_plot", height = "1000px"),
                            downloadButton("download_bar", "Download Bar Plot")
                   )
                 )
        ),
        
        #Resutls table tab
        tabPanel("Results",
                 br(),
                 h4("Differential Gene Expression Tabular Results"),
                 tabsetPanel(
                   tabPanel("Annotated",
                            br(),
                            DT::DTOutput("results_ann")),
                   tabPanel("Raw DESeq2 Results",
                            br(),
                            DT::DTOutput("results_raw")),
                   tabPanel("Normalized Counts",
                            br(),
                            DT::DTOutput("results_norm"))
                 ),
                 br(),
                 downloadButton("download_ann", "Download Annotated Results"),
                 downloadButton("download_counts", "Download DESeq2 Results"),
                 downloadButton("download_norm", "Download Normalized Counts")
        )
      )
      
    )
    
    
  )
  
)



########################################################################
###########################Defining Server##############################
########################################################################

server<- function(input, output, session){
  loadlib()  #load required libraries
  
  #makes a reactive to sotre uploaded files
  data_store<- reactiveValues(
    dat=NULL
  )
  
  #Ensure that files are uploaded before user can click upload button
  observeEvent(input$upload_data,{
    validate(
      need(input$Count, "Please upload count matrix file"),
      need(input$Metadata, "Please upload metadata file")
    )
    
    #reads count and metadata
    dat<- read_file(
      cnts= input$Count$datapath,
      met= input$Metadata$datapath,
      na_to_zero = input$NA_to_zero)
    
    #basic quality checks
    qc_check(dat$counts, dat$metadata)
    
    
    #uploaded data are stored in the reactive container initially created to be unsed in pipeline
    data_store$dat<- dat
    
    #Autodetect experimental desing from metadata (User has option to change these)
     meta<- dat$metadata
      factor_cols<- names(meta)
      
      candidate<-NULL
      for (col in factor_cols){
        u<- unique(meta[[col]])
        if (length(u)>= 2 && length(u)<=10){
          candidate <- col
          break
        }
      }
      if (!is.null(candidate)){
        f<- factor(meta[[candidate]])
        lv<- levels(f)
        
          ref_level<- lv[1]    #Choose first level as the reference level
          case_level<- lv[length(lv)]    #Choose last level as the case
          
          #Prefils the UI fiels with what was guessed, may not be want the user wants so user can change to desired settings
          updateTextInput(session, "Exp_groups", value= candidate)
          updateTextInput(session, "ref_col", value= candidate)
          updateTextInput(session, "ref_levels", value= ref_level)
          updateTextInput(session, "Exp_design", value= paste0("~", candidate))
          updateTextInput(session, "contrast", 
                      value= paste(candidate, case_level, ref_level, sep=","))
    
        
      }
      
      showNotification("Data uploaded and design guessed. You can review/edit settings and then click Run.", type="message")
  })
   
  
  ###building ref levels from UI
  ref_levels<- reactive({
    ref_list<- list()
    if(nzchar(input$ref_col)&&nzchar(input$ref_levels)){
      ref_list[[input$ref_col]] <- input$ref_levels
    }
    ref_list
  })
  
  ###chooisng organism data base
  orgdb <- reactive({
    if(input$Orgdb != "other"){
      db<- switch(input$Orgdb,
                  human = org.Hs.eg.db,
                  mouse = org.Mm.eg.db,
                  rat = org.Rn.eg.db,
                  fly = org.Dm.eg.db,
                  zfish = org.Dr.eg.db,
                  yeast = org.Sc.sgd.db      
      )
      return(db)
    }
    
    #if the user selected "other", the user is required to provide the organism annotation package name
    validate(need(nzchar(input$OtherOrg),
                  "Please enter the organism annotation package name.")) 
    pkg_name<- input$OtherOrg
    
    #informs user if the package is not installed
    if(!requireNamespace(pkg_name, quietly=TRUE)){
      validate(
        need(FALSE,
             paste0("Package '", pkg_name,
                    "'is not installed. Install it with BiocManager::install('", pkg_name, "')"))
      )
    }
    
    #loads package and return the annotation object
    ns<- loadNamespace(pkg_name)
    get(pkg_name, envir = ns)
  })
  
  #choose gene ID type
  id_type <- reactive({
    if(input$idtype != "other"){
      id<- switch(input$idtype,
                  ensembl = "ENSEMBL",
                  entrezid= "ENTREZID",
                  symbol= "SYMBOL"
      )
      return(id)
    }
    
    #similarly, user need to enter the custom ID type if the select other
    validate(need(nzchar(input$OtherIDtype),
                  "Please enter custom ID type.")) 
    toupper(trimws(input$OtherIDtype))
  })
  
  
  ##reactive to run RNA pipeline
  pipeline<- eventReactive(input$run, {
    validate(   #ensures count and metadata have been uploaded
      need(!is.null(data_store$dat),
           "Please upload data and detect design first (click 'Upload')")
    )
    
    #Retrieve stored data
    dat<- data_store$dat
    
    
    #Prepare experiemtnal groping vector 
    exp_groups<- strsplit(input$Exp_groups, ",")[[1]]
    exp_groups<- trimws(exp_groups)
    exp_groups<- exp_groups[nzchar(exp_groups)]
    
    #prepare metadata
    meta<-prepare_meta(dat$metadata, Exp_groups = exp_groups, ref_levels = ref_levels())
    
    #give error if desing formula is invalid
    Design_formula<- tryCatch(as.formula(input$Exp_design),
                              error=function(e){
                                stop("Design formula is not valid. Example: ~Activation or ~Activation + Infection")
                              })
    
    #Run DESeq2 pipeline
    dseq<- run_deseq(dat$counts, meta, Exp_design = Design_formula)
    
    #tries to compute VST and rlog transforms
    vsd=NULL
    rld=NULL
    try({
      vsd<- DESeq2::vst(dseq, blind = TRUE)
      }, silent = TRUE)
    
    try({
      rld<- DESeq2::rlog(dseq, blind = TRUE)
    }, silent = TRUE)
    
    
    #Building Contrast Vector, gives error if contrast length is equal to 3
    contrast_vec<- strsplit(input$contrast, ",")[[1]]
    contrast_vec<- trimws(contrast_vec)
    if(length(contrast_vec) !=3){
      stop("Invalid contrast entered")
    }
    
    #extract DESeq2 results
    res_list<- extract_results(dseq,  contrast = contrast_vec, alpha = input$alpha,  normalized = TRUE)
    
    #orgdb selection
    db<- orgdb()
    id<- id_type()
    
    ##annoate gene names
    
    ann<- annotate_results(res_list$count_results, orgdb = db, id_type = id, to_symbol = "SYMBOL")
    
    ##top genese
    
    top__genes<- top_genes(ann, padj= "padj", n= input$n)
    
    list(dat=dat,
         meta= meta,
         dseq=dseq,
         vsd= vsd,
         rld= rld,
         res_list= res_list,
         ann= ann,
         top__genes= top__genes)
    
  })
  
  #### PLOTS OUTPUTS
  ###Rendering plots and including download options
  output$volcano_plot<- renderPlot({
    p<- pipeline()
    plot_volcano(
      p$ann, p$top__genes, out_png = NULL,
      lfc_thr = input$lfc_thr,
      padj_thr = input$padj_thr_plot
    )
  })  
  
  output$volcano_top15<- DT::renderDT({
    req(input$run>0)
    p<- pipeline()
    df<- p$ann
    keep<- !is.na(df$padj)& !is.na(df$log2FoldChange) &
      df$padj<= input$padj_thr_plot &
      abs(df$log2FoldChange)>= input$lfc_thr
    df<- df[keep, , drop=FALSE]
    top15<- top_genes(df, padj="padj", n=15 )
    top15[, c("GeneName", "log2FoldChange", "padj"), drop=FALSE]
    }, options= list(pageLength=15, searching=FALSE, ordering=TRUE))
  

  
  output$download_volcano<- downloadHandler(
    filename = function() "Volcano_plot.png",
    content = function(file){
      p<- pipeline()
      plot_volcano(p$ann, p$top__genes, out_png=file, 
                   lfc_thr =input$lfc_thr, padj_thr = input$padj_thr_plot )
    }
  )
  output$heatmap <- renderPlot({
    p<- pipeline()
    
    mat<- NULL
    if(input$heatmap_transform== "vst" && !is.null(p$vsd)){
      mat<- SummarizedExperiment::assay(p$vsd)
    }else if(input$heatmap_transform=="rlog"&& !is.null(p$rld)){
        mat<- SummarizedExperiment::assay(p$rld)
    }else {
        mat<- p$res_list$normalized_counts
      }
    
    plot_Heatmap(p$dseq,p$ann,out_png = NULL, n=input$n, main = paste("Top DE Genes (", input$heatmap_transform,")", sep=""),
                 mat=mat)
  })
  
  
  output$download_heatmap<- downloadHandler(
    filename = function() "Heatmap.png",
    content = function(file){
      p<- pipeline()
      mat<- NULL
      if(input$heatmap_transform=="vst" && !is.null(p$vsd)){
        mat<- SummarizedExperiment::assay(p$vsd)
      }else if (input$heatmap_transform =="rlog" && !is.null(p$rld)){
        mat<- SummarizedExperiment::assay(p$rld)
      }else{
        mat<- p$res_list$normalized_counts
      }
      plot_Heatmap(p$dseq, p$ann, out_png=file, n=input$n, 
                   main = paste("Top DE Genes (", input$heatmap_transform,")", sep=""),
                   mat=mat)
    }
  )
  output$bar_plot <- renderPlot({
    p<- pipeline()
    plot_bar(p$ann, n=input$n, direction="all",out_png = NULL)
  })
  output$download_bar<- downloadHandler(
    filename = function() "Barplot.png",
    content = function(file){
      p<- pipeline()
      plot_bar(p$ann, direction="all", out_png=file, n=input$n)
    }
  )
  
  
  #####Tables and downloads#########
  
  #annotated results preview
  output$results_ann<- DT::renderDT({
    p<- pipeline()
    head(p$ann, 100)
  }, options =list(pageLength=25, scrollX=TRUE))
  #raw DESeq2 results preview
  output$results_raw<- DT::renderDT({
    p<- pipeline()
    head(as.data.frame(p$res_list$count_results), 100)
  }, options =list(pageLength=25, scrollX=TRUE))
  
  #normalized cunt preview
  output$results_norm<- DT::renderDT({
    p<- pipeline()
    head(as.data.frame(p$res_list$normalized_counts), 100)
  }, options =list(pageLength=25, scrollX=TRUE))
  
  #download annotated results
  output$download_ann<- downloadHandler(
    filename = function() "Annotated_results.csv",
    content = function(file){
      p<- pipeline()
      write.csv(p$ann, file)
    }
  )
  
  #download DESeq2 results
  output$download_counts<- downloadHandler(
    filename = function() "DESeq2_results.csv",
    content = function(file){
      p<- pipeline()
      write.csv(as.data.frame(p$res_list$count_results), file)
    }
  )
  
  #download norm counts
  output$download_norm<- downloadHandler(
    filename = function() "Normalized_counts.csv",
    content = function(file){
      p<- pipeline()
      write.csv(as.data.frame(p$res_list$normalized_counts), file)
    }
  )
   
}

##########################Run SHiny APP##############

shinyApp(ui=ui, server= server)

