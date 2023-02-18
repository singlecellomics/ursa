############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Alioth: scRNASEQ
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last Update Date: 2022-08-31
############################################################################################################
#' @include ini.R
#' @include common.R
#' @import celltalker
#' @import ComplexHeatmap
#' @import cowplot
#' @import data.table
#' @import dplyr
#' @import ggbeeswarm
#' @import ggplot2
#' @import ggpubr
#' @import ggrepel
#' @import ggridges
#' @import ggthemes
#' @import gplots
#' @import gridExtra
#' @import HGNChelper
#' @import patchwork
#' @import plot3D
#' @import plyr
#' @import RColorBrewer
#' @import reshape2
#' @import scales
#' @import tidyverse
#' @import viridis
#' @import akmedoids
#' @import celldex
#' @import celltalker
#' @import clusterProfiler
#' @import ggupset
#' @import DOSE
#' @import enrichplot
#' @import ensembldb
#' @import harmony
#' @import HGNChelper
#' @import htmlwidgets
#' @import monocle3
#' @import org.Hs.eg.db
#' @import plotly
#' @import Seurat
#' @import SeuratWrappers
#' @import SingleR
#' @import findPC
#' @import DoubletFinder
#' @import scubi
#' @import randomcoloR
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format scRNASEQ pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification scRNASEQ
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_scRNASEQ' by default.
#' @param input_dir Directory to all input files. Default directory
#' is the current working directory.
#' @param output_dir Output directory. Default directory is the
#' current working directory.
#' A new folder with the given project name with time stamp as suffix
#' will be created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt
#' format files.
#' @param num_genes_lower_bound Lower bound to select for cells with
#' number of genes more than this value. Default is 200.
#' @param num_genes_upper_bound Upper bound to select for cells with
#' number of genes less than or equals to this value Default is 25000.
#' @param mito_cutoff Threshold percentage to select for cells with
#' percentage of mitochondria genes less than this value.
#' This is remove cells which are debris or artefacts and doublets.
#' Default is 5.
#' @param integration_method Integration method.
#' Accepts 'Seurat' or "Harmony' integration. Default is Seurat.
#' @param cc_regression Cell cycle regression method.
#' Accepts 0 for no regression,
#' 1 for regression with two-phases, 2 for regression with phase
#' difference (See Seurat). Default is 0.
#' @param pc_selection_method Select a method to run an automatic selection
#' of the number of principal components (PCs) or Harmonys (if harmony was
#' chosen for integration). Methods include 'none', 'all', 'piecewise linear
#' model', 'first derivative', 'second derivative', 'preceding residual',
#' 'perpendicular line', and 'k-means clustering'. Default is set to 'none',
#' which will ignore this option and use the stated number of PCs in the
#' num_pcs option. If 'all' is used, all methods will be assessed and the
#' final PC number will be determined via the mean of the output of all
#' methods. For more information, please visit:
#' https://github.com/haotian-zhuang/findPC
#' (Zhuang et al., Bioinformatics, 2022)
#' @param num_pcs Number of PCs or Harmonys (if harmony was chosen for
#'  integration) to be used for the analysis. Default to 30. Please
#'  state a reasonable number which is no more than the total number
#'  of cells submitted to avoid running into errors. This option will
#'  only be considered if option pc_selection_method is set to 'none'.
#' @param to_impute Pass 'YES' to perform imputation on individual
#' samples or 'NO' to skip imputation process. Default is set to 'NO'.
#' @param find_doublet Pass 'YES' to perform doublets removal using
#'  DoubletFinder (Christopher S. M., Cell Systems, 2019). Default
#' is set to 'NO'. If set to 'YES', doublets will be removed using
#' the true doublet rate (TDR).
#' @param run_unbias_vis Plot additional unbias visualizations for top
#' features using SCUBI (Wenpin H. & Zhicheng J., Cell Reports Methods,
#' 2021) for dimensionally reduced plots such as UMAP. Pass 'YES' to
#' plot, default to 'NO'.
#' @export

scRNASEQPip <- function(project_name = "Ursa_scRNASEQ",
                      input_dir = "./",
                      output_dir = "./",
                      pheno_file = "study_meta.csv",
                      num_genes_lower_bound = 200,
                      num_genes_upper_bound = 25000,
                      mito_cutoff = 5,
                      integration_method = "Seurat",
                      cc_regression = 0,
                      pc_selection_method = "none",
                      num_pcs = 30,
                      to_impute = "NO",
                      find_doublet = "NO",
                      run_unbias_vis = "NO"){
  print("Initialising pipeline environment..")
  hs <- org.Hs.eg.db
  data("hgnc.table", package="HGNChelper")
  head(hgnc.table)
  hpca.se <- HumanPrimaryCellAtlasData()
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  pheno_data <- pheno_ini(pheno_file, pipeline = "scRNASEQ", isDir = T)
  pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, sep = "_")
  color_conditions <- color_ini()
  ctime <- time_ini()
  sample_colors <- gen_colors(n = length(unique(pheno_data$SID)))
  names(sample_colors) <- unique(pheno_data$SID)
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir)

  sample_files <- list.files(input_dir, pattern = ".*h5", full.names = T, ignore.case = T, recursive = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]

  # Imputation
  source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")

  #######################################################################################################################################
  print("Preparing..")
  data <- NULL
  data_current <- NULL
  annot_names <- NULL
  results <- NULL

  overall_pcs <- NULL
  if(pc_selection_method == "none"){
    overall_pcs <- as.numeric(as.character(num_pcs))
  }else{
    overall_pcs <- pc_selection_method
  }

  for(j in 1:nrow(pheno_data)){
    print(paste("Running ", pheno_data[j,"FILE"], "..", sep = ""))
    annot_names <- c(annot_names, pheno_data[j,"SID"])
    data_current[[j]] <- NULL
    data_current[[j]] <- sample_files[grep(pheno_data[j,"FILE"], sample_files, ignore.case = T)]
    data_current[[j]] <- Read10X_h5(data_current[[j]])
    if((length(data_current[[j]]) == 1) | (length(data_current[[j]]) > 10)){
      data_current[[j]] <- CreateSeuratObject(counts = data_current[[j]], project = annot_names[j], min.cells = 3, min.features = 200)
    }else{
      data_current[[j]] <- CreateSeuratObject(counts = data_current[[j]]$`Gene Expression`, project = annot_names[j], min.cells = 3, min.features = 200)
    }

    data_current[[j]][["Percent_Mito"]] <- PercentageFeatureSet(data_current[[j]], pattern = "^MT-")
    if(max(data_current[[j]][["Percent_Mito"]]) == 0){
      plot_mito <- FALSE
    }else{
      plot_mito <- TRUE
    }

    names(data_current)[j] <- annot_names[j]
    data_current[[j]] <- add_names(data_current[[j]], sample_name = annot_names[j], current_ident = "orig.ident")
    data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, pheno_data[j,which(colnames(pheno_data) != "SAMPLE_ID")])
    plotx <- data_current[[j]]@meta.data
    colnames(plotx)[grep("orig.ident", colnames(plotx), ignore.case = T)] <- "SAMPLE_ID"

    p1 <- NULL
    p2 <- NULL
    p <- NULL
    p1 <- own_violin(plotx, feature = "nFeature_RNA", plotx_title = "No. Genes Detected / Cell", col = color_conditions$tenx[1], title.size = 18)
    p2 <- own_violin(plotx, feature = "nCount_RNA", plotx_title = "No. Molecules Detected / Cell", col = color_conditions$tenx[2], title.size = 18)
    if(plot_mito == TRUE){
      p3 <- own_violin(plotx, feature = "Percent_Mito", plotx_title = "Mitochondria Percent / Cell", col = color_conditions$tenx[3], title.size = 18)
      p <- p1+p2+p3
    }else{p <- p1+p2}

    somePDFPath <- paste(cdir,"1URSA_PLOT_scRNASEQ_PREFILTERING_RNA_INFO_VIOLIN_",pheno_data[j,"SID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=8,pointsize=12)
    print(p+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    p1 <- NULL
    p2 <- NULL
    p <- NULL
    p2 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                      title_name = annot_names[j], col = color_conditions$tenx[5],
                      xlabel = "No. Molecules Detected/Cell", ylabel = "No. Genes Detected/Cell")
    if(plot_mito == FALSE){
      p <- p2 }else{
        p1 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "Percent_Mito",
                          title_name = annot_names[j], col = color_conditions$tenx[4],
                          xlabel = "No. Molecules Detected/Cell", ylabel = "Mitochondria Percent/Cell")
        p <- p1+p2
      }

    somePDFPath <- paste(cdir,"2URSA_PLOT_scRNASEQ_PREFILTERING_RNA_INFO_SCATTERED_",pheno_data[j,"SID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=5,pointsize=12)
    print(p + plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
    dev.off()

    data_current[[j]] <- NormalizeData(data_current[[j]], verbose = TRUE)
    data_current[[j]]@assays$RNA@data@x[is.na(data_current[[j]]@assays$RNA@data@x)] <- 0
    data_current[[j]] <- FindVariableFeatures(data_current[[j]],selection.method = "vst")
    data_current[[j]] <- ScaleData(data_current[[j]], features = row.names(data_current)[j])

    if(nrow(data_current[[j]]) > 2000){
      orig_gene_names <- NULL
      scale_orig_gene_names <- NULL
      if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
        orig_gene_names <- row.names(data_current)[j]
        scale_orig_gene_names <- row.names(data_current[[j]]@assays$RNA@scale.data)
        row.names(data_current[[j]]@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@counts), ignore.case = T)
        row.names(data_current[[j]]@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@data), ignore.case = T)
        row.names(data_current[[j]]@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@scale.data), ignore.case = T)
      }

      data_current[[j]] <- CellCycleScoring(data_current[[j]], g2m.features=g2m.genes, s.features=s.genes, set.ident = TRUE)
      data_current[[j]]@meta.data$Phase <- factor(data_current[[j]]@meta.data$Phase, levels = c("G1","S","G2M"))
      data_current[[j]] <- RunPCA(data_current[[j]], features = c(s.genes, g2m.genes))
      data_current[[j]]@reductions$pca_selected <- data_current[[j]]@reductions$pca

      if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
        row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
        row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
        row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
      }

      data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))

      # G2/M and S Phase Markers: Tirosh et al, 2015
      p1 <- own_2d_scatter(data_current[[j]], "pca", "Phase", "PCA Based on Variable Features")
      p2 <- own_2d_scatter(data_current[[j]], "pca_selected", "Phase", "PCA Based on G2/M and S Phase Markers")

      somePDFPath = paste(cdir,"3URSA_PLOT_scRNASEQ_CELL_CYCLE_PHASES_PCA_BEFORE_CC_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=7,pointsize=12)
      print(p1 / p2 + plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      selected_pcs <- NULL
      if(class(overall_pcs) == "numeric"){
        selected_pcs <- overall_pcs
      }else if (toupper(overall_pcs) == toupper("all")){
        selected_pcs <- findPC(sdev = data_current[[j]]@reductions$pca@stdev, number = 50, method = 'all',aggregate = 'mean')
      }else if(toupper(overall_pcs) %in% toupper(c('piecewise linear model', 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line', 'k-means clustering'))){
        selected_pcs <- findPC(sdev = data_current[[j]]@reductions$pca@stdev, number = 50, method = tolower(overall_pcs))
      }

      print(paste("Selected number of PCs to use for sample ", annot_names[j], ": ", selected_pcs, sep = ""))

      data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca_selected", dims = 1:ifelse(ncol(data_current[[j]]@reductions$pca_selected@cell.embeddings) < selected_pcs, ncol(data_current[[j]]@reductions$pca_selected@cell.embeddings), selected_pcs))
      data_current[[j]]@reductions$umap_selected <- data_current[[j]]@reductions$umap
      data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:selected_pcs)

      p1 <- own_2d_scatter(data_current[[j]], "umap", "Phase", "UMAP Based on Variable Features")
      p2 <- own_2d_scatter(data_current[[j]], "umap_selected", "Phase", "UMAP Based on G2/M and S Phase Markers")

      somePDFPath = paste(cdir,"4URSA_PLOT_scRNASEQ_CELLCYCLE_UMAP_BEFORE_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=7,pointsize=12)
      print(p1 / p2 + plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

    }else{
      data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
      selected_pcs <- NULL
      if(class(overall_pcs) == "numeric"){
        selected_pcs <- overall_pcs
      }else if (toupper(overall_pcs) == toupper("all")){
        selected_pcs <- findPC(sdev = data_current[[j]]@reductions$pca@stdev, number = 50, method = 'all',aggregate = 'mean')
      }else if(toupper(overall_pcs) %in% toupper(c('piecewise linear model', 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line', 'k-means clustering'))){
        selected_pcs <- findPC(sdev = data_current[[j]]@reductions$pca@stdev, number = 50, method = tolower(overall_pcs))
      }
      print(paste("Selected number of PCs to use for sample ", annot_names[j], ": ", selected_pcs, sep = ""))

      data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:selected_pcs)
    }

    metrics <- colnames(data_current[[j]]@meta.data)[grep("nCount_RNA|nFeature_RNA|S.Score|G2M.Score|Percent_Mito", colnames(data_current[[j]]@meta.data), ignore.case = T)]

    somePDFPath = paste(cdir,"5URSA_PLOT_scRNASEQ_CELL_CYCLE_PHASES_UMAP_BEFORE_CC_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=ceiling(length(metrics))/3*5,pointsize=12)
    print(FeaturePlot(data_current[[j]],
                      reduction = "umap",
                      features = metrics,
                      pt.size = 0.05,
                      order = TRUE,
                      min.cutoff = 'q10',
                      label = FALSE)+
            plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""),
                            theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    ########################################################################################################################
    # Filtering
    data_current[[j]] <- subset(data_current[[j]], subset = nFeature_RNA > num_genes_lower_bound & nFeature_RNA <= num_genes_upper_bound & Percent_Mito < mito_cutoff)

    if(ncol(data_current[[j]]) > 200){

      data_current[[j]] <- NormalizeData(data_current[[j]], verbose = TRUE)

      # Imputation
      if(toupper(to_impute) == "YES"){
        print("Performing imputation using ALRA..")
        cimpute_norm <- NULL
        cimpute_norm <- t(as.matrix(data_current[[j]]@assays$RNA@data))
        k_choice <- choose_k(cimpute_norm)
        cimpute_norm <- alra(cimpute_norm,k=k_choice$k)[[3]]
        cimpute_norm <- t(cimpute_norm)
        colnames(cimpute_norm) <- colnames(data_current[[j]])
        data_current[[j]]@assays$RNA@data <- cimpute_norm
      }

      data_current[[j]]@assays$RNA@data[is.na(data_current[[j]]@assays$RNA@data)] <- 0
      data_current[[j]] <- FindVariableFeatures(data_current[[j]],selection.method = "vst")
      data_current[[j]] <- ScaleData(data_current[[j]], features = rownames(data_current)[j])

      if(cc_regression != 0 & (nrow(data_current[[j]][['RNA']]) > 2000)){

        orig_gene_names <- NULL
        scale_orig_gene_names <- NULL

        if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
          orig_gene_names <- row.names(data_current)[j]
          scale_orig_gene_names <- row.names(data_current[[j]]@assays$RNA@scale.data)
          row.names(data_current[[j]]@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@counts), ignore.case = T)
          row.names(data_current[[j]]@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@data), ignore.case = T)
          row.names(data_current[[j]]@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@scale.data), ignore.case = T)
        }

        data_current[[j]] <- CellCycleScoring(data_current[[j]], g2m.features=g2m.genes, s.features=s.genes, set.ident = TRUE)

        data_current[[j]]@meta.data$Phase <- factor(data_current[[j]]@meta.data$Phase, levels = c("G1","S","G2M"))

        if((cc_regression == 1)){
          cc_method <- c("S.Score", "G2M.Score")
          cc_type <- "Phase Scores"

        }else if (cc_regression == 2){
          data_current[[j]]$CC.Difference <- data_current[[j]]$S.Score - data_current[[j]]$G2M.Score
          cc_method <- "CC.Difference"
          cc_type <- "Difference"
        }

        data_current[[j]] <- ScaleData(data_current[[j]], vars.to.regress = cc_method, features = rownames(data_current)[j])
        data_current[[j]] <- RunPCA(data_current[[j]], features = c(s.genes, g2m.genes))
        data_current[[j]]@reductions$pca_selected <- data_current[[j]]@reductions$pca

        if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
          row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
          row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
          row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
        }

        data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
        data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca_selected", dims = 1:ifelse(ncol(data_current[[j]]@reductions$pca_selected@cell.embeddings) < selected_pcs, ncol(data_current[[j]]@reductions$pca_selected@cell.embeddings), selected_pcs))
        data_current[[j]]@reductions$umap_selected <- data_current[[j]]@reductions$umap
        data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:selected_pcs)

        metrics <- colnames(data_current[[j]]@meta.data)[grep("nCount_RNA|nFeature_RNA|S.Score|G2M.Score|Percent_Mito", colnames(data_current[[j]]@meta.data), ignore.case = T)]

        p1 <- NULL
        p2 <- NULL
        p1 <- DimPlot(data_current[[j]], reduction = "pca",split.by = "Phase", cols = color_conditions$tenx) +
          ggtitle(paste("PCA Based on Variable Features", sep = "")) +
          theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
        p2 <- DimPlot(data_current[[j]], reduction = "pca_selected", split.by = "Phase", cols = color_conditions$tenx) +
          ggtitle(paste("PCA Based on G2/M and S Phase Markers", sep = "")) +
          theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

        somePDFPath = paste(cdir,"6URSA_PLOT_scRNASEQ_CELL_CYCLE_PHASES_PCA_POST_CC_",gsub("\\s+","_",toupper(cc_type)),"_REGRESSION_", annot_names[j], "_", project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=10, height=7,pointsize=12)
        print(p1 / p2 + plot_annotation(title = paste("Post Cell Cycle Regression: ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
        dev.off()

        p1 <- NULL
        p2 <- NULL
        p1 <- DimPlot(data_current[[j]], reduction = "umap", cols = color_conditions$tenx, split.by = "Phase") +
          ggtitle(paste("UMAP Based on Variable Features", sep = "")) + theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
        p2 <- DimPlot(data_current[[j]], reduction = "umap_selected", cols = color_conditions$tenx, split.by = "Phase") +
          ggtitle(paste("UMAP Based on G2/M and S Phase Markers", sep = "")) + theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

        somePDFPath = paste(cdir,"7URSA_PLOT_scRNASEQ_CELLCYCLE_UMAP_POST_CC_",gsub("\\s+","_",toupper(cc_type)),"_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=10, height=7,pointsize=12)
        print(p1 / p2 + plot_annotation(title = paste("Post Cell Cycle Regression: ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
        dev.off()

        somePDFPath = paste(cdir,"8URSA_PLOT_scRNASEQ_CELL_CYCLE_PHASES_UMAP_POST_CC_",gsub("\\s+","_",toupper(cc_type)),"_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=10, height=ceiling(length(metrics))/3*5,pointsize=12)
        print(FeaturePlot(data_current[[j]],
                          reduction = "umap",
                          features = metrics,
                          pt.size = 0.05,
                          order = TRUE,
                          min.cutoff = 'q10',
                          label = FALSE)+
                plot_annotation(title = paste("UMAP Post Cell Cycle Regression - ",annot_names[j], sep = ""),
                                theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))))
        dev.off()

      }else{
        DefaultAssay(data_current[[j]]) <- 'RNA'
        data_current[[j]] <- RunPCA(data_current[[j]], features = VariableFeatures(data_current[[j]]))
        data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:selected_pcs)

        metrics <- colnames(data_current[[j]]@meta.data)[grep("nCount_RNA|nFeature_RNA|S.Score|G2M.Score|Percent_Mito", colnames(data_current[[j]]@meta.data), ignore.case = T)]

        somePDFPath = paste(cdir,"8URSA_PLOT_scRNASEQ_NO_CELL_CYCLE_REGRESSION_UMAP_POST_CC_NO_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=10, height=ceiling(length(metrics))/3*5,pointsize=12)
        print(FeaturePlot(data_current[[j]],
                          reduction = "umap",
                          features = metrics,
                          pt.size = 0.05,
                          order = TRUE,
                          min.cutoff = 'q10',
                          label = FALSE)+
                plot_annotation(title = paste("UMAP Features - ",annot_names[j], sep = ""),
                                theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))))
        dev.off()
}

      colnames(data_current[[j]]@assays$RNA@counts) <- gsub("_|-|\\s+|\\t","",colnames(data_current[[j]]@assays$RNA@counts))
      colnames(data_current[[j]]@assays$RNA@data) <- gsub("_|-|\\s+|\\t","",colnames(data_current[[j]]@assays$RNA@data))
      colnames(data_current[[j]]@assays$RNA@scale.data) <- gsub("_|-|\\s+|\\t","",colnames(data_current[[j]]@assays$RNA@scale.data))
      row.names(data_current[[j]]@meta.data) <- gsub("_|-|\\s+|\\t","",row.names(data_current[[j]]@meta.data))

      plotx <- data_current[[j]]@meta.data
      colnames(plotx)[grep("orig.ident", colnames(plotx), ignore.case = T)] <- "SAMPLE_ID"

      if(is.null(data_current[[j]]@commands$NormalizeData.RNA)){
        temp <- NULL
        temp <- NormalizeData(data_current[[j]])
        data_current[[j]]@commands$NormalizeData.RNA <- temp@commands$NormalizeData.RNA
      }
      temp <- NULL
      temp <- data_current[[j]]
      sweep.res.list <- paramSweep_v3(temp, PCs = 1:selected_pcs, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
      cpK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),"pK"]))
      nExp_poi <- round(0.075*nrow(temp@meta.data))
      temp <- doubletFinder_v3(temp, PCs = 1:selected_pcs, pN = 0.25, pK = cpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
      data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, temp@meta.data[,grep("^pANN_|^DF\\.classifications_", colnames(temp@meta.data), ignore.case = T)])
      cdoublets <- NULL
      cclassname <- NULL
      cclassname <- colnames(data_current[[j]]@meta.data)[grep("^DF\\.classifications_", colnames(data_current[[j]]@meta.data), ignore.case = T)]
      cdoublets <- table(data_current[[j]]@meta.data[,cclassname])
      cdoublets <- cdoublets[grep("Doublet", names(cdoublets), ignore.case = T)]
      print(paste(annot_names[j],": found ",cdoublets, " doublets out of ", ncol(data_current[[j]])," cells..removing..", sep = ""))
      Idents(data_current[[j]]) <- cclassname
      data_current[[j]] <- subset(data_current[[j]], cells = row.names(data_current[[j]]@meta.data[which(data_current[[j]]@meta.data[,cclassname] == "Singlet"),]))

      p1 <- NULL
      p2 <- NULL
      p <- NULL
      p1 <- own_violin(plotx, feature = "nFeature_RNA", plotx_title = "No. Genes Detected / Cell", col = color_conditions$tenx[1], title.size = 15)
      p2 <- own_violin(plotx, feature = "nCount_RNA", plotx_title = "No. Molecules Detected / Cell", col = color_conditions$tenx[2], title.size = 15)
      if(plot_mito == TRUE){
        p3 <- own_violin(plotx, feature = "Percent_Mito", plotx_title = "Mitochondria Percent / Cell", col = color_conditions$tenx[3], title.size = 15)
        p <- p1+p2+p3
      }else{
        p <- p1+p2
      }

      somePDFPath = paste(cdir,"9URSA_PLOT_scRNASEQ_POSTFILTERING_RNA_INFO_VIOLIN_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=14, height=8,pointsize=12)
      print(p+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      p2 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        title_name = annot_names[j], col = color_conditions$tenx[5],
                        xlabel = "No. Molecules Detected/Cell", ylabel = "No. Genes Detected/Cell")
      if(plot_mito == TRUE){
        p1 <- own_feature(plotx, feature1 = "nCount_RNA", feature2 = "Percent_Mito",
                          title_name = annot_names[j], col = color_conditions$tenx[4],
                          xlabel = "No. Molecules Detected/Cell", ylabel = "Mitochondria Percent/Cell")
        p <- p1+p2
      }else{
        p <- p2
      }

      somePDFPath = paste(cdir,"10URSA_PLOT_scRNASEQ_POSTFILTERING_RNA_INFO_SCATTERED_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=5,pointsize=12)
      print(p + plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
      dev.off()

      n <- 10
      p1 <- NULL
      p2 <- NULL
      p1 <- VariableFeaturePlot(data_current[[j]], cols = color_conditions$bright[1:2])
      p2 <- LabelPoints(plot = p1, points = VariableFeatures(data_current[[j]])[1:n], repel = TRUE, xnudge = 0, ynudge = 0)

      somePDFPath = paste(cdir,"11URSA_PLOT_scRNASEQ_POSTNORMALIZED_VARIABLE_FEATURE_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=7,pointsize=12)
      print(p2+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
      dev.off()

      Idents(data_current[[j]]) <- "orig.ident"
      data_current[[j]] <- RunPCA(data_current[[j]])

      somePDFPath = paste(cdir,"12URSA_PLOT_scRNASEQ_TOP_FEATURES_ASSOC_PC12_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=6,pointsize=12)
      print(VizDimLoadings(data_current[[j]], dims = 1:2, reduction = "pca",
                           col = color_conditions$manycolors[1])+
              plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
      dev.off()

      somePDFPath = paste(cdir,"13SCA_PLOT_VISUALIZE_PC12_", annot_names[j], "_", project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=8, height=6,pointsize=12)
      print(DimPlot(data_current[[j]], reduction = "pca",cols = color_conditions$tenx) + theme(legend.position = "none")+
              plot_annotation(title = paste(annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      somePDFPath = paste(cdir,"14SCA_PLOT_HEATMAP_PC1_TO_15_TOP_FEATURES_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=8, height=16,pointsize=12)
      print(DimHeatmap(data_current[[j]], dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE, assays = "RNA")+
              plot_annotation(title = paste(annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      data_current[[j]] <- FindNeighbors(data_current[[j]], dims = 1:selected_pcs)
      data_current[[j]] <- FindClusters(data_current[[j]], resolution = 0.8)
      data_current[[j]] <- RunTSNE(data_current[[j]], reduction = "pca", dims = 1:selected_pcs, check_duplicates = FALSE)
      data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:selected_pcs)

      plotx <- gen10x_plotx(data_current[[j]])
      plotx$CLUSTER <- Idents(data_current[[j]])
      plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(as.numeric(as.character(unique(plotx$CLUSTER)))))
      ccolors <- gen_colors(n = length(levels(plotx$CLUSTER)))
      names(ccolors) <- levels(plotx$CLUSTER)

      p1 <- NULL
      p2 <- NULL
      p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER", plot_title = "UMAP",point_size = 1,
                         col = color_conditions$colorful, annot = TRUE, legend_position = "right", numeric = TRUE)
      p2 <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CLUSTER", plot_title = "tSNE",point_size = 1,
                         col = color_conditions$colorful, annot = TRUE, legend_position = "right", numeric = TRUE)

      somePDFPath = paste(cdir,"15URSA_PLOT_scRNASEQ_UMAP_AUTOCLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=9, height=7,pointsize=12)
      print(p1+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      somePDFPath = paste(cdir,"16URSA_PLOT_scRNASEQ_TSNE_AUTOCLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=9, height=7,pointsize=12)
      print(p2+plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      current_clusters <- sort(as.numeric(as.character(unique(plotx$CLUSTER))))
      Idents(data_current[[j]]) <- "seurat_clusters"
      current_out <- NULL
      current_out <- deanalysis(data_current[[j]], current_clusters, plot_title = annot_names[j],group=NULL,de_analysis = "findallmarkers")

      current_data_markers <- data.frame(SAMPLE = annot_names[j], current_out$current_data_markers[order(current_out$current_data_markers$p_val_adj, decreasing = F),])
      write.table(current_data_markers, paste(cdir, "17URSA_TABLE_scRNASEQ_DEGS_NO_FILTER_",annot_names[j],"_",project_name,".txt", sep = ""), quote = F, row.names = F, sep = "\t")

      current_data_markers <- NULL
      de_type <- current_out$de_type
      de_name <- current_out$de_name
      top1 <- current_out$top1
      topn <- current_out$topn
      current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
      ######################### MEDIAN EXPRESSIONS #####################################
      current <- group_medianexpr(current_out$current_data_markers, data = data_current[[j]], group = "seurat_clusters", cell_type = F)
      plot_median <- current$plot_median
      top_markers <- current$top_markers

      plot_median[is.na(plot_median)] <- 0
      somePDFPath = paste(cdir,"18URSA_PLOT_scRNASEQ_HEATMAP_MEDIAN_EXPR_CLUSTER_TOP_", n, "_GENES_", annot_names[j], "_",project_name, ".pdf", sep = "")
      pdf(file=somePDFPath, width=6, height=9,pointsize=12)
      print(complex_heatmap(plot_median, col = color_conditions$BlueYellowRed))
      dev.off()

      #################################### TOP MARKERS + TSNE / UMAP ###################################################################
      wt <- colnames(current_out$current_data_markers)[grep("log.*FC", colnames(current_out$current_data_markers), ignore.case = T)]
      current_clusters <- sort(as.numeric(as.character(unique(top_markers$cluster))), decreasing = F)
      logFC_list <- NULL
      p1 <- NULL
      P2 <- NULL
      P3 <- NULL
      p1 <- ggplot(data=data.frame(top_markers), aes(x=gene, y=data.frame(top_markers)[,wt], fill=cluster)) +
        geom_bar(stat="identity", position=position_dodge())+
        scale_fill_manual(values = ccolors)+
        theme_classic() +
        coord_flip()+
        ylab(wt)+
        facet_wrap(~cluster, scales = "free_y") +
        theme(axis.text.y=element_text(size = 10),
              strip.text.x = element_text(size = 15, face = "bold"),
              legend.title = element_text(size =20, face = "bold"),
              legend.text = element_text(size = 15))

      plotx <- gen10x_plotx(data_current[[j]], groups = NULL)
      plotx$CLUSTER <- data_current[[j]]$seurat_clusters

      p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CLUSTER",
                         plot_title = annot_names[j],col = ccolors,numeric = T,
                         annot = T, legend_position = "right", point_size = 1)
      p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CLUSTER",
                         plot_title = annot_names[j],col = ccolors,numeric = T,
                         annot = T, legend_position = "right", point_size = 1)

      somePDFPath = paste(cdir,"19URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_TSNE_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=16, height=7,pointsize=12)
      grid.arrange(p1, p2, ncol = 2)
      dev.off()

      somePDFPath = paste(cdir,"20URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_UMAP_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=16, height=7,pointsize=12)
      grid.arrange(p1, p3, ncol = 2)
      dev.off()

      top_1_markers <- current_out$current_data_markers %>% group_by(cluster) %>% top_n(n = 1, eval(parse(text = wt)))

      somePDFPath = paste(cdir,"21URSA_PLOT_scRNASEQ_RIDGE_TOP_FC_SIG_GENES_IN_CLUSTERS_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=14, height=length(unique(unlist(top_1_markers$gene)))/4*6,pointsize=12)
      print(RidgePlot(data_current[[j]], features = unique(unlist(top_1_markers$gene)), ncol = 4,
                      cols = ccolors)+
              plot_annotation(title = paste("TOP1: ", annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      somePDFPath = paste(cdir,"22URSA_PLOT_scRNASEQ_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=ceiling(length(top1$gene)/4)*5,pointsize=12)
      print(current_out$featureplot +plot_layout(ncol = 4) +
              plot_annotation(title = paste("TOP1: ",de_name,annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

if(run_unbias_vis == "YES"){
      cumap <- NULL
      cumap <- data_current[[j]]@reductions$umap@cell.embeddings
      cumapplot <- NULL
      cumapplot <- scubi_expression(dim1 = cumap[,"UMAP_1"], dim2 = cumap[,"UMAP_2"], count = data_current[[j]]@assays$RNA@counts, gene = as.list(current_out$top1$gene))

      somePDFPath = paste(cdir,"22URSA_PLOT_scRNASEQ_SCUBI_UNBIASED_VIS_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=20, height=ceiling(length(top1$gene)/4)*5,pointsize=12)
      print(cumapplot + theme_classic(base_size = 20) + plot_annotation(title = paste("TOP1: ",de_name,annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()
}

      somePDFPath = paste(cdir,"23URSA_PLOT_scRNASEQ_VIOLIN_TOP_",n,"_MARKERS_IN_CLUSTERS_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=16, height=10,pointsize=12)

      for(k in 1:length(current_clusters)){
        current <- topn[which(topn$cluster == current_clusters[k]),]
        current <- current[order(current$p_val_adj, decreasing = F),]
        print(VlnPlot(data_current[[j]], features = current$gene, pt.size = 0,
                      ncol = 4, cols = gen_colors(n = length(unique(data_current[[j]]$seurat_clusters))))&
                xlab("CLUSTERS")&
                plot_annotation(title = paste("TOP",n,": ",current_out$de_name, annot_names[j], " - CLUSTER ", current_clusters[k], sep = ""),
                                theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))
      }
      dev.off()

      somePDFPath = paste(cdir,"24URSA_PLOT_scRNASEQ_HEATMAP_TOP_",n,"_MARKERS_IN_CLUSTERS_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=length(unique(topn$gene))*0.1,pointsize=12)
      print(DoHeatmap(data_current[[j]], features = topn$gene,
                      group.colors = ccolors, size = 8) +
              ggtitle(annot_names[j])+
              NoLegend()+theme(axis.text.x = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.title.x = element_blank(),
                               axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                               legend.title = element_blank(),
                               legend.text = element_text(size = 15),
                               plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
      dev.off()

      Idents(data_current[[j]]) <- "seurat_clusters"
      DefaultAssay(data_current[[j]]) <- "RNA"

      orig_gene_names <- NULL
      scale_orig_gene_names <- NULL
      if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
        orig_gene_names <- row.names(data_current)[j]
        scale_orig_gene_names <- row.names(data_current[[j]]@assays$RNA@scale.data)
        row.names(data_current[[j]]@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@counts), ignore.case = T)
        row.names(data_current[[j]]@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@data), ignore.case = T)
        row.names(data_current[[j]]@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data_current[[j]]@assays$RNA@scale.data), ignore.case = T)
      }

      clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data_current[[j]])),
                         clusters =data_current[[j]]$seurat_clusters,
                         ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.fine)
      data_current[[j]]$CELL_TYPE <- clu_ann$labels[match(data_current[[j]]$seurat_clusters,row.names(clu_ann))]
      clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data_current[[j]])),
                         clusters =data_current[[j]]$seurat_clusters,
                         ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.main)
      data_current[[j]]$CELL_TYPE_LEVEL2 <- clu_ann$labels[match(data_current[[j]]$seurat_clusters,row.names(clu_ann))]
      data_current[[j]]@meta.data[which(is.na(data_current[[j]]$CELL_TYPE)),"CELL_TYPE"] <- "Unidentifiable"
      if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
        row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
        row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
        row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
      }

      Idents(data_current[[j]]) <- data_current[[j]]$CELL_TYPE
      plotx <- gen10x_plotx(data_current[[j]], include_meta = T)

      write.table(plotx, paste(cdir, "25URSA_TABLE_scRNASEQ_UMAP_TSNE_PCA_COORDINATES_CLUSTERS_",annot_names[j],"_",project_name,".txt", sep = ""), quote = F, row.names = F, sep = "\t")

      cct_colors <- gen_colors(n = length(unique(plotx$CELL_TYPE)))
      names(cct_colors) <- sort(unique(plotx$CELL_TYPE))

      somePDFPath = paste(cdir,"26URSA_PLOT_scRNASEQ_UMAP_AUTO_CELL_TYPE_IDENTIFICATION_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=9, height=6,pointsize=12)
      print(plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = annot_names[j],
                         col = cct_colors, annot = TRUE, label_size = 3, legend_position = "right", point_size = 0.1, legendsize = 10))
      dev.off()

      somePDFPath = paste(cdir,"27URSA_PLOT_scRNASEQ_TSNE_AUTO_CELL_TYPE_IDENTIFICATION_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=9, height=6,pointsize=12)
      print(plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CELL_TYPE", plot_title = annot_names[j],
                         col = cct_colors, annot = TRUE, label_size = 3, legend_position = "right", point_size = 0.1, legendsize = 10))
      dev.off()

      ################ Pathway Analysis #########################################################################
      filtered_markers <- current_out$current_data_markers[grep("AC[0-9]+\\.[0-9]+|AL[0-9]+\\.[0-9]+",current_out$current_data_markers$gene, ignore.case = T, invert = T),]
      topn <- split(filtered_markers, filtered_markers$cluster)
      topn <- lapply(topn, function(x){
      x <- x[which(abs(x$avg_log2FC) > 0.25),]
      })
      topn <- do.call(rbind.data.frame, topn)
      topn <- topn[order(topn$avg_log2FC, decreasing = T),]
      if(length(grep("ENS.*-.*", topn$gene)) > length(topn$gene)/2){
        topn$ENSEMBL_ID <- gsub("(ENS.*?[0-9]+)[-|_|\\s+|\\.].*","\\1",topn$gene)
        mapped_id <- ensembldb::select(hs, keys = topn$ENSEMBL_ID, columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
        topn$ENTREZ_ID <- mapped_id[match(topn$ENSEMBL_ID, mapped_id$ENSEMBL),"ENTREZID"]
      }else{
        mapped_id <- ensembldb::select(hs, keys = topn$gene, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
        topn$ENTREZ_ID <- mapped_id[match(topn$gene, mapped_id$SYMBOL),"ENTREZID"]
      }

      if(length(which(is.na(topn$ENTREZ_ID))) > 0 & !(length(grep("ENS.*-.*", topn$gene)) > length(topn$gene)/2)){
        current <- topn[which(is.na(topn$ENTREZ_ID)),]
        current$alternate_gene <- hgnc.table[match(current$gene, hgnc.table$Symbol),"Approved.Symbol"]
        current <- current[which(current$gene != current$alternate_gene),]
        if(nrow(current) > 0){
          current_id <- ensembldb::select(hs, keys = current$alternate_gene[!is.na(current$alternate_gene)], columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
          current$ENTREZ_ID <- current_id[match(current$alternate_gene, current_id$SYMBOL),"ENTREZID"]
          topn[which(is.na(topn$ENTREZ_ID)),"ENTREZ_ID"] <- current[match(as.character(unlist(topn[which(is.na(topn$ENTREZ_ID)),"gene"])), current$gene),"ENTREZ_ID"]
        }
      }
      topn <- topn[!is.na(topn$ENTREZ_ID),]

      p_threshold <- 0.05
      current <- split(topn, topn$cluster)
      pathway_EA_result <- NULL
      pathway_EA_result <- lapply(current, function(x){
        x <- enrichDGN(gene=unique(as.numeric(as.character(x$ENTREZ_ID))),pvalueCutoff=0.05, qvalueCutoff = 0.05,readable=T)
      })

      for(i in 1:length(pathway_EA_result)){
        if(!is.null(pathway_EA_result[[i]])){
          pathway_EA_result[[i]]@result <- pathway_EA_result[[i]]@result[which(pathway_EA_result[[i]]@result$p.adjust < p_threshold),]
        }
      }

      results[['pathway_EA_result']][[j]] <- pathway_EA_result
      names(results[['pathway_EA_result']])[j] <- pheno_data[j,"SID"]

      summary_pathways <- NULL
      n <- 15
      somePDFPath = paste(cdir,"28URSA_PLOT_scRNASEQ_TOP_",n,"_PATHWAYS_IN_EACH_CLUSTER_TOPFC200GENES_CLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=12, height=6,pointsize=12)
      for(i in 1:length(pathway_EA_result)){
        current <- pathway_EA_result[[i]]
        if(!is.null(current)){
          if(nrow(current@result) > 0){
            summary_pathways <- rbind(summary_pathways, data.frame(CLUSTER = names(pathway_EA_result)[i],current))
            print(barplot(current, showCategory=n, orderBy = "y")+
              ggtitle(paste(annot_names[j],"\nEnriched Terms for ORA: Cluster ", names(pathway_EA_result)[i],sep = "")))
          }
        }
      }

      dev.off()

      write.table(summary_pathways, paste(cdir, "29URSA_TABLE_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_ENRICHMENT_TERMS_",annot_names[j], "_",project_name,".txt",sep = ""), quote = F, row.names = F, sep = "\t")

      somePDFPath = paste(cdir,"30URSA_PLOT_scRNASEQ_DOT_PLOT_TOP_ABS0.25_AVG_LOG2FC_PATHWAYS_IN_EACH_CLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=12, height=6,pointsize=12)
      for(i in 1:length(pathway_EA_result)){
        current <- pathway_EA_result[[i]]
        if(!is.null(current)){
          if(nrow(current@result) > 0){
            print(dotplot(current, showCategory=n, orderBy = "x")+ggtitle(paste(annot_names[j], "\nEnriched Terms for ORA: Cluster ", names(pathway_EA_result)[i],sep = "")))
          }
        }
      }
      dev.off()

      topn_cluster <- split(topn, topn$cluster)
      somePDFPath = paste(cdir,"31URSA_PLOT_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_GENES_CONCEPT_NETWORK_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=35, height=9,pointsize=12)

      for(i in 1:length(pathway_EA_result)){
        current <- pathway_EA_result[[i]]
        if(!is.null(current)){
          if(nrow(current@result) > 0){
            current_gene <- topn_cluster[[i]]
            if(toupper(wt) == toupper("avg_log2FC")){
              gene_list <- current_gene$avg_log2FC
            }else{
              gene_list <- current_gene$avg_logFC
            }
            names(gene_list) <- current_gene$ENTREZ_ID
            gene_list <- sort(gene_list, decreasing = TRUE)
            gene_list <- gene_list[!duplicated(names(gene_list))]

            p1 <- NULL
            p2 <- NULL
            p3 <- NULL
            p1 <- cnetplot(current, foldChange=gene_list) +
              ggtitle(paste(annot_names[j], "\nA: Gene-Concept Network, ",wt,": Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
            p2 <- cnetplot(current, categorySize="pvalue", foldChange=gene_list) +
              ggtitle(paste(annot_names[j], "\nB: Gene-Concept Network, ",wt," & P-Value: Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
            p3 <- cnetplot(current, foldChange=gene_list, circular = TRUE, colorEdge = TRUE) +
              ggtitle(paste(annot_names[j], "\nC: Gene-Concept Network, ",wt,": Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))

            print(cowplot::plot_grid(p1, p2, p3, ncol=3, rel_widths=c(1.2, 1.2, 1.2)))
          }
        }
      }
      dev.off()

      somePDFPath = paste(cdir,"32URSA_PLOT_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_NETWORK_ENRICHMENT_TERMS_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=22, height=8,pointsize=12)
      for(i in 1:length(pathway_EA_result)){
        current <- pathway_EA_result[[i]]
        if(!is.null(current)){
          if(nrow(current@result) > 0){
            p1 <- NULL
            p2 <- NULL
            p1 <- cnetplot(current, node_label="category") +
              ggtitle(paste(annot_names[j], "\nA: Network Terms: Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))
            p2 <- cnetplot(current, node_label="gene") +
              ggtitle(paste(annot_names[j], "\nA: Network Genes: Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm"))

            print(cowplot::plot_grid(p1, p2, ncol=2, rel_widths=c(1.2, 1.2, 1.2)))
          }
        }
      }
      dev.off()

      somePDFPath = paste(cdir,"33URSA_PLOT_scRNASEQ_HEATMAP_TOP_ABS0.25_AVG_LOG2FC_ENRICHED_NETWORK_TERMS_BY_CLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=22, height=6,pointsize=12)
      for(i in 1:length(pathway_EA_result)){
        current <- pathway_EA_result[[i]]
        if(!is.null(current)){
          if(nrow(current@result) > 0){
            current_gene <- topn_cluster[[i]]
            if(toupper(wt) == toupper("avg_log2FC")){
              gene_list <- current_gene$avg_log2FC
            }else{
              gene_list <- current_gene$avg_logFC
            }
            names(gene_list) <- current_gene$ENTREZ_ID
            gene_list <- sort(gene_list, decreasing = TRUE)
            gene_list <- gene_list[!duplicated(names(gene_list))]
            # currentR <- setReadable(current, 'org.Hs.eg.db', 'ENTREZID')
            print(heatplot(current, foldChange = gene_list, label_format = 20) +
              ggtitle(paste(annot_names[j], "\nHeatmap of Enrichment Terms: Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold")))
          }
        }
      }
      dev.off()

      somePDFPath = paste(cdir,"34URSA_PLOT_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_ENRICHMENT_PROPORTIONS_GENES_CLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=12, height=10,pointsize=12)
      for(i in 1:length(pathway_EA_result)){
        current <- pathway_EA_result[[i]]
        if(!is.null(current)){
          if(nrow(current@result) > 0){
            currentR <- pairwise_termsim(current)
            if(nrow(currentR@result) > 10){
              print(emapplot(currentR, cex_line = 0.1) +
                ggtitle(paste(annot_names[j], "\nEnrichment Map: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                scale_colour_gradient_tableau(palette = "Red-Gold") +
                theme(plot.title = element_text(face = "bold")))
            }
          }
        }
      }
      dev.off()

      current <- NULL
      for(i in 1:length(topn_cluster)){
        if(!is.null(topn_cluster[[i]])){
          current[[i]] <- topn_cluster[[i]]$ENTREZ_ID
        }
      }

      names(current) <- names(topn_cluster)
      topn_cluster <- topn_cluster[!unlist(lapply(topn_cluster,is.null))]
      current <- current[!unlist(lapply(current,is.null))]
      plotx <- compareCluster(current, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)
      currentR <- pairwise_termsim(plotx)

      somePDFPath = paste(cdir,"35URSA_PLOT_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_ENRICHMENT_PROPORTIONS_GENES_CLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=9, height=8,pointsize=12)
      print(emapplot(currentR, cex_line = 0.1, pie="count", pie_scale=1.5, layout="kk") +
        ggtitle(paste(annot_names[j], "\nEnrichment Map with Proportions", sep = "")) +
        scale_fill_tableau(palette = "Tableau 20") +
        theme(plot.title = element_text(face = "bold")))
      dev.off()

      somePDFPath = paste(cdir,"36URSA_PLOT_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_KEGG_ENRICHMENT_MAP_PROPORTION_CLUSTER_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=6,pointsize=12)
      for(i in 1:length(pathway_EA_result)){
        if(!is.null(pathway_EA_result[[i]])){
          if(nrow(pathway_EA_result[[i]]@result) > 0){
            print(upsetplot(pathway_EA_result[[i]]) +
              ggtitle(paste(annot_names[j], "\nUpset Plot: Cluster ", names(pathway_EA_result)[i],sep = "")) +
              theme(plot.title = element_text(face = "bold"),plot.margin = unit(c(1,1,1,1), "cm")))
          }
        }
      }
      dev.off()

      GSEA_result <- NULL
      for(i in 1:length(topn_cluster)){
        if(toupper(wt) == toupper("avg_logFC")){
          gene_list <- topn_cluster[[i]]$avg_logFC
        }else{
          gene_list <- topn_cluster[[i]]$avg_log2FC
        }
        if(nrow(topn_cluster[[i]])>10){
          names(gene_list) <- topn_cluster[[i]]$ENTREZ_ID
          gene_list <- sort(gene_list, decreasing = TRUE)
          gene_list <- gene_list[!duplicated(names(gene_list))]
          GSEA_result[[i]] <- gseNCG(gene=gene_list, pvalueCutoff = p_threshold)
        }else{
          GSEA_result[[i]] <- NULL
        }
      }

      if_plot <- NULL

      for(i in 1:length(GSEA_result)){
        if(!is.null(GSEA_result[[i]])){
          if(nrow(GSEA_result[[i]]@result) > 0){
            GSEA_result[[i]]@result <- GSEA_result[[i]]@result[GSEA_result[[i]]@result$p.adjust < p_threshold,]
            if_plot[[i]] <- ifelse(nrow(GSEA_result[[i]]@result) ==0, FALSE, TRUE)
          }else{
            if_plot[[i]] <- FALSE
          }
        }else{
          if_plot[[i]] <- FALSE
        }
      }

      if(!all(if_plot == FALSE)){

        somePDFPath = paste(cdir,"37URSA_PLOT_scRNASEQ_TOP_ABS0.25_AVG_LOG2FC_GSEA_GENE_ENRICHMENT_SCORE_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=10, height=7,pointsize=12)
        for(i in 1:length(GSEA_result)){
          current <- GSEA_result[[i]]
          if(!is.null(current)){
            if(nrow(current@result[which(current@result$ID != "-"),]) > 0){
              print(gseaplot2(current, geneSetID = which(current@result$ID != "-"), title = current$Description[which(current$Description != "-")], pvalue_table = TRUE,
ES_geom = "line")+
ggtitle(paste(annot_names[j], "\nGSEA Plot: Cluster ", names(pathway_EA_result)[i],sep = "")) +
                theme(plot.title = element_text(face = "bold")))
            }
          }
        }

        dev.off()
      }

      ##################### CELL TYPES ##########################################################################

      if(length(unique(data_current[[j]]$CELL_TYPE)) > 1){
        current <- group_medianexpr(current_out$current_data_markers, data = data_current[[j]], group = "CELL_TYPE", cell_type = T)
        plot_median_cell_type <- current$plot_median
        top_markers_cell_type <- current$top_markers
        n <- length(unique(top_markers_cell_type$gene))

        plot_median_cell_type[is.na(plot_median_cell_type)] <- 0
        somePDFPath = paste(cdir,"38SCA_PLOT_HEATMAP_MEDIAN_EXPR_CELL_TYPE_TOP_FC_", annot_names[j], "_", project_name, ".pdf", sep = "")
        pdf(file=somePDFPath, width=6, height=10,pointsize=12)
        print(complex_heatmap(plot_median_cell_type, col = color_conditions$BlueYellowRed))
        dev.off()

        cell_types <- sort(as.character(unique(top_markers_cell_type$CELL_TYPE)))
        colors <- gen_colors(n = length(cell_types))

        p1 <- NULL
        p2 <- NULL
        p3 <- NULL
        p1 <- ggplot(data=data.frame(top_markers_cell_type), aes(x=gene, y=data.frame(top_markers_cell_type)[,wt], fill=CELL_TYPE)) +
          geom_bar(stat="identity", position=position_dodge())+
          scale_fill_manual(values = as.character(colors))+
          theme_classic() +
          coord_flip()+
          ylab(wt)+
          facet_wrap(~CELL_TYPE, scales = "free_y") +
          theme(axis.text.y=element_text(size = 10),
                strip.text.x = element_text(size = 15, face = "bold"),
                legend.title = element_text(size =20, face = "bold"),
                legend.text = element_text(size = 15))

        plotx <- gen10x_plotx(data_current[[j]], groups = NULL)
        plotx$CELL_TYPE <- data_current[[j]]$CELL_TYPE

        p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CELL_TYPE",
                           plot_title = annot_names[j],col = cct_colors,numeric = F,
                           annot = T, legend_position = "right", point_size = 1, legendsize = 10)
        p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CELL_TYPE",
                           plot_title = annot_names[j],col = cct_colors,numeric = F,
                           annot = T, legend_position = "right", point_size = 1, legendsize = 10)

        somePDFPath = paste(cdir,"39URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_TSNE_CELL_TYPES_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=24, height=ceiling(length(unique(plotx$CELL_TYPE))/3)*5,pointsize=12)
        grid.arrange(p1, p2, ncol = 2)
        dev.off()

        somePDFPath = paste(cdir,"40URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_UMAP_CELL_TYPES_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=24, height=ceiling(length(unique(plotx$CELL_TYPE))/3)*5,pointsize=12)
        grid.arrange(p1, p3, ncol = 2)
        dev.off()

        current_out$current_data_markers$CELL_TYPE <- data_current[[j]]@meta.data[match(current_out$current_data_markers$cluster, data_current[[j]]$seurat_clusters),"CELL_TYPE"]
        top_1_markers <- current_out$current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = 1, eval(parse(text = wt)))
        Idents(data_current[[j]]) <- "CELL_TYPE"
        somePDFPath = paste(cdir,"41URSA_PLOT_scRNASEQ_RIDGE_TOP_FC_SIG_GENES_CELL_TYPES_",annot_names[j],"_",project_name,".pdf", sep = "")
        pdf(file=somePDFPath, width=16, height=length(unique(unlist(top_1_markers$gene)))/2*6,pointsize=12)
        print(RidgePlot(data_current[[j]], features = unique(unlist(top_1_markers$gene)), ncol = 2,
                        cols = cct_colors)+
                plot_annotation(title = annot_names[j], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
        dev.off()

        Idents(data_current[[j]]) <- "orig.ident"
        DefaultAssay(data_current[[j]]) <- 'RNA'
        rm(p)
      }
    }
  }
  ########## Integration Steps #####################################################################
  if(nrow(pheno_data) == 1){

    data <- data_current[[1]]
    DefaultAssay(data) <- "RNA"
    n <- 10
    top_features <- data.frame(TOP_PCA_POS_NEG_GENES = paste(capture.output(print(data[["pca"]], dims = 1:5, nfeatures = n))))

  }else{
    data_current <- lapply(data_current, function(x){
      # x <- ScaleData(x, verbose=F, features = data_features, vars.to.regress = c("nCount_RNA", "Percent_Mito"))
      # x <- RunPCA(x, npcs = 30, verbose = FALSE, features = data_features)
      # x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      # x <- RunPCA(x, npcs = 30, verbose = FALSE, features = VariableFeatures(x))
    })
    integrative_features <- SelectIntegrationFeatures(object.list = data_current)
    data_anchors <- FindIntegrationAnchors(object.list = data_current,
                                           reduction = "rpca", anchor.features = integrative_features)
    data <- IntegrateData(anchorset = data_anchors, k.weight = 50)
    rm(data_anchors)
    DefaultAssay(data) <- "integrated"
    data <- ScaleData(data, verbose = FALSE)

    DefaultAssay(data) <- "RNA"
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    data <- ScaleData(data, verbose = FALSE)
    data <- RunPCA(data, verbose = FALSE)

    selected_pcs <- NULL
    if(class(overall_pcs) == "numeric"){
      selected_pcs <- overall_pcs
    }else if (toupper(overall_pcs) == toupper("all")){
      selected_pcs <- findPC(sdev = data@reductions$pca@stdev, number = 50, method = 'all',aggregate = 'mean')
    }else if(toupper(overall_pcs) %in% toupper(c('piecewise linear model', 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line', 'k-means clustering'))){
      selected_pcs <- findPC(sdev = data@reductions$pca@stdev, number = 50, method = tolower(overall_pcs))
    }

    data <- RunUMAP(data, reduction = "pca", dims = 1:selected_pcs)
    data <- RunTSNE(data, reduction = "pca", dims = 1:selected_pcs, check_duplicates = FALSE)
    data_dim <- data.frame(gen10x_plotx(data), DATA_TYPE = "BEFORE_INTEGRATION", SAMPLE_ID = data$orig.ident)

    write.table(data_dim, paste(cdir, "42URSA_TABLE_scRNASEQ_DIM_PARAMETERS_BEFORE_INTEGRATION_", project_name,".txt", sep = ""), quote = F, row.names = T, sep = "\t")

    if((toupper(integration_method) == "SEURAT") | (toupper(integration_method) == "NULL")){
      # integration_method <- "SEURAT"
      integration_name <- "SEURAT_INTEGRATED"
      DefaultAssay(data) <- "integrated"
      reduction_method <- "pca"

      data <- RunPCA(data, verbose = FALSE)
      data <- RunUMAP(data, reduction = "pca", dims = 1:selected_pcs)
      data <- RunTSNE(data, reduction = "pca", dims = 1:selected_pcs, check_duplicates = FALSE)
      current <- cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident))
    }else if(toupper(integration_method) == "HARMONY"){
      integration_name <- "HARMONY_INTEGRATED"
      DefaultAssay(data) <- "RNA"
      reduction_method <- "harmony"
      data <- RunHarmony(data, group.by.vars = "BATCH")

      selected_pcs <- NULL
      if(class(overall_pcs) == "numeric"){
        selected_pcs <- overall_pcs
      }else if (toupper(overall_pcs) == toupper("all")){
        selected_pcs <- findPC(sdev = sort(data@reductions$harmony@stdev, decreasing = T), number = 50, method = 'all',aggregate = 'mean')
      }else if(toupper(overall_pcs) %in% toupper(c('piecewise linear model', 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line', 'k-means clustering'))){
        selected_pcs <- findPC(sdev = sort(data@reductions$harmony@stdev, decreasing = T), number = 50, method = tolower(overall_pcs))
      }

      data <- RunUMAP(data, reduction = "harmony", dims = 1:selected_pcs)
      data <- RunTSNE(data, reduction = "harmony", dims = 1:selected_pcs, check_duplicates = FALSE)
      current <- cbind(data.frame(gen10x_plotx(data, selected = c("HARMONY","UMAP","TSNE"), include_meta = T), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident))
    }

    write.table(current, paste(cdir, "43URSA_TABLE_scRNASEQ_DIM_PARAMETERS_AFTER_", integration_method, "_INTEGRATION_", project_name,".txt", sep = ""), quote = F, row.names = T, sep = "\t")

    data_dim <- rbind(data_dim, cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident)))

    somePDFPath = paste(cdir,"44URSA_PLOT_scRNASEQ_UMAP_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=10,pointsize=12)
    print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                              data_dim[data_dim$DATA_TYPE == integration_name,],
                              dim1= "UMAP_1", dim2 = "UMAP_2", group = "SAMPLE_ID",
                              subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                              maintitle = project_name, titlesize = 35, col = sample_colors))
    dev.off()

    somePDFPath = paste(cdir,"45URSA_PLOT_scRNASEQ_TSNE_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=10,pointsize=12)
    print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                              data_dim[data_dim$DATA_TYPE == integration_name,],
                              dim1= "tSNE_1", dim2 = "tSNE_2", group = "SAMPLE_ID",
                              subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                              maintitle = project_name, titlesize = 35, col = sample_colors))
    dev.off()

    if(toupper(integration_method) == "SEURAT"){

      somePDFPath = paste(cdir,"46URSA_PLOT_scRNASEQ_PCA_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=30, height=15,pointsize=12)
      print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                   data_dim[data_dim$DATA_TYPE == integration_name,],
                                                   dim1= "PC_1", dim2 = "PC_2", group = "SAMPLE_ID",
                                                   subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                   maintitle = project_name, titlesize = 35, col = sample_colors))
      dev.off()
    }else if(toupper(integration_method) == "HARMONY"){
      p1 <- NULL
      p2 <- NULL
      p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",], x = "PC_1", y = "PC_2", group = "SAMPLE_ID", plot_title = "BEFORE_INTEGRATION",
                         col = sample_colors, annot = FALSE, legend_position = "bottom")
      p2 <- plot_bygroup(current, x = "HARMONY_1", y = "HARMONY_2", group = "SAMPLE_ID", plot_title = integration_name,
                         col = sample_colors, annot = FALSE, legend_position = "bottom")

      somePDFPath = paste(cdir,"46URSA_PLOT_scRNASEQ_PCA_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=30, height=15,pointsize=12)
      print(p1+p2+
              plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5))))
      dev.off()
    }

    data_dim$GROUP <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"GROUP"]
    data_dim$BATCH <- data@meta.data[match(data_dim$SAMPLE_ID, data$orig.ident),"BATCH"]

    cgroup_colors <- gen_colors(n = length(unique(data_dim$GROUP)))
    names(cgroup_colors) <- sort(unique(data_dim$GROUP))

    somePDFPath = paste(cdir,"47URSA_PLOT_scRNASEQ_BY_GROUP_UMAP_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=30, height=15,pointsize=12)
    print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                              data_dim[data_dim$DATA_TYPE == integration_name,],
                              dim1= "UMAP_1", dim2 = "UMAP_2",group = "GROUP",
                              subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                              maintitle = project_name, titlesize = 35, col = cgroup_colors))
    dev.off()

    somePDFPath = paste(cdir,"48URSA_PLOT_scRNASEQ_BY_GROUP_TSNE_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=30, height=15,pointsize=12)
    print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                 data_dim[data_dim$DATA_TYPE == integration_name,],
                                                 dim1= "tSNE_1", dim2 = "tSNE_2",group = "GROUP",
                                                 subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                 maintitle = project_name, titlesize = 35, col = cgroup_colors))
    dev.off()

    p1 <- NULL
    p2 <- NULL
    p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",], x = "PC_1", y = "PC_2", group = "GROUP", plot_title = "BEFORE_INTEGRATION",
                       col = cgroup_colors, annot = FALSE, legend_position = "bottom", point_size = 1)
    if(toupper(integration_method) == "SEURAT"){
      p2 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == integration_name,], x = "PC_1", y = "PC_2", group = "GROUP", plot_title = integration_name,
                         col = cgroup_colors, annot = FALSE, legend_position = "bottom", point_size = 1)
    }else{
      p2 <- plot_bygroup(current, x = "HARMONY_1", y = "HARMONY_2", group = "GROUP", plot_title = integration_name,
                         col = cgroup_colors, annot = FALSE, legend_position = "bottom", point_size = 1)

    }

    somePDFPath = paste(cdir,"49URSA_PLOT_scRNASEQ_BY_GROUP_PCA_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=30, height=15,pointsize=12)
    print(p1+p2+
            plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5))))
    dev.off()

    cbatch_colors <- gen_colors(n = length(unique(data_dim$BATCH)))
    names(cbatch_colors) <- sort(unique(data_dim$BATCH))

    somePDFPath = paste(cdir,"50URSA_PLOT_scRNASEQ_BY_BATCH_UMAP_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=30, height=15,pointsize=12)
    print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                              data_dim[data_dim$DATA_TYPE == integration_name,],
                              dim1= "UMAP_1", dim2 = "UMAP_2", group = "BATCH",
                              subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                              maintitle = project_name, titlesize = 35, col = cbatch_colors))
    dev.off()

    somePDFPath = paste(cdir,"51URSA_PLOT_scRNASEQ_BY_BATCH_TSNE_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=30, height=15,pointsize=12)
    print(beforeafter_dimplot(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",],
                                                 data_dim[data_dim$DATA_TYPE == integration_name,],
                                                 dim1= "tSNE_1", dim2 = "tSNE_2", group = "BATCH",
                                                 subtitle1 = "BEFORE_INTEGRATION", subtitle2 = integration_name,
                                                 maintitle = project_name, titlesize = 35, col = cbatch_colors))
    dev.off()

    p1 <- NULL
    p2 <- NULL
    p1 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == "BEFORE_INTEGRATION",], x = "PC_1", y = "PC_2", group = "BATCH", plot_title = "BEFORE_INTEGRATION",point_size = 1,
                       col = cbatch_colors, annot = FALSE, legend_position = "bottom")
    if(toupper(integration_method) == "SEURAT"){
      p2 <- plot_bygroup(data_dim[data_dim$DATA_TYPE == integration_name,], x = "PC_1", y = "PC_2", group = "BATCH", plot_title = integration_name,point_size = 1,
                         col = cbatch_colors, annot = FALSE, legend_position = "bottom")
    }else{
      p2 <- plot_bygroup(current, x = "HARMONY_1", y = "HARMONY_2", group = "BATCH", plot_title = integration_name,
                         col = cbatch_colors, point_size = 1,annot = FALSE, legend_position = "bottom")
    }

    somePDFPath = paste(cdir,"52URSA_PLOT_scRNASEQ_BY_BATCH_PCA_BEFORE_AFTER_",integration_method, "_INTEGRATION_", project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=30, height=15,pointsize=12)
    print(p1+p2+
            plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5))))
    dev.off()

    DefaultAssay(data) <- "integrated"

    n <- 10
    top_features <- data.frame(TOP_PCA_POS_NEG_GENES = paste(capture.output(print(data[[ifelse(toupper(integration_method) == "SEURAT", "pca", "harmony")]], dims = 1:5, nfeatures = n))))

    somePDFPath = paste(cdir,"53URSA_PLOT_scRNASEQ_INTEGRATED_DATA_TOP_FEATURES_ASSOC_PC12_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=7,pointsize=12)
    print(VizDimLoadings(data, dims = 1:2, reduction = ifelse(toupper(integration_method) == "SEURAT", "pca", "harmony"),col = color_conditions$mark[1])+ plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))))
    dev.off()

    ############## Integration Plots #################################################################
    DefaultAssay(data) <- "integrated"
    somePDFPath = paste(cdir,"54URSA_PLOT_scRNASEQ_HEATMAP_PC1_TO_15_TOP_FEATURES_INTEGRATED_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=8, height=16,pointsize=12)
    print(DimHeatmap(data, assays = ifelse(toupper(integration_method) == "SEURAT", "integrated", "RNA"),
                                        reduction = ifelse(toupper(integration_method) == "SEURAT", "pca", "harmony"),
                                        dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)+
      plot_annotation(title = paste(project_name, sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    if(toupper(integration_method) == "HARMONY"){
      DefaultAssay(data) <- "RNA"
    }

    data <- FindNeighbors(data, reduction = reduction_method, dims = 1:selected_pcs)
    data <- FindClusters(data, resolution = 0.8)

    if(toupper(integration_method) == "HARMONY"){
      integration_cluster <- "RNA_snn_res.0.8"
    }else{
      integration_cluster <- "integrated_snn_res.0.8"
    }

    plotx <- gen10x_plotx(data)
    plotx$CLUSTER <- Idents(data)
    plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(as.numeric(as.character(unique(plotx$CLUSTER)))))
    plotx$GROUP <- data$GROUP
    current_clusters <- sort(as.numeric(as.character(unique(Idents(data)))),decreasing = F)

    groups <- sort(unique(plotx$GROUP))
    cluster_colors <- gen_colors(n = length(levels(plotx$CLUSTER)))
    names(cluster_colors) <- levels(plotx$CLUSTER)

    p1 <- NULL
    p2 <- NULL
    p1 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", title = project_name,isfacet = T,
                            col=cluster_colors, color_by = "CLUSTER", group_by = "GROUP",
                            xlabel = "UMAP_1", ylabel = "UMAP_2")
    p2 <- own_facet_scatter(plotx, "tSNE_1", "tSNE_2", title = project_name,isfacet = T,
                            col=cluster_colors, color_by = "CLUSTER", group_by = "GROUP",
                            xlabel = "UMAP_1", ylabel = "UMAP_2")

    somePDFPath = paste(cdir,"55URSA_PLOT_scRNASEQ_UMAP_AUTOCLUSTER_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=10,pointsize=12)
    print(p1)
    dev.off()

    somePDFPath = paste(cdir,"56URSA_PLOT_scRNASEQ_UMAP_AUTOCLUSTER_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=10,pointsize=12)
    print(p2)
    dev.off()

    Idents(data) <- integration_cluster
    DefaultAssay(data) <- "RNA"
    current_out <- deanalysis(data, current_clusters, plot_title = project_name,group=NULL,de_analysis = "findallmarkers", cluster_name = integration_cluster)
    de_type <- current_out$de_type
    de_name <- current_out$de_name
    top1 <- current_out$top1
    topn <- current_out$topn
    current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)
    current_data_markers <- current_out$current_data_markers
    current_data_markers <- current_data_markers[order(current_data_markers$avg_log2FC, decreasing = T),]
    write.table(current_data_markers, paste(cdir, "57URSA_TABLE_scRNASEQ_AUTO_CLUSTER_TOP_MARKERS_",integration_name,"_",project_name,".txt", sep = ""), quote = F, row.names = F, sep = "\t")

    somePDFPath = paste(cdir,"58URSA_PLOT_scRNASEQ_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=ceiling(length(top1$gene)/4)*5,pointsize=12)
    print(current_out$featureplot+plot_layout(ncol = 4) +
      plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5))))
    dev.off()

    somePDFPath = paste(cdir,"59URSA_PLOT_scRNASEQ_VIOLIN_TOP_",n,"_MARKERS_IN_CLUSTERS_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=12,pointsize=12)
    for(k in 1:length(current_clusters)){
      current <- topn[which(topn$cluster == current_clusters[k]),]
      current <- current[order(current$p_val_adj, decreasing = F),]
      print(VlnPlot(data, features = unique(current$gene), pt.size = 0, split.by = "GROUP",
                                                                          stack = T,flip = T,
                                                                          cols = cgroup_colors)&
        xlab("CLUSTERS")&
        plot_annotation(title = paste(project_name, ": CLUSTER ", current_clusters[k], sep = ""),
                        theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))
    }
    dev.off()

    somePDFPath = paste(cdir,"60URSA_PLOT_scRNASEQ_HEATMAP_TOP_",n,"_MARKERS_IN_CLUSTERS_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=16*length(current_clusters)/(7+2),pointsize=12)
    print(DoHeatmap(data, features = topn$gene,
                                       group.colors = cluster_colors, size = 8) +
      ggtitle(project_name)+
      NoLegend()+theme(axis.text.x = element_blank(),
                       axis.text.y = element_text(size = 10),
                       axis.title.x = element_blank(),
                       axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                       legend.title = element_blank(),
                       legend.text = element_text(size = 15),
                       plot.title = element_text(size = 25, face = "bold", hjust = 0.5)))
    dev.off()

    DefaultAssay(data) <- "RNA"
    orig_gene_names <- NULL
    scale_orig_gene_names <- NULL
    if(length(grep("ENSG[0-9]+",row.names(data), ignore.case = T)) > nrow(data)/2){
      orig_gene_names <- row.names(data)
      scale_orig_gene_names <- row.names(data@assays$RNA@scale.data)
      row.names(data@assays$RNA@counts) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data@assays$RNA@counts), ignore.case = T)
      row.names(data@assays$RNA@data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data@assays$RNA@data), ignore.case = T)
      row.names(data@assays$RNA@scale.data) <- gsub("^ENSG[0-9]+[-|\\||\\s+\\||\\t](.*)","\\1",row.names(data@assays$RNA@scale.data), ignore.case = T)
    }

    Idents(data) <- "seurat_clusters"
    DefaultAssay(data) <- "RNA"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
                       clusters =  data$seurat_clusters,
                       ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.fine)
    data$INTEGRATED_CELL_TYPE <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
    data@meta.data[which(is.na(data$INTEGRATED_CELL_TYPE)),"INTEGRATED_CELL_TYPE"] <- "Unidentifiable"

    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
                       clusters =  data$seurat_clusters,
                       ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.main)
    data$INTEGRATED_CELL_TYPE_LEVEL2 <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
    data@meta.data[which(is.na(data$INTEGRATED_CELL_TYPE_LEVEL2)),"CELL_TYPE"] <- "Unidentifiable"

    if(length(grep("ENSG[0-9]+",row.names(data), ignore.case = T)) > nrow(data)/2){
      row.names(data@assays$RNA@counts) <- orig_gene_names
      row.names(data@assays$RNA@data) <- orig_gene_names
      row.names(data@assays$RNA@scale.data) <- scale_orig_gene_names
    }

    Idents(data) <- data$INTEGRATED_CELL_TYPE
    plotx <- gen10x_plotx(data, include_meta = T)
    write.table(data.frame(plotx), paste(cdir, "61URSA_TABLE_scRNASEQ_AUTO_ANNOTATION_CELL_TYPE_META_",integration_name,"_",project_name,".txt", sep = ""), quote = F, row.names = T, sep = "\t")

    ct_colors <- gen_colors(n = length(unique(plotx$INTEGRATED_CELL_TYPE)))
    names(ct_colors) <- sort(unique((plotx$INTEGRATED_CELL_TYPE)))

    somePDFPath = paste(cdir,"62URSA_PLOT_scRNASEQ_UMAP_AUTO_CELL_TYPE_IDENTIFICATION_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=12,pointsize=12)
    print(plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "INTEGRATED_CELL_TYPE", plot_title = project_name,
                                          col = ct_colors, annot = TRUE, legend_position = "right", point_size = 1))
    dev.off()

    somePDFPath = paste(cdir,"63URSA_PLOT_scRNASEQ_TSNE_AUTO_CELL_TYPE_IDENTIFICATION_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=12,pointsize=12)
    print(plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "INTEGRATED_CELL_TYPE", plot_title = project_name,
                                          col = ct_colors, annot = TRUE, legend_position = "right", point_size = 1))
    dev.off()

    somePDFPath = paste(cdir,"64URSA_PLOT_scRNASEQ_UMAP_BY_GROUP_AUTO_CELL_TYPE_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=8,pointsize=12)
    print(own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", isfacet = T,
                                               color_by = "CELL_TYPE", xlabel = "UMAP_1", ylabel = "UMAP_2",
                                               group_by = "GROUP", title = project_name,
                                               col = ct_colors))
    dev.off()

    #################################### RECEPTOR-LIGAND #########################################################
    ligs <- as.character(unique(ramilowski_pairs$ligand))
    recs <- as.character(unique(ramilowski_pairs$receptor))
    ligs.present <- rownames(data)[rownames(data) %in% ligs]
    recs.present <- rownames(data)[rownames(data) %in% recs]
    genes.to.use <- union(ligs.present,recs.present)

    Idents(data) <- "INTEGRATED_CELL_TYPE"
    data_interactions <- celltalk(input_object=data,
                                  metadata_grouping="INTEGRATED_CELL_TYPE",
                                  ligand_receptor_pairs=ramilowski_pairs,
                                  number_cells_required=10,
                                  min_expression=10,
                                  max_expression=20000,
                                  scramble_times=100)

    top_stats <- data_interactions %>%
      filter(p_val<0.05) %>%
      group_by(cell_type1) %>%
      top_n(3,interact_ratio) %>%
      ungroup()

    somePDFPath = paste(cdir,"64URSA_PLOT_scRNASEQ_RECEPTOR_LIGAND_TOP_INTERACTIONS_INTEGRATED_DATA_CELL_TYPES_", integration_name, project_name, ".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=10,pointsize=10)
    circos.clear()
    circos_plot(ligand_receptor_frame=top_stats,
                cell_group_colors=ct_colors,
                ligand_color=color_conditions$general[3],
                receptor_color=color_conditions$general[1],
                cex_outer=1.6,
                cex_inner=0.8)
    lgd_celltypes = Legend(at = c("Receptor","Ligand"),
                           labels = c("Receptor","Ligand"),
                           type = "grid",
                           legend_gp = gpar(fill =c(color_conditions$general[c(1,3)])), title_position = "topleft",
                           title = "Legend")
    lgd_list_vertical = packLegend(lgd_celltypes)
    lgd_list_vertical
    circle_size = unit(1.25, "snpc")
    draw(lgd_list_vertical, x = circle_size, just = "left")
    dev.off()

    data_interactions <- data_interactions[which(data_interactions$p_val < 0.05),]
    data_interactions <- data_interactions[order(data_interactions$interact_ratio, decreasing = T),]
    write.table(data_interactions, paste(cdir, "64URSA_PLOT_scRNASEQ_RECEPTOR_LIGAND_P005_INTERACTIONS_INTEGRATED_DATA_CELL_TYPES_",integration_name,"_",project_name,".txt", sep = ""), quote = F, row.names = T, sep = "\t")

    #################################### MEDIAN EXPRESSION #########################################################

    current <- group_medianexpr(current_out$current_data_markers, data, ref_group = integration_cluster, group = "seurat_clusters", cell_type = F)
    plot_median <- current$plot_median
    top_markers <- current$top_markers
    plot_median[is.na(plot_median)] <- 0

    somePDFPath = paste(cdir,"65URSA_PLOT_scRNASEQ_HEATMAP_MEDIAN_EXPR_CLUSTER_TOP_FOLDCHANGE_SIG_GENES_", integration_name, project_name, ".pdf", sep = "")
    pdf(file=somePDFPath, width=8, height=10,pointsize=12)
    print(complex_heatmap(plot_median, col = color_conditions$BlueYellowRed))
    dev.off()

    #################################### TOP MARKERS + TSNE / UMAP ###################################################################
    current_clusters <- sort(as.numeric(as.character(unique(top_markers$cluster))), decreasing = F)

    p1 <- NULL
    p2 <- NULL
    p3 <- NULL
    p1 <- ggplot(data=data.frame(top_markers), aes(x=gene, y=data.frame(top_markers)[,wt], fill=cluster)) +
      geom_bar(stat="identity", position=position_dodge())+
      scale_fill_manual(values = cluster_colors)+
      theme_classic() +
      coord_flip()+
      ylab(wt)+
      facet_wrap(~cluster, scales = "free_y") +
      theme(axis.text.y=element_text(size = 10),
            strip.text.x = element_text(size = 15, face = "bold"),
            legend.title = element_text(size =20, face = "bold"),
            legend.text = element_text(size = 15))

    plotx <- gen10x_plotx(data, include_meta = T)
    plotx$CLUSTER <- data@meta.data[,integration_cluster]
    plotx$GROUP <- data$GROUP

    p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CLUSTER",
                       plot_title = project_name,col = cluster_colors,numeric = T,
                       annot = T, legend_position = "right", point_size = 1)
    p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CLUSTER",
                       plot_title = project_name,col = cluster_colors,numeric = T,
                       annot = T, legend_position = "right", point_size = 1)

    somePDFPath = paste(cdir,"66URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_TSNE_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=9,pointsize=12)
    grid.arrange(p1, p2, ncol = 2)
    dev.off()

    somePDFPath = paste(cdir,"67URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_UMAP_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=9,pointsize=12)
    grid.arrange(p1, p3, ncol = 2)
    dev.off()

    current_clusters <- sort(as.numeric(as.character(unique(top_markers$cluster))), decreasing = F)
    top_1_markers <- current_out$current_data_markers %>% group_by(cluster) %>% top_n(n = 1, eval(parse(text = wt)))

    Idents(data) <- "seurat_clusters"
    somePDFPath = paste(cdir,"68URSA_PLOT_scRNASEQ_RIDGE_TOP_FC_SIG_GENES_IN_CLUSTERS_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=18, height=length(unique(unlist(top_1_markers$gene)))/2*4,pointsize=12)
    print(RidgePlot(data, features = unique(unlist(top_1_markers$gene)), ncol = 4,
                                       cols = cluster_colors)+
      plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    current <- data@meta.data
    total_counts <- data.frame(table(current[,integration_cluster]))

    current <- data.frame(table(current[,c("SID",integration_cluster)]))
    current <- current[current$Freq != 0,]
    colnames(current) <- c("SID","CLUSTER","COUNT")
    current$CLUSTER_COUNT <- total_counts[match(current$CLUSTER, total_counts$Var1),"Freq"]
    current$PROPORTION <- current$COUNT/current$CLUSTER_COUNT

    node_proportion <- current

    somePDFPath = paste(cdir,"69URSA_PLOT_scRNASEQ_NODE_SUMMARY_SAMPLE_PROPORTION_",integration_name, "_",project_name,".pdf", sep = "")
    pdf(somePDFPath, width = 25, height = 20, pointsize = 10)
      print(ggplot(node_proportion, aes(CLUSTER, PROPORTION, fill = SID))+
      geom_bar(stat="identity", alpha=0.8)+
      coord_polar()+
      scale_fill_viridis(discrete = T, option = "C")+
      ggtitle(paste("Frequency of Samples in Each Node\n", project_name, sep = ""))+
      theme_bw(base_size = 28)+
      theme(plot.margin = unit(c(3,3,3,3), "cm"),
            plot.title = element_text(size=20, face = "bold", hjust = 0.5),
            legend.title=element_text(size=20, face = "bold"),
            legend.text=element_text(size=15),
            legend.key.size = unit(2, 'lines'),
            axis.title.x = element_text(colour="black", size = 20, face = "bold", vjust = -10),
            axis.title.y = element_text(colour="black", size = 20, face = "bold", vjust = 10),
            strip.text = element_text(size = 20, face = "bold"),
            axis.text.x=element_text(colour="black", size = 20),
            axis.text.y=element_text(colour="black", size = 20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()))
      dev.off()

    current <- data@meta.data
    total_counts <- data.frame(table(current$CELL_TYPE))

    current <- data.frame(table(current[,c("SID","CELL_TYPE")]))
    current <- current[current$Freq != 0,]
    colnames(current) <- c("SID","CELL_TYPE","COUNT")
    current$CELL_TYPE_COUNT <- total_counts[match(current$CELL_TYPE, total_counts$Var1),"Freq"]
    current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT

    somePDFPath = paste(cdir,"70URSA_PLOT_scRNASEQ_CELL_TYPE_SUMMARY_SAMPLE_PROPORTION_",integration_name, "_",project_name,".pdf", sep = "")
    pdf(somePDFPath, width = 12, height = 8, pointsize = 10)
    print(ggplot(current,aes(SID,PROPORTION,fill=CELL_TYPE))+
      geom_bar(stat = "identity", position = "fill", color = "black", size = 1)+ #, color = "black", size = 1
      theme_classic()+
      # facet_wrap(~GROUP)+
      theme(plot.margin = unit(c(1,1,1,1), "cm"),
            axis.text.x = element_text(angle = 45, size = 15, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1, "cm"),
            legend.position = "right",
            strip.text.x = element_text(size = 25),
            strip.text.y = element_text(size = 25),
            plot.title = element_text(size = 20, face=1, hjust = 0.5))+
      scale_fill_manual(values = ct_colors)+
      guides(color=guide_legend(title="CELL TYPES", ncol = 1),fill = guide_legend(title="CELL TYPES", ncol = 1)))
    dev.off()

    #####################################################################################################################
    current <- group_medianexpr(current_out$current_data_markers, data, group = "CELL_TYPE", cell_type = T)
    plot_median_cell_type <- current$plot_median
    top_markers_cell_type <- current$top_markers
    plot_median_cell_type[is.na(plot_median_cell_type)] <- 0
    somePDFPath = paste(cdir,"71URSA_PLOT_scRNASEQ_HEATMAP_MEDIAN_EXPR_CELL_TYPE_TOP_FOLDCHANGE_SIG_GENES_", integration_name, project_name, ".pdf", sep = "")
    pdf(file=somePDFPath, width=8, height=10,pointsize=12)
    print(complex_heatmap(plot_median_cell_type, col = color_conditions$BlueYellowRed))
    dev.off()

    p1 <- NULL
    p2 <- NULL
    p3 <- NULL
    cell_types <- sort(as.character(unique(top_markers_cell_type$CELL_TYPE)))

    p1 <- ggplot(data=data.frame(top_markers_cell_type), aes(x=gene, y=data.frame(top_markers_cell_type)[,wt], fill=CELL_TYPE)) +
      geom_bar(stat="identity", position=position_dodge())+
      scale_fill_manual(values = ct_colors)+
      theme_classic() +
      coord_flip()+
      ylab(wt)+
      facet_wrap(~CELL_TYPE, scales = "free_y") +
      theme(axis.text.y=element_text(size = 7))

    plotx <- gen10x_plotx(data, groups = NULL)
    plotx$CELL_TYPE <- data$CELL_TYPE

    p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "CELL_TYPE",
                       plot_title = integration_name,col = ct_colors,numeric = F,
                       annot = T, legend_position = "right", point_size = 1)
    p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "CELL_TYPE",
                       plot_title = integration_name,col = ct_colors,numeric = F,
                       annot = T, legend_position = "right", point_size = 1)

    somePDFPath = paste(cdir,"72URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_TSNE_CELL_TYPES_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=24, height=ceiling(length(unique(plotx$CELL_TYPE))/3)*5,pointsize=12)
    grid.arrange(p1, p2, ncol = 2)
    dev.off()

    somePDFPath = paste(cdir,"73URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_UMAP_CELL_TYPES_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=24, height=ceiling(length(unique(plotx$CELL_TYPE))/3)*5,pointsize=12)
    grid.arrange(p1, p3, ncol = 2)
    dev.off()

    current_out$current_data_markers$CELL_TYPE <- data@meta.data[match(current_out$current_data_markers$cluster, data@meta.data[,integration_cluster]),"CELL_TYPE"]
    top_1_markers <- current_out$current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = 1, eval(parse(text = wt)))
    Idents(data) <- "CELL_TYPE"

    somePDFPath = paste(cdir,"74URSA_PLOT_scRNASEQ_RIDGE_TOP_FC_SIG_GENES_CELL_TYPES_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=length(unique(unlist(top_1_markers$gene)))/2*6,pointsize=12)
    print(RidgePlot(data, features = unique(unlist(top_1_markers$gene)), ncol = 2,
                                       cols = ct_colors)+
      plot_annotation(title = integration_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    current <- data@meta.data
    total_counts <- data.frame(table(current$CELL_TYPE))

    current <- data.frame(table(current[,c("SID","CELL_TYPE")]))
    current <- current[current$Freq != 0,]
    colnames(current) <- c("SID","CELL_TYPE","COUNT")
    current$CELL_TYPE_COUNT <- total_counts[match(current$CELL_TYPE, total_counts$Var1),"Freq"]
    current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT

    node_proportion <- current

    somePDFPath = paste(cdir,"75URSA_PLOT_scRNASEQ_CELL_TYPE_SAMPLE_PROPORTION_",integration_name, "_",project_name,".pdf", sep = "")
    pdf(somePDFPath, width = 25, height = 20, pointsize = 10)
    print(ggplot(node_proportion, aes(CELL_TYPE, PROPORTION, fill = SID))+
      geom_bar(stat="identity", alpha=0.8)+
      coord_polar()+
      scale_fill_viridis(discrete = T)+
      ggtitle(paste("Frequency of Samples in Each Cell Type\n", project_name, sep = ""))+
      theme_bw(base_size = 25)+
      theme(plot.margin = unit(c(3,3,3,3), "cm"),
            plot.title = element_text(size=20, face = "bold", hjust = 0.5),
            legend.title=element_text(size=20, face = "bold"),
            legend.text=element_text(size=15),
            legend.key.size = unit(2, 'lines'),
            axis.title.x = element_text(colour="black", size = 20, face = "bold", vjust = -10),
            axis.title.y = element_text(colour="black", size = 20, face = "bold", vjust = 10),
            strip.text = element_text(size = 20, face = "bold"),
            axis.text.x=element_text(colour="black", size = 20),
            axis.text.y=element_text(colour="black", size = 20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()))
    dev.off()

    Idents(data) <- integration_cluster

    ############ GROUP COMPARISON #########################################################

    Idents(data) <- integration_cluster
    current_clusters <- sort(as.numeric(as.character(unique(Idents(data)))), decreasing = F)

    if(length(unique(data$GROUP)) > 1){
      current_out <- deanalysis(data, current_clusters, plot_title = project_name,group="GROUP",de_analysis = "finddegroup")
      current_data_markers <- current_out$current_data_markers
      current_data_markers <- current_data_markers[order(current_data_markers$avg_log2FC, decreasing = T),]
      write.table(current_data_markers, paste(cdir, "76URSA_TABLE_scRNASEQ_BY_CLUSTERS_TOP_MARKERS_DE_GROUPS_EACH_CLUSTER_",integration_name,"_",project_name,".txt", sep = ""), quote = F, row.names = F, sep = "\t")

      de_type <- current_out$de_type
      de_name <- current_out$de_name
      top1 <- current_out$top1
      topn <- current_out$topn
      current_clusters <- sort(as.numeric(as.character(unique(topn$cluster))), decreasing = F)

      somePNGPath <- paste(cdir,"76URSA_PLOT_scRNASEQ_CLUSTERS_BY_GROUP_TOP1_FEATURE_DIFFERENCE_COLOR_EXPR_",integration_name,"_",project_name,".png", sep = "")
      png(somePNGPath, width = 2000, height = length(current_clusters)*length(unique(data$GROUP))*50, units = "px", res = 150)
      print(DotPlot(data, features = unique(top1$gene), cols = cgroup_colors, dot.scale = 8,
                                       split.by = "GROUP", cluster.idents = F) + RotatedAxis()+ylab("CLUSTER+GROUP")+
        plot_annotation(title = project_name, theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))
      dev.off()

      somePDFPath = paste(cdir,"77URSA_PLOT_scRNASEQ_VIOLIN_TOP_",n,"_MARKERS_DE_GROUPS_EACH_CLUSTER_",integration_name,"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=16, height=10,pointsize=12)
      for(k in 1:length(current_clusters)){
        current <- topn[which(topn$cluster == current_clusters[k]),]
        current <- current[order(current$p_val_adj, decreasing = F),]
        print(VlnPlot(data, features = unique(current$gene), pt.size = 0, split.by = "GROUP",stack = T,flip = T,
                                                                            cols = cgroup_colors)&
          xlab("CLUSTERS")&
          plot_annotation(title = paste(project_name, ": CLUSTER ", current_clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))))
      }
      dev.off()

      somePDFPath = paste(cdir,"78URSA_PLOT_scRNASEQ_UMAP_DENSITY_TOP_1_MARKER_BY_DE_GROUP_EACH_CLUSTER_",integration_name,"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=5,pointsize=12)
      for(i in 1:length(current_out$featureplot)){
        print(current_out$featureplot[[i]])
        }
      dev.off()
    }
  }

  #################### PSEUDOTIME TRAJECTORY #######################################################
  mono3_current <- NULL
  sce_current <- NULL

  for(j in 1:length(data_current)){
    mono3_current[[j]] <- as.cell_data_set(data_current[[j]])

    ##### UMAP ######
    reducedDim(mono3_current[[j]], type = "PCA") <- data_current[[j]]@reductions$pca@cell.embeddings
    mono3_current[[j]]@reduce_dim_aux$prop_var_expl <- data_current[[j]]@reductions$pca@stdev
    mono3_current[[j]]@reduce_dim_aux$gene_loadings <- data_current[[j]]@reductions[["pca"]]@feature.loadings
    mono3_current[[j]]@int_colData@listData$reducedDims$UMAP <- data_current[[j]]@reductions$umap@cell.embeddings
    mono3_current[[j]]@clusters$UMAP$clusters <- data_current[[j]]$seurat_clusters
    mono3_current[[j]]@clusters@listData[["UMAP"]][["clusters"]] <- data_current[[j]]$seurat_clusters
    mono3_current[[j]]@clusters@listData[["UMAP"]]$cluster_result$optim_res$membership <- data_current[[j]]$seurat_clusters
    mono3_current[[j]]@clusters@listData[["UMAP"]][["louvain_res"]] <- NULL
    rownames(mono3_current[[j]]@principal_graph_aux$UMAP$dp_mst) <- NULL
    colnames(mono3_current[[j]]@int_colData@listData$reducedDims$UMAP) <- NULL
    mono3_current[[j]] <- cluster_cells(mono3_current[[j]], reduction_method = "UMAP", resolution = 1e-3)
    mono3_current[[j]]@clusters@listData[["UMAP"]][["clusters"]] <- data_current[[j]]$seurat_clusters
    mono3_current[[j]]@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
    mono3_current[[j]] <- learn_graph(mono3_current[[j]])

    ccolors <- gen_colors(n = length(unique(mono3_current[[j]]@clusters$UMAP$clusters)))
    names(ccolors) <- sort(as.numeric(as.character(unique(mono3_current[[j]]@clusters$UMAP$clusters))))

    p <- NULL
    somePDFPath = paste(cdir,"79URSA_PLOT_scRNASEQ_UMAP_TRAJECTORY_CLUSTERS_",annot_names[j],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=10,pointsize=12)
    p <- plot_pseudo(mono3_current[[j]], reduction = "UMAP",
                      group = "seurat_clusters", label_size = 7,
                      paste(annot_names[j],"\nTRAJECTORY - COLOR BY CLUSTERS"),
                      ccolors, length(ccolors))
    print(p)
    dev.off()

    p <- ggplotly(p+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                             guides(color=guide_legend(title="CLUSTER"))+
                                             ggtitle(paste("<b>",names(data_current)[j],"</b>",sep = "")))
    htmlwidgets::saveWidget(p, paste(cdir, "80URSA_PLOT_scRNASEQ_UMAP_TRAJECTORY_CLUSTERS_",annot_names[j],"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_pseudo(mono3_current[[j]], reduction = "UMAP",
                                              group = "CELL_TYPE", label_size = 6,
                                              paste(annot_names[j],"\nTRAJECTORY - COLOR BY CELL TYPE"),
                                              color_conditions$tenx, length(unique(mono3_current[[j]]$CELL_TYPE)))

    somePDFPath = paste(cdir,"81URSA_PLOT_scRNASEQ_UMAP_TRAJECTORY_CELL_TYPE_",annot_names[j],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=10,pointsize=12)
    print(p)
    dev.off()

    p <- ggplotly(p+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                             # guides(color=guide_legend(title=""))+
                                             ggtitle(paste("<b>",names(data_current)[j],"</b>",sep = "")))
    htmlwidgets::saveWidget(p, paste(cdir, "82URSA_PLOT_scRNASEQ_UMAP_TRAJECTORY_CELL_TYPE_",annot_names[j],"_",project_name,".html", sep = ""))

    current_clusters <- sort(as.numeric(as.character(unique(mono3_current[[j]]@clusters$UMAP$clusters))), decreasing = F)

    somePDFPath = paste(cdir,"83URSA_PLOT_scRNASEQ_UMAP_PSEUDOTIME_EACH_CLUSTER_AS_ROOT_", annot_names[j],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=10,pointsize=12)
    for(i in 1:length(current_clusters)){
      mono3_current[[j]] <- order_cells(mono3_current[[j]],
                                        root_pr_nodes=get_earliest_principal_node(mono3_current[[j]],
                                                                                  group = "seurat_clusters", group_element = current_clusters[i]))
      p <- NULL
      p <- plot_cells(mono3_current[[j]],
                                                                             color_cells_by = "pseudotime",
                                                                             group_label_size = 7,
                                                                             graph_label_size = 5,
                                                                             cell_size = 1,
                                                                             cell_stroke = I(2/2),
                                                                             alpha = 1,
                                                                             trajectory_graph_segment_size = 1.5,
                                                                             label_cell_groups=FALSE,
                                                                             label_leaves=FALSE,
                                                                             label_branch_points=FALSE)

      p <- adjust_plot(p, col = color_conditions$colorful, n = length(current_clusters),
                                                                            plot_title = paste(annot_names[j],"\nPSEUDOTIME - ROOT: CLUSTER ", current_clusters[i], sep = ""),
                                                                            fill = T)
      print(p)

    }
dev.off()

pseudo_para <- c("PSEUDOTIME_PC1","PSEUDOTIME_UMAP1","PSEUDOTIME_tSNE1")
chosen_para <- c("seurat_clusters", "CELL_TYPE")
chosen_para_names <- c("CLUSTERS", "CELL_TYPE")

    plotx <- data.frame(
      PSEUDOTIME_PC1 = rank(data_current[[j]]@reductions$pca@cell.embeddings[,"PC_1"]),
      PSEUDOTIME_UMAP1 = rank(data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_1"]),
      PSEUDOTIME_tSNE1 = rank(data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_1"]),
      data_current[[j]]@meta.data)

    somePDFPath = paste(cdir,"84URSA_PLOT_scRNASEQ_PSEUDOTIME_RANKED_PCA_UMAP_TSNE_",annot_names[j],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=8, height=6,pointsize=12)
    for(m in 1:length(chosen_para)){
      for(n in 1:length(pseudo_para)){
        p <- NULL
        p <- ggplot(plotx, aes_string(x = pseudo_para[n], y = chosen_para[m],
                          colour = chosen_para[m])) +
          geom_quasirandom(groupOnX = FALSE) +
          scale_color_manual(values = color_conditions$manycolors) + theme_classic() +
          xlab(pseudo_para[n]) + ylab(chosen_para_names[m]) +
          ggtitle(annot_names[j])
        print(p)
      }
    }
dev.off()

somePDFPath = paste(cdir,"85URSA_PLOT_scRNASEQ_UMAP_PSEUDOTIME_EVERY_CELL_TYPE_AS_ROOT_", annot_names[j],"_",project_name,".pdf", sep = "")
pdf(file=somePDFPath, width=12, height=10,pointsize=12)
    current_cell_types <- sort(as.character(unique(colData(mono3_current[[j]])[,"CELL_TYPE"])))
    for(i in 1:length(current_cell_types)){
      mono3_current[[j]] <- order_cells(mono3_current[[j]],
                                        root_pr_nodes=get_earliest_principal_node(mono3_current[[j]],
                                                                                  group = "CELL_TYPE", group_element = current_cell_types[i]))

      p <- NULL
      p <- plot_cells(mono3_current[[j]],
                                                                             color_cells_by = "pseudotime",
                                                                             group_label_size = 7,
                                                                             graph_label_size = 5,
                                                                             cell_size = 2,
                                                                             cell_stroke = I(2/2),
                                                                             alpha = 0.7,
                                                                             trajectory_graph_segment_size = 1.5,
                                                                             label_cell_groups=FALSE,
                                                                             label_leaves=FALSE,
                                                                             label_branch_points=FALSE)

      p <- adjust_plot(p, col = color_conditions$tenx, n = length(current_cell_types),
                                                                            plot_title = paste(annot_names[j],"\nPSEUDOTIME - ROOT CELL TYPE: ", toupper(current_cell_types[i]), sep = ""),
                                                                            fill = T)
      print(p)
    }
    dev.off()

    selected_pcs <- NULL
    if(class(overall_pcs) == "numeric"){
      selected_pcs <- overall_pcs
    }else if (toupper(overall_pcs) == toupper("all")){
      selected_pcs <- findPC(sdev = data_current[[j]]@reductions$pca@stdev, number = 50, method = 'all',aggregate = 'mean')
    }else if(toupper(overall_pcs) %in% toupper(c('piecewise linear model', 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line', 'k-means clustering'))){
      selected_pcs <- findPC(sdev = data_current[[j]]@reductions$pca@stdev, number = 50, method = tolower(overall_pcs))
    }

    data_current[[j]] <- RunUMAP(data_current[[j]], n.components = 3, reduction = "pca", dims = 1:selected_pcs)

    p <- NULL
    p <- plot_3d(data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                                          data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_2"],
                                          data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_3"],
                                          color_group = data_current[[j]]@meta.data[,"CELL_TYPE"],
                                          plot_title = names(data_current)[j],
                                          col = color_conditions$tenx,
                                          n = length(unique(data_current[[j]]$seurat_clusters)),
                                          lt = "CELL TYPE", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                          current_text= paste("Cluster ",data_current[[j]]$seurat_clusters,"\n",
                                                              "Cell: ",row.names(data_current[[j]]@meta.data),"\n",
                                                              data_current[[j]]@meta.data[,"CELL_TYPE"], sep = ""))

    htmlwidgets::saveWidget(p, paste(cdir, "86URSA_PLOT_scRNASEQ_UMAP_INTERACTIVE_CELL_TYPES_",annot_names[j],"_",project_name,".html", sep = ""))

    data_current[[j]] <- RunTSNE(data_current[[j]], dim.embed = 3, reduction = "pca", dims = 1:selected_pcs, check_duplicates = FALSE)

    p <- NULL
    p <- plot_3d(data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                          data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                          data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                          color_group = data_current[[j]]@meta.data[,"CELL_TYPE"],
                                          plot_title = names(data_current)[j],
                                          col = color_conditions$tenx,
                                          n = length(unique(data_current[[j]]$seurat_clusters)),
                                          lt = "CELL TYPE", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                          current_text= paste("Cluster ",data_current[[j]]$seurat_clusters,"\n",
                                                              "Cell: ",row.names(data_current[[j]]@meta.data),"\n",
                                                              data_current[[j]]@meta.data[,"CELL_TYPE"], sep = ""))

    htmlwidgets::saveWidget(p, paste(cdir, "87URSA_PLOT_scRNASEQ_TSNE_INTERACTIVE_CELL_TYPES_",annot_names[j],"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_3d(data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_1"],
                                          data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_2"],
                                          data_current[[j]]@reductions$umap@cell.embeddings[,"UMAP_3"],
                                          color_group = data_current[[j]]$seurat_clusters,
                                          plot_title = names(data_current)[j],
                                          col = color_conditions$colorful,
                                          n = length(unique(data_current[[j]]$seurat_clusters)),
                                          lt = "CLUSTER", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                          current_text= paste("Cell: ",row.names(data_current[[j]]@meta.data),"\nCluster: ",
                                                              data_current[[j]]$seurat_clusters, sep = ""))
    htmlwidgets::saveWidget(p, paste(cdir, "88URSA_PLOT_scRNASEQ_UMAP_INTERACTIVE_CLUSTERS_",annot_names[j],"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_3d(data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                          data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                          data_current[[j]]@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                          color_group = data_current[[j]]$seurat_clusters,
                                          plot_title = names(data_current)[j],
                                          col = color_conditions$colorful,
                                          n = length(unique(data_current[[j]]$seurat_clusters)),
                                          lt = "CLUSTER", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                          current_text= paste("Cell: ",row.names(data_current[[j]]@meta.data),"\nCluster: ",
                                                              data_current[[j]]$seurat_clusters, sep = ""))
    htmlwidgets::saveWidget(p, paste(cdir, "89URSA_PLOT_scRNASEQ_TSNE_INTERACTIVE_CLUSTERS_",annot_names[j],"_",project_name,".html", sep = ""))

    data_current[[j]] <- RunTSNE(data_current[[j]], reduction = "pca", dims = 1:selected_pcs, check_duplicates = FALSE)
    data_current[[j]] <- RunUMAP(data_current[[j]], reduction = "pca", dims = 1:selected_pcs)
  }

  if(length(data_current) > 1){

    pseudo_para <- c("PSEUDOTIME_PC1","PSEUDOTIME_UMAP1","PSEUDOTIME_tSNE1")
    chosen_para <- c(integration_cluster, "CELL_TYPE","orig.ident")
    chosen_para_names <- c("CLUSTERS", "CELL_TYPE","SAMPLE_ID")
    plotx <- data.frame(
      PSEUDOTIME_PC1 = rank(data@reductions$pca@cell.embeddings[,"PC_1"]),
      PSEUDOTIME_UMAP1 = rank(data@reductions$umap@cell.embeddings[,"UMAP_1"]),
      PSEUDOTIME_tSNE1 = rank(data@reductions$tsne@cell.embeddings[,"tSNE_1"]),
      data@meta.data)

    somePDFPath = paste(cdir,"90URSA_PLOT_scRNASEQ_INTEGRATED_DATA_PSEUDOTIME_RANKED_PCA_UMAP_TSNE_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=8, height=6,pointsize=12)
    for(m in 1:length(chosen_para)){
      for(n in 1:length(pseudo_para)){
        p <- NULL
        p <- ggplot(plotx, aes_string(x = pseudo_para[n], y = chosen_para[m],
                                      colour = chosen_para[m])) +
          geom_quasirandom(groupOnX = FALSE) +
          scale_color_manual(values = color_conditions$manycolors) + theme_classic() +
          xlab(pseudo_para[n]) + ylab(chosen_para_names[m]) +
          ggtitle(project_name)
        print(p)
      }
    }
    dev.off()

    if(integration_name == "SEURAT_INTEGRATED"){
      DefaultAssay(data) <- "integrated"
      mono3data <- as.cell_data_set(data)
    }else if(integration_name == "HARMONY_INTEGRATED"){
      DefaultAssay(data) <- "RNA"
      mono3data <- as.cell_data_set(data)
    }

    ##### UMAP ######
    named_clusters <- data@meta.data[,integration_cluster]
    names(named_clusters) <- row.names(data@meta.data)
    reducedDim(mono3data, type = "PCA") <- data@reductions[[reduction_method]]@cell.embeddings
    mono3data@reduce_dim_aux$prop_var_expl <- data@reductions$pca@stdev
    mono3data@reduce_dim_aux$gene_loadings <- data@reductions[[reduction_method]]@feature.loadings
    mono3data@int_colData@listData$reducedDims$UMAP <- data@reductions$umap@cell.embeddings
    mono3data@clusters$UMAP$clusters <- named_clusters
    mono3data@clusters@listData[["UMAP"]][["clusters"]] <- named_clusters
    mono3data@clusters@listData[["UMAP"]]$cluster_result$optim_res$membership <- named_clusters
    mono3data@clusters@listData[["UMAP"]][["louvain_res"]] <- NULL
    rownames(mono3data@principal_graph_aux$UMAP$dp_mst) <- NULL
    colnames(mono3data@int_colData@listData$reducedDims$UMAP) <- NULL
    mono3data <- cluster_cells(mono3data, reduction_method = "UMAP", resolution = 1e-3)
    mono3data@clusters@listData[["UMAP"]][["clusters"]] <- named_clusters
    mono3data@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
    mono3data <- learn_graph(mono3data)

    p <- NULL
    p <- plot_pseudo(mono3data, reduction = "UMAP",
                                         group = integration_cluster, label_size = 7,
                                         paste(project_name,"\nINTEGRATED TRAJECTORY - COLOR BY CLUSTERS"),
                                         color_conditions$colorful, length(unique(mono3data@clusters$UMAP$clusters)))+facet_wrap(~GROUP)

    somePDFPath = paste(cdir,"91URSA_PLOT_scRNASEQ_INTEGRATED_DATA_UMAP_TRAJECTORY_CLUSTERS_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=9*length(unique(mono3data@colData$GROUP)), height=8,pointsize=12)
    print(p)
    dev.off()

    p <- ggplotly(p+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                        guides(color=guide_legend(title="CLUSTER"))+
                                        ggtitle(paste("<b>",project_name,"</b>",sep = "")))

    htmlwidgets::saveWidget(p, paste(cdir, "92URSA_PLOT_scRNASEQ_INTEGRATED_DATA_UMAP_TRAJECTORY_CLUSTERS_",integration_name,"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_pseudo(mono3data, reduction = "UMAP",
                                         group = "CELL_TYPE", label_size = 6,
                                         # cell_size = 0.5,traj_size = 0.8,
                                         paste(project_name,"\nINTEGRATED TRAJECTORY - COLOR BY CELL TYPE"),
                                         color_conditions$tenx, length(unique(mono3data$CELL_TYPE)))+
      facet_wrap(~GROUP)

    somePDFPath = paste(cdir,"93URSA_PLOT_scRNASEQ_INTEGRATED_DATA_UMAP_TRAJECTORY_CELL_TYPE_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=9*length(unique(mono3data@colData$GROUP)), height=8,pointsize=12)
    print(p)
    dev.off()

    p <- ggplotly(p+theme(legend.position = "right", plot.margin=unit(c(1,5,1,1),"cm"))+
                                        ggtitle(paste("<b>",project_name,"</b>",sep = "")))
    htmlwidgets::saveWidget(p, paste(cdir, "94URSA_PLOT_scRNASEQ_INTEGRATED_DATA_UMAP_TRAJECTORY_CELL_TYPE_",integration_name,"_",project_name,".html", sep = ""))

    current_clusters <- sort(as.numeric(as.character(unique(mono3data@clusters$UMAP$clusters))), decreasing = F)

    somePDFPath = paste(cdir,"95URSA_PLOT_scRNASEQ_INTEGRATED_DATA_UMAP_PSEUDOTIME_EVERY_CLUSTER_AS_ROOT_", integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=9*length(unique(mono3data@colData$GROUP)), height=8,pointsize=12)
    for(i in 1:length(current_clusters)){
      mono3data <- order_cells(mono3data,
                               root_pr_nodes=get_earliest_principal_node(mono3data,
                                                                         group = "seurat_clusters", group_element = current_clusters[i]))

      p <- plot_cells(mono3data,
                                                                             color_cells_by = "pseudotime",
                                                                             group_label_size = 7,
                                                                             graph_label_size = 5,
                                                                             cell_size = 1,
                                                                             cell_stroke = I(2/2),
                                                                             alpha = 0.7,
                                                                             trajectory_graph_segment_size = 1.5,
                                                                             label_cell_groups=FALSE,
                                                                             label_leaves=FALSE,
                                                                             label_branch_points=FALSE)+
        facet_wrap(~GROUP)

      p <- adjust_plot(p, col = color_conditions$colorful, n = length(current_clusters),
                                                                            plot_title = paste(project_name,"\nPSEUDOTIME - ROOT: CLUSTER ", current_clusters[i], sep = ""),
                                                                            fill = T)
    print(p)
    }
    dev.off()

    if(toupper(reduction_method) == "HARMONY"){
      DefaultAssay(data) <- 'RNA'
      selected_pcs <- NULL
      if(class(overall_pcs) == "numeric"){
        selected_pcs <- overall_pcs
      }else if (toupper(overall_pcs) == toupper("all")){
        selected_pcs <- findPC(sdev = sort(data@reductions$harmony@stdev, decreasing = T), number = 50, method = 'all',aggregate = 'mean')
      }else if(toupper(overall_pcs) %in% toupper(c('piecewise linear model', 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line', 'k-means clustering'))){
        selected_pcs <- findPC(sdev = sort(data@reductions$harmony@stdev, decreasing = T), number = 50, method = tolower(overall_pcs))
      }

      data <- RunTSNE(data, dim.embed = 3, reduction = "harmony", dims = 1:selected_pcs, check_duplicates = FALSE)
      data <- RunUMAP(data, n.components = 3, reduction = "harmony", dims = 1:selected_pcs)

    }else if(reduction_method == "pca"){
      DefaultAssay(data) <- 'integrated'
      data <- RunTSNE(data, dim.embed = 3, reduction = "pca", dims = 1:selected_pcs, check_duplicates = FALSE)
      data <- RunUMAP(data, n.components = 3, reduction = "pca", dims = 1:selected_pcs)
    }

    p <- NULL
    p <- plot_3d(data@reductions$umap@cell.embeddings[,"UMAP_1"],
                                     data@reductions$umap@cell.embeddings[,"UMAP_2"],
                                     data@reductions$umap@cell.embeddings[,"UMAP_3"],
                                     color_group = data@meta.data[,"CELL_TYPE"],
                                     plot_title = names(data_current)[j],
                                     col = color_conditions$tenx,
                                     n = length(unique(data$CELL_TYPE)),
                                     lt = "CELL TYPE", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                     current_text= paste("Cluster ",data@meta.data[,integration_cluster],"\n",
                                                         "Cell: ",row.names(data@meta.data),"\n",
                                                         data@meta.data[,"CELL_TYPE"], sep = ""))

    htmlwidgets::saveWidget(p, paste(cdir, "96URSA_PLOT_scRNASEQ_UMAP_INTERACTIVE_CELL_TYPES_",integration_name,"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_3d(data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                     data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                     data@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                     color_group = data@meta.data[,"CELL_TYPE"],
                                     plot_title = names(data_current)[j],
                                     col = color_conditions$tenx,
                                     n = length(unique(data$CELL_TYPE)),
                                     lt = "CELL TYPE", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                     current_text= paste("Cluster ",data@meta.data[,integration_cluster],"\n",
                                                         "Cell: ",row.names(data@meta.data),"\n",
                                                         data@meta.data[,"CELL_TYPE"], sep = ""))
    htmlwidgets::saveWidget(p, paste(cdir, "97URSA_PLOT_scRNASEQ_TSNE_INTERACTIVE_CELL_TYPES_",integration_name,"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_3d(data@reductions$umap@cell.embeddings[,"UMAP_1"],
                                     data@reductions$umap@cell.embeddings[,"UMAP_2"],
                                     data@reductions$umap@cell.embeddings[,"UMAP_3"],
                                     color_group = data@meta.data[,integration_cluster],
                                     plot_title = names(data_current)[j],
                                     col = color_conditions$colorful,
                                     n = length(unique(data@meta.data[,integration_cluster])),
                                     lt = "CLUSTER", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                     current_text= paste("Cell: ",row.names(data@meta.data),"\nCluster: ",
                                                         data@meta.data[,integration_cluster], sep = ""))

    htmlwidgets::saveWidget(p, paste(cdir, "98URSA_PLOT_scRNASEQ_UMAP_INTERACTIVE_CLUSTERS_",integration_name,"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_3d(data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                     data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                     data@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                     color_group = data@meta.data[,integration_cluster],
                                     plot_title = names(data_current)[j],
                                     col = color_conditions$colorful,
                                     n = length(unique(data@meta.data[,integration_cluster])),
                                     lt = "CLUSTER", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                     current_text= paste("Cell: ",row.names(data@meta.data),"\nCluster: ",
                                                         data@meta.data[,integration_cluster], sep = ""))
    htmlwidgets::saveWidget(p, paste(cdir, "99URSA_PLOT_scRNASEQ_TSNE_INTERACTIVE_CLUSTERS_",integration_name,"_",project_name,".html", sep = ""))

  }
  ##########################################################################################
  ##########################################################################################

  print("Completed!")
  saveRDS(data, paste(cdir, "100URSA_DATA_scRNASEQ_",project_name, ".RDS", sep = ""))

}
