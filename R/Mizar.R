############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Mizar: Spatial
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Date: 2022-02-16
############################################################################################################
#' @include ini.R
#' @include common.R
#' @import ComplexHeatmap
#' @import cowplot
#' @import data.table
#' @import dplyr
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
#' @import celldex
#' @import Seurat
#' @import SingleR
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format Spatial Transcriptomics pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification Spatial
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_Spatial' by default.
#' @param input_dir Directory to all input files. Current working directory by
#' default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format
#' files.
#' @param rnaseq_dir scRNA-seq data files directory. We recommend putting the scRNA-seq
#' data in a separate directory separating them from spatial input files.
#' scRNA-seq files should be filtered .h5 format files from 10X Genomics and their file
#' names should have the same sample ID as file prefix to their correspond spatial files.
#' @param run_rnaseq If scRNA-seq integration with spatial data should be run. Default to FALSE.
#' There is no need to prepare scRNA-seq folder if run_rnaseq is set to FALSE.
#' @export
#'
SpatialPip <- function(project_name = "Ursa_Spatial",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file,
                     rnaseq_dir = "./",
                     run_rnaseq = FALSE){
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "SPATIAL", isDir = T)
  color_conditions <- color_ini()
  ctime <- time_ini()
  hpca.se <- HumanPrimaryCellAtlasData()
  hs <- org.Hs.eg.db
  hgnc.table <- data("hgnc.table", package="HGNChelper")
  p_val_adj <- 0.1
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir)

  for(i in 1:nrow(pheno_data)){
    if(length(grep("\\.csv$|\\.txt$", pheno_data$FILE[i], ignore.case = T)) == 0){
      if(length(grep("\\.h5$", pheno_data$FILE[i], ignore.case = T)) == 0){
        pheno_data$FILE[i] <- paste(pheno_data$FILE[i],".h5", sep = "")
      }
    }
  }

  sample_files <- list.files(input_dir, recursive = T, full.names = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep(gsub(".*\\/(.*)","\\1",pheno_file, ignore.case = T),sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep("\\.h5$", sample_files, ignore.case = T)]

  data <- NULL
  sample_features <- NULL
  data_current <- NULL
  data_markers <- NULL
  plotx <- NULL
  results <- NULL
  annot_names <- NULL
  scrna_data <- NULL

  for(i in 1:length(sample_files)){
    current_sample <- gsub(".*\\/(.*\\.h5)$","\\1",sample_files[i], ignore.case = T)
    cpheno <- pheno_data[which(toupper(pheno_data$FILE) == toupper(current_sample)),]
    cname <- cpheno$SAMPLE_ID
    annot_names <- c(annot_names, cname)
    print(current_sample)
    print(paste("Running sample: ", current_sample, "..", sep = ""))
    current <- Load10X_Spatial(data.dir = gsub("(.*\\/).*","\\1",sample_files[i], ignore.case = T), filename = current_sample)
    current <- SCTransform(current, assay = "Spatial", verbose = FALSE)
    current@project.name <- project_name
    current$orig.ident <- cname
    names(current@images) <- cname
    # current@images[[1]]@key <- cname
    # current@images[[1]]@key <- paste(gsub("\\.","",cname),"_", sep = "")
    colnames(current@meta.data)[grep("nCount_Spatial", colnames(current@meta.data), ignore.case = T)] <- "Spatial_Counts"
    colnames(current@meta.data)[grep("nFeature_Spatial", colnames(current@meta.data), ignore.case = T)] <- "Feature_Counts"
    levels(current@active.ident) <- cname
    current <- add_names(current, cname, cname)
    current$CELL_ID <- row.names(current@meta.data)
    print(current)
    DefaultAssay(current) <- "SCT"
    current <- FindVariableFeatures(current)
    sample_features[[i]] <- VariableFeatures(current)
    names(sample_features)[i] <- cname

    plotx <- current
    plotx$orig.ident <- ifelse(nchar(plotx$orig.ident) > 15, substr(plotx$orig.ident, 1, 15), plotx$orig.ident)
    p1 <- NULL
    p2 <- NULL
    p1 <- VlnPlot(plotx, group.by = "orig.ident", features = "Spatial_Counts", pt.size = 0.1, cols = color_conditions$general) + NoLegend() + xlab("SAMPLE_NAME")
    p2 <- SpatialFeaturePlot(plotx, features = "Spatial_Counts") + theme(legend.position = "right")
    somePDFPath = paste(cdir,"1URSA_PLOT_VIOLIN_IMAGE_FEATURES_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
    print((p1|p2)+plot_annotation(title = gsub("\\."," - ",annot_names[i]), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    DefaultAssay(plotx) <- "SCT"
    p1 <- SpatialFeaturePlot(plotx, features = head(VariableFeatures(plotx), 5), ncol = 5)
    p2 <- SpatialFeaturePlot(plotx, features = head(VariableFeatures(plotx), 5), alpha = c(0.1, 1), ncol = 5)
    somePDFPath = paste(cdir,"2URSA_PLOT_SPATIALLY_VARIABLE_FEATURES_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
    print(p1/p2)
    dev.off()

    print(paste("Running dimension reduction and clustering..", sep = ""))
    DefaultAssay(current) <- "SCT"
    current <- RunPCA(current, assay = "SCT", verbose = FALSE)
    current <- RunUMAP(current, reduction = "pca", dims = 1:30)
    current <- RunTSNE(current, reduction = "pca", dims = 1:30, check_duplicates = FALSE)
    current <- FindNeighbors(current, reduction = "pca", dims = 1:30)
    current <- FindClusters(current, verbose = FALSE)

    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "SCT"
    current_de_markers <- NULL
    current_de_markers <- FindAllMarkers(current, min.pct = 0.25, logfc.threshold = 0.25)
    current_de_markers <- data.frame(SAMPLE = annot_names[i], current_de_markers)
    current_de_markers <- current_de_markers[which(current_de_markers$p_val_adj < 0.1),]
    current_de_markers <- current_de_markers[order(current_de_markers$avg_log2FC, decreasing = T),]

    data_markers[[i]] <- current_de_markers
    names(data_markers)[i] <- annot_names[i]

    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "SCT"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(current)),
                       clusters =  current$seurat_clusters,
                       ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.fine)
    current$Cell_Type <- clu_ann$labels[match(current$seurat_clusters,row.names(clu_ann))]
    current@meta.data[which(is.na(current$Cell_Type)),"Cell_Type"] <- "Unidentifiable"
    Idents(current) <- "Cell_Type"

    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(current$seurat_clusters)))
    names(cluster_colors) <- sort(unique(current$seurat_clusters), decreasing = F)
    ct_colors <- gen_colors(color_conditions$colorful, length(unique(current$Cell_Type)))
    names(ct_colors) <- sort(unique(current$Cell_Type))
    sample_colors <- gen_colors(color_conditions$bright, length(unique(current$orig.ident)))
    names(sample_colors) <- sort(unique(current$orig.ident), decreasing = F)

    plotx <- gen10x_plotx(current, include_meta = T)
    plotx$CLUSTERS <- plotx$seurat_clusters

    p <- NULL
    p <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTERS", point_size = 1, label_size = 6, numeric = T, annot = T,
                       plot_title = paste(project_name, ": UMAP CLUSTERS", sep = ""), col = cluster_colors)
    p <- adjust_theme(p, legend_text = 10, legend_title = 15)

    somePDFPath = paste(cdir,"3URSA_PLOT_UMAP_ALL_SAMPLES_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
    print(p)
    dev.off()

    p <- NULL
    p <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "CLUSTERS", point_size = 1, label_size = 6, numeric = T, annot = T,
                      plot_title = paste(project_name, ": tSNE CLUSTERS", sep = ""), col = cluster_colors)
    p <- adjust_theme(p, legend_text = 10, legend_title = 15)

    somePDFPath = paste(cdir,"4URSA_PLOT_tSNE_ALL_SAMPLES_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
    print(p)
    dev.off()

    Idents(current) <- "seurat_clusters"
    p <- NULL
    p <- SpatialDimPlot(current,ncol = 1,pt.size.factor = 1.6,
                        images = annot_names[i], cols = cluster_colors)+
      ggtitle(annot_names[i]) +
      labs(fill= "CLUSTERS")+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.size = unit(0.5, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(fill=guide_legend(override.aes = list(size = 5)))

    somePDFPath = paste(cdir,"5URSA_PLOT_SLIDE_CLUSTERS_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
    print(p)
    dev.off()

    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "SCT"
    somePDFPath = paste(cdir,"6URSA_PLOT_SLIDE_IMAGE_BY_EACH_CLUSTER_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=7,pointsize=12)
    print(SpatialDimPlot(current, cells.highlight = CellsByIdentities(object = current,
                                                                      idents = sort(unique(current$seurat_clusters))), facet.highlight = TRUE, ncol = 5))
    dev.off()

    Idents(current) <- "seurat_clusters"
    plotx <- data.frame(Cluster = current$seurat_clusters,
                        X = current@images[[annot_names[i]]]@coordinates$imagerow,
                        Y = current@images[[annot_names[i]]]@coordinates$imagecol)
    p1 <- NULL
    p2 <- NULL
    p1 <- plot_slice(plotx, pt_size = 2, plot_title = "SLICE CLUSTER VIEW", col = cluster_colors, is.facet = FALSE, annot = FALSE)
    p1 <- p1+theme(legend.title = element_text(size = 15),
                   legend.text = element_text(size = 15),
                   legend.key.size = unit(0.2, "cm"))
    p2 <- plot_slice(plotx, pt_size = 1, plot_title = "SLICE CLUSTER VIEW", col = cluster_colors, is.facet = TRUE, annot = FALSE, strip_size = 15)

    somePDFPath = paste(cdir,"8URSA_PLOT_UMAP_SAMPLE_IMAGE_CLUSTER_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=5,pointsize=12)
    print((p1|p2)+plot_annotation(title = annot_names[i], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    somePDFPath = paste(cdir,"9URSA_PLOT_SPATIALLY_TOP_FEATURES_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=5,pointsize=12)
    cclusters <- sort(as.numeric(as.character(unique(current_de_markers$cluster))))
    for(j in 1:length(cclusters)){
      ccurrent_de <- current_de_markers[which(current_de_markers$cluster == cclusters[j]),]
      ccurrent_de <- ccurrent_de[order(ccurrent_de$avg_log2FC, decreasing = T),]
      if(nrow(ccurrent_de) > 0){
        ccurrent_de <- ccurrent_de[1:ifelse(nrow(ccurrent_de) < 5, nrow(ccurrent_de), 5),]
        p <- NULL
        p <- SpatialFeaturePlot(object = current, features = ccurrent_de$gene, images = annot_names[i], ncol = ifelse(nrow(ccurrent_de) < 5, nrow(ccurrent_de), 5))+plot_annotation(title = paste(annot_names[i],": CLUSTER ", cclusters[j], sep = ""))
        print(p)
      }
    }
    dev.off()

    write.table(current_de_markers, paste(cdir,"10URSA_PLOT_TABLE_DATA_CLUSTER_DEG_",annot_names[i],"_",project_name, ".txt", sep = ""), row.names = F, quote = F, sep = "\t")

    current_de <- current_de_markers[current_de_markers$p_val_adj < 0.1,]
    current_de <- current_de[order(current_de$avg_log2FC, decreasing = T),]
    current_de <- split(current_de, current_de$cluster)
    top_markers <- NULL
    for(j in 1:length(current_de)){
      if(nrow(current_de[[j]]) > 0){
        top_markers <- rbind(top_markers,current_de[[j]][1,])
      }
    }

    DefaultAssay(current) <- "Spatial"
    p <- NULL
    p <- VlnPlot(current, group.by = "seurat_clusters", features = unique(top_markers$gene),
                 pt.size = 0, ncol = 2, cols = cluster_colors, log = T)
    p <- p+plot_annotation(title = paste("TOP GENE IN EACH CLUSTER\nLog(Average Expression)", sep = ""), theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

    somePDFPath = paste(cdir,"11URSA_PLOT_VIOLINPLOT_TOP_MARKERS_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=10, height=2*ceiling(length(unique(top_markers$gene))/2),pointsize=12)
    print(p&xlab("CLUSTERS"))
    dev.off()

    Idents(current) <- current$Cell_Type
    plotx <- gen10x_plotx(current, include_meta = T)
    p1 <- NULL
    p2 <- NULL
    p <- NULL
    p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "Cell_Type", plot_title = "SPATIAL UMAP - CELL TYPE",
                       col = ct_colors, annot = TRUE, legend_position = "right", point_size = 1)
    p2 <- SpatialDimPlot(current,ncol = 1,pt.size.factor = 1.6,images = annot_names[i],
                         cols = ct_colors)+
      ggtitle("SPATIAL SLIDE - CELL TYPE") +
      labs(fill= "CELL TYPES")+
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1, "cm"))+
      guides(fill=guide_legend(override.aes = list(size = 5)))
    p <- (p1/p2)+plot_annotation(title = paste(project_name,": ",annot_names[i],sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

    somePDFPath = paste(cdir,"12URSA_PLOT_UMAP_AUTO_CELLTYPE_ANNOTATIONS_",annot_names[i],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=15, height=10,pointsize=12)
    print(p)
    dev.off()
    data_current[[i]] <- current
    names(data_current)[i] <- cname
}
    if(run_rnaseq == TRUE){
    scrna_files <- list.files(rnaseq_dir, pattern = "\\.h5$", ignore.case = T, full.names = T, recursive = T)
    if(length(scrna_files) > 0){
      for(i in 1:nrow(pheno_data)){
      cfile <- NULL
      cfile <- scrna_files[grep(pheno_data[i,"SAMPLE_ID"], scrna_files, ignore.case = T)]
      if((length(cfile) > 0)){
      scrna_data <- Read10X_h5(filename = cfile)
      if((length(scrna_data) == 1) | (length(scrna_data) > 10)){
        scrna_data <- CreateSeuratObject(counts = scrna_data, project = pheno_data[i,"SAMPLE_ID"], min.cells = 3, min.features = 200)
      }else{
        scrna_data <- CreateSeuratObject(counts = scrna_data$`Gene Expression`, project = pheno_data[i,"SAMPLE_ID"], min.cells = 3, min.features = 200)
      }
      scrna_data[["Percent_Mito"]] <- PercentageFeatureSet(scrna_data, pattern = "^MT-")
      scrna_data <- subset(scrna_data, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & Percent_Mito < 5)
    scrna_data <- SCTransform(scrna_data, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE)
    scrna_data <- RunUMAP(scrna_data, reduction = "pca", dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30))
    scrna_data <- RunTSNE(scrna_data, reduction = "pca", dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30), check_duplicates = FALSE)
    scrna_data <- FindNeighbors(scrna_data, reduction = "pca", dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30))
    scrna_data <- FindClusters(scrna_data, verbose = FALSE)

    plotx <- gen10x_plotx(scrna_data, include_meta = T)
    plotx$Cluster <- scrna_data$seurat_clusters

    somePDFPath = paste(cdir,"13URSA_PLOT_scRNASEQ_UMAP_BY_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=10,pointsize=10)
    p <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", "Cluster", paste("scRNASEQ - UMAP: ", pheno_data[i,"SAMPLE_ID"], sep = ""), col = NULL, numeric = TRUE, point_size = 1)
    print(p)
    dev.off()

    somePDFPath = paste(cdir,"14URSA_PLOT_scRNASEQ_tSNE_BY_CLUSTERS_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=10,pointsize=10)
    p <- plot_bygroup(plotx, "tSNE_1", "tSNE_2", "Cluster", paste("scRNASEQ- tSNE: ", pheno_data[i,"SAMPLE_ID"], sep = ""), col = NULL, numeric = TRUE, point_size = 2)
    print(p)
    dev.off()

    Idents(scrna_data) <- "seurat_clusters"
    DefaultAssay(scrna_data) <- "RNA"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(scrna_data)),
                       clusters =  scrna_data$seurat_clusters,
                       ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.fine)

    scrna_data$Cell_Type <- clu_ann$labels[match(scrna_data$seurat_clusters,row.names(clu_ann))]
    scrna_data@meta.data[which(is.na(scrna_data$Cell_Type)),"Cell_Type"] <- "Unidentifiable"
    Idents(scrna_data) <- scrna_data$Cell_Type
    plotx <- gen10x_plotx(scrna_data, include_meta = T)

    p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "Cell_Type", plot_title = paste(pheno_data[i,"SAMPLE_ID"], ": scRNASEQ UMAP - CELL TYPE", sep = ""), col = color_conditions$tenx, annot = TRUE, legend_position = "bottom", point_size = 1)
    p2 <- plot_bygroup(plotx, x = "tSNE_1", y = "tSNE_2", group = "Cell_Type", plot_title = paste(pheno_data[i,"SAMPLE_ID"], ": scRNASEQ tSNE - CELL TYPE", sep = ""),
                       col = color_conditions$tenx, annot = TRUE, legend_position = "bottom", point_size = 1)

    somePDFPath = paste(cdir,"15URSA_PLOT_scRNASEQ_UMAP_TSNE_AUTO_CELL_TYPE_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=10,pointsize=12)
    print((p1|p2)+plot_annotation(title = paste(project_name,sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
    )
    dev.off()

    Idents(scrna_data) <- "Cell_Type"
    DefaultAssay(scrna_data) <- "SCT"
    DefaultAssay(current) <- "SCT"
    Idents(current) <- "Cell_Type"

    data_anchors <- FindTransferAnchors(reference = scrna_data, query = current, normalization.method = "SCT")
    predictions_assay <- TransferData(anchorset = data_anchors, refdata = scrna_data$Cell_Type,
                                      prediction.assay = TRUE,
                                      weight.reduction = current[["pca"]], dims = 1:ifelse(length(scrna_data@reductions$pca) < 30, length(scrna_data@reductions$pca), 30))
    current[["predictions"]] <- predictions_assay

    DefaultAssay(current) <- "predictions"
    Idents(current) <- "Cell_Type"
    cell_types <- unique(scrna_data$Cell_Type)
    p <- SpatialFeaturePlot(current, features = row.names(current)[grep("^max$", row.names(current), ignore.case = T, invert = T)], pt.size.factor = 1.6, ncol = 2, crop = TRUE)
    p <- p+plot_annotation(title = paste("CELL TYPE PREDICTION SCORE: ", pheno_data[i,"SAMPLE_ID"], sep = ""), theme = theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)))

    somePDFPath = paste(cdir,"16URSA_PLOT_scRNASEQ_IMAGE_BY_scRNASEQ_PREDICTION_SCORES_",pheno_data[i,"SAMPLE_ID"],"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=3*length(cell_types),pointsize=10)
    print(p)
    dev.off()

    DefaultAssay(current) <- "predictions"
    p <- VlnPlot(current, group.by = "seurat_clusters", features = row.names(current)[grep("^max$", row.names(current), ignore.case = T, invert = T)],
                 pt.size = 0, ncol = 2, cols = gen_colors(color_conditions$tenx, length(unique(current$seurat_clusters))))
    p <- p+plot_annotation(title = paste(pheno_data[i,"SAMPLE_ID"], ": CELL TYPE PREDICTION SCORES PER CLUSTER", sep = ""), theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

    somePDFPath = paste(cdir,"17URSA_PLOT_scRNASEQ_VIOLIN_PREDICTION_SCORES_CLUSTERS_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=14, height=6*ceiling(length(cell_types)/2),pointsize=10)
    print(p&xlab("CLUSTERS"))
    dev.off()
      }
      }
    }
    }
saveRDS(data_current, paste(cdir,"13URSA_SPATIAL_DATA_",annot_names[i],"_",project_name,".RDS", sep = ""))
print("Completed!")
}
