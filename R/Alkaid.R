############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Alkaid: scATAC
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
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import celldex
#' @import chromVAR
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v86
#' @import GenomeInfoDb
#' @import JASPAR2020
#' @import motifmatchr
#' @import Seurat
#' @import Signac
#' @import SingleR
#' @import TFBSTools
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format scATAC Transcriptomics pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification scATAC
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_scATAC' by default.
#' @param input_dir Directory to all input files. Current working directory by
#' default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format
#' files.
#' @export
#'
scATACPip <- function(project_name = "Ursa_scATAC",
                       input_dir = "./",
                       output_dir = "./",
                       pheno_file){
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "scATAC", isDir = T)
  pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, sep = "_")
  color_conditions <- color_ini()
  ctime <- time_ini()
  sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SID)))
  names(sample_colors) <- unique(pheno_data$SID)
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir)

  pheno_peak <- list.files(input_dir, pattern = ".*peak.*h5", full.names = T, ignore.case = T, recursive = T)
  pheno_singlecell <- list.files(input_dir, pattern = ".*singlecell\\.csv$|.*per_barcode_metrics\\.csv$", full.names = T, ignore.case = T, recursive = T)
  pheno_fragments <- list.files(input_dir, pattern = ".*fragment.*tsv.*", full.names = T, ignore.case = T, recursive = T)
  pheno_fragments <- pheno_fragments[grep("\\.tbi", pheno_fragments, ignore.case = T, invert = T)]
  error_files <- list.files(input_dir, pattern = "gz\\.tsv$", full.names = T, ignore.case = T)

  if(length(error_files) > 0){
    for(i in 1:length(error_files)){
      current_corrected_file <- gsub("gz\\.tsv","gz",error_files[i])
      system(paste("mv ", error_files[i], " ",current_corrected_file, sep = ""))
    }
  }

  data <- NULL
  sample_names <- NULL
  results <- NULL
  for(i in 1:nrow(pheno_data)){
    print(paste("Running ", pheno_data[i,"FILE"], "..", sep = ""))
    sample_names <- c(sample_names, pheno_data[i,"SID"])
    current <- NULL
    current <- pheno_peak[grep(pheno_data[i,"FILE"], pheno_peak, ignore.case = T)]
    current <- Read10X_h5(current)
    current_singlecell <- NULL
    current_singlecell <- pheno_singlecell[grep(pheno_data[i,"FILE"], pheno_singlecell, ignore.case = T)]

    if(length(current_singlecell) > 0){
      current_meta <- "YES"
      current_singlecell <- read.csv(file = current_singlecell, header = TRUE, row.names = 1)
    }else{
      current_meta <- "NO"
    }

    options(download.file.method = "curl")
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      ref_genome <- "hg19"
      ref_annot <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
      seqlevelsStyle(ref_annot) <- "UCSC"
      genome(ref_annot) <- "hg19"
    }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      ref_genome <- "GRCh38.p12"
      ref_annot <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      seqlevelsStyle(ref_annot) <- "UCSC"
      genome(ref_annot) <- "GRCh38.p12"
    }else{
      stop("You are supplying non-human samples. Please try again.")
    }


    print(paste("Constructing chromatin assay of ", pheno_data[i,"FILE"], "..", sep = ""))

    current_assay <- CreateChromatinAssay(
      counts = current,
      sep = c(unique(gsub(".*?([-.:])[0-9]+([-.:])[0-9]+","\\1",row.names(current))),
              unique(gsub(".*?([-.:])[0-9]+([-.:])[0-9]+","\\2",row.names(current)))),
      genome = ref_genome,
      fragments = pheno_fragments[grep(pheno_data[i,"FILE"], pheno_fragments, ignore.case = T)],
      min.cells = 10,
      min.features = 200)

    if(current_meta == "YES"){
      current_seurat <- CreateSeuratObject(counts = current_assay, assay = "peaks",meta.data = current_singlecell)
      rm(current_singlecell)
    }else{
      current_seurat <- CreateSeuratObject(counts = current_assay, assay = "peaks")
    }

    current_seurat <- add_names(current_seurat, pheno_data[i,"SID"], current_ident = "orig.ident")

    rm(current_assay)
    rm(current)


    print(paste("Finding nucelosome signal for ", pheno_data[i,"FILE"], "..", sep = ""))

    DefaultAssay(current_seurat) <- 'peaks'
    Annotation(current_seurat) <- ref_annot
    current_seurat <- NucleosomeSignal(object = current_seurat)

    print(paste("Calculating TSS enrichment for ", pheno_data[i,"FILE"], "..", sep = ""))
    current_seurat <- TSSEnrichment(object = current_seurat, fast = FALSE)

    print(paste("Calculating fragments for ", pheno_data[i,"FILE"], "..", sep = ""))
    if(length(grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat$pct_reads_in_peaks <- current_seurat@meta.data[,grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)] / current_seurat@meta.data[,grep("passed_filters", colnames(current_seurat@meta.data),ignore.case = T)] * 100
      if(length(grep("blacklist_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat$blacklist_ratio <- current_seurat@meta.data[,grep("blacklist_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)] / current_seurat@meta.data[,grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)]
      }
    }
    current_seurat$high.tss <- ifelse(current_seurat@meta.data[,grep("TSS.enrichment",colnames(current_seurat@meta.data), ignore.case = T)] > 2, 'High', 'Low')

    p <- NULL
    p <- TSSPlot(current_seurat, group.by = 'high.tss') + NoLegend()+ scale_color_manual(values = color_conditions$alternate, length(current_seurat$high.tss))

    results[['p1plots']][[i]] <- p
    names(results[['p1plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"1URSA_PLOT_scATAC_NORM_TSS_ENRICHMENT_SCORES_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
    print(results[['p1plots']][[i]]+ggtitle(paste("TSS Enrichment: ", pheno_data[i,"SID"], sep = "")))
    dev.off()

    p <- NULL
    current_col <- grep("nucleosome_signal",colnames(current_seurat@meta.data), ignore.case = T)
    current_seurat$nucleosome_group <- ifelse(current_seurat@meta.data[,current_col] > 4, 'NS > 4', 'NS < 4')
    current_groups <- unique(current_seurat$nucleosome_group)
    if(length(current_groups) == 1){
      p <- FragmentHistogram(object = current_seurat) +
        scale_fill_manual(values = color_conditions$manycolors, length(current_seurat$nucleosome_group))+ggtitle(paste("FRAGMENT LENGTH PERIODICITY: ", current_groups, sep = ""))
    }else{
      p <- FragmentHistogram(object = current_seurat, group.by = 'nucleosome_group') +
        scale_fill_manual(values = color_conditions$manycolors, length(current_seurat$nucleosome_group))+ggtitle("FRAGMENT LENGTH PERIODICITY")
    }

    results[['p2plots']][[i]] <- p
    names(results[['p2plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"2URSA_PLOT_scATAC_FRAGMENT_LENGTH_PERIODICITY_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
    print(results[['p2plots']][[i]]+ggtitle(paste("Fragment Length Periodicity: ", pheno_data[i,"SID"], sep = "")))
    dev.off()

    colnames(current_seurat@meta.data)[grep("pct_reads_in_peaks", colnames(current_seurat@meta.data), ignore.case = T)] <- "PCT_Reads_in_Peaks"
    colnames(current_seurat@meta.data)[grep("blacklist_ratio", colnames(current_seurat@meta.data), ignore.case = T)] <- "Blacklist_Ratio"
    colnames(current_seurat@meta.data)[grep("peak_region_fragments", colnames(current_seurat@meta.data), ignore.case = T)] <- "Peak_Region_Fragments"
    colnames(current_seurat@meta.data)[grep("TSS.enrichment", colnames(current_seurat@meta.data), ignore.case = T)] <- "TSS_Enrichment"
    colnames(current_seurat@meta.data)[grep("nucleosome_signal", colnames(current_seurat@meta.data), ignore.case = T)] <- "Nucleosome_Signal"

    current_features <- colnames(current_seurat@meta.data)[grep("PCT_Reads_in_Peaks|Peak_Region_Fragments$|TSS_Enrichment$|Blacklist_Ratio$|Nucleosome_Signal$", colnames(current_seurat@meta.data), ignore.case = T)]

    p <- NULL

    results[['p3plots']][[i]] <- violin_plot(current_seurat, features = current_features, ncol = length(current_features),
                                             col = sample_colors, x_lab = "SID", log_status = TRUE)
    names(results[['p3plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"3URSA_PLOT_scATAC_QUALITY_CONTROL_PEAKS_LOGSCALE_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 150)
    print(results[['p3plots']][[i]])
    dev.off()

    if(length(grep("Peak_Region_Fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Peak_Region_Fragments > 3000 & Peak_Region_Fragments < 20000)
    }

    if(length(grep("^nCount_peaks$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = nCount_peaks < 100000 & nCount_peaks > 1000)
    }

    if(length(grep("^PCT_Reads_in_Peaks$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = PCT_Reads_in_Peaks > 15)
    }

    if(length(grep("^Blacklist_Ratio$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Blacklist_Ratio < 0.05)
    }

    if(length(grep("^Nucleosome_Signal$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Nucleosome_Signal < 4)
    }

    if(length(grep("^TSS_Enrichment$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = TSS_Enrichment > 2)
    }

    if(length(grep("^nCount_RNA$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = nCount_RNA < 25000 & nCount_RNA > 1000)
    }

    if(length(grep("^nFeature_RNA$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
    }

    if(length(grep("^Percent_Mito$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Percent_Mito < 5)
    }


    print(paste("Running dimension reduction and clustering for ", pheno_data[i,"FILE"], "..", sep = ""))
    DefaultAssay(current_seurat) <- "peaks"
    current_seurat <- RunTFIDF(current_seurat)
    current_seurat <- FindTopFeatures(current_seurat, min.cutoff = 'q0')
    current_seurat <- RunSVD(current_seurat)
    results[['p4plots']][[i]] <- DepthCor(current_seurat)
    names(results[['p4plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"4URSA_PLOT_scATAC_CORRELATION_SEQUENCING_DEPTH_LSI_COMPONENTS_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
    print(results[['p4plots']][[i]])
    dev.off()

    dim_i <- 2 # Skip first component with high correlation
    m <- ifelse((length(current_seurat@assays$peaks@var.features)-1) < 30, (length(current_seurat@assays$peaks@var.features) - 1), 30)
    current_seurat <- RunUMAP(object = current_seurat, reduction = 'lsi', dims = dim_i:m)
    current_seurat <- FindNeighbors(object = current_seurat, reduction = 'lsi', dims = dim_i:m)
    current_seurat <- FindClusters(object = current_seurat, verbose = FALSE, algorithm = 3)

    plotx <- NULL
    plotx <- gen10x_plotx(current_seurat, selected = c("UMAP"), include_meta = T)
    plotx$CLUSTER <- plotx$seurat_clusters
    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(plotx$CLUSTER)))
    names(cluster_colors) <- sort(unique(plotx$CLUSTER), decreasing = F)

    p <- NULL
    p <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "CLUSTER", plot_title = project_name, col = cluster_colors, annot = TRUE,
                      numeric = T, legend_position = "right", point_size = 1, label_size = 8)
    results[['p5plots']][[i]] <- p
    names(results[['p5plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"5SCA_scATAC_UMAP_scATAC_AUTOCLUST_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 3000, units = "px", res = 300)
    print(results[['p5plots']][[i]])
    dev.off()


    print(paste("Finding top peaks for ", pheno_data[i,"FILE"], "..", sep = ""))
    p_thresh <- 0.05
    DefaultAssay(current_seurat) <- "peaks"
    da_peaks <- NULL
    da_peaks <- FindAllMarkers(
      object = current_seurat,
      min.pct = 0.2,
      test.use = 'LR',
      return.thresh = p_thresh,
      latent.vars = if(is.null(current_seurat@meta.data$Peak_Region_Fragments)) newobj <- NULL else newobj <- "Peak_Region_Fragments")

    if(nrow(da_peaks) > 0){
      da_peaks <- da_peaks[da_peaks$p_val_adj < p_thresh,]
      wt <- colnames(da_peaks)[grep("log.*FC", colnames(da_peaks), ignore.case = T)]

      da_peaks <- da_peaks[order(da_peaks[,wt], decreasing = TRUE),]
      results[['p6data']][[i]] <- da_peaks
      names(results[['p6data']])[i] <- pheno_data[i,"SID"]

      write.csv(results[['p6data']][[i]], paste(cdir, "6URSA_TABLE_scATAC_TOP_PEAKS_IN_CLUSTERS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

      n <- 6
      topn <- split(da_peaks, da_peaks$cluster)
      topn <- lapply(topn, function(x){
        x <- x[order(x[,wt]),]
        x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
      })
      topn <- do.call(rbind.data.frame, topn)
      results[['p7data']][[i]] <- topn
      names(results[['p7data']])[i] <- pheno_data[i,"SID"]
      write.csv(results[['p7data']][[i]], paste(cdir, "7URSA_TABLE_scATAC_TOP_",n,"_PEAKS_IN_CLUSTERS_EASY_TABLE_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

      closest_genes <- NULL
      cclusters <- unique(topn$cluster)
      for(k in 1:length(cclusters)){
        current <- unique(topn[which(topn$cluster == cclusters[k]),"gene"])
        current <- data.frame(cluster = cclusters[k], ClosestFeature(current_seurat, regions = current))
        closest_genes <- rbind(closest_genes, current)
      }
      results[['p8data']][[i]] <- closest_genes
      names(results[['p8data']])[i] <- pheno_data[i,"SID"]
      write.csv(results[['p8data']][[i]], paste(cdir, "8URSA_TABLE_scATAC_CLOSEST_GENE_NEAR_TOP_PEAKS_IN_CLUSTERS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)
    }

    print(paste("Done with top peaks search for ", pheno_data[i,"FILE"], "..", sep = ""))

    cref <- NULL
    DefaultAssay(current_seurat) <- "peaks"
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      cref <- BSgenome.Hsapiens.UCSC.hg19
      main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)
    }else if(length(grep("hg38|grch38|38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      cref <- BSgenome.Hsapiens.UCSC.hg38
      main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
    }

    keep.peaks <- which(as.character(seqnames(granges(current_seurat))) %in% main.chroms)
    current_seurat <- subset(current_seurat, features = rownames(current_seurat)[keep.peaks])

    print(paste("Deciphering gene activities for ", pheno_data[i,"FILE"], "..", sep = ""))
    gene_activities <- GeneActivity(current_seurat)
    current_seurat[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
    DefaultAssay(current_seurat) <- "ACTIVITY"
    current_seurat <- NormalizeData(current_seurat)
    current_seurat <- ScaleData(current_seurat, features = rownames(current_seurat))
    Idents(current_seurat) <- "seurat_clusters"

    print(paste("Finding top genes for ", pheno_data[i,"FILE"], "..", sep = ""))
    current_data_markers <- FindAllMarkers(current_seurat, min.pct = 0.25, logfc.threshold = 0.25)
    p_thresh <- 0.05
    current_data_markers <- current_data_markers[current_data_markers$p_val_adj < p_thresh,]
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]
    results[['p9data']][[i]] <- current_data_markers
    names(results[['p9data']])[i] <- pheno_data[i,"SID"]
    write.csv(results[['p9data']][[i]], paste(cdir, "9URSA_TABLE_scATAC_TOP_GENE_ACTIVITY_IN_CLUSTERS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

    DefaultAssay(current_seurat) <- "ACTIVITY"
    wt <- colnames(da_peaks)[grep("log.*FC", colnames(da_peaks), ignore.case = T)]
    top1 <- split(da_peaks, da_peaks$cluster)
    top1 <- lapply(top1, function(x){
      x <- x[order(x[,wt]),]
      if(nrow(x) > 0){
        x <- x[1:ifelse(nrow(x)>1, 1, nrow(x)),]
      }
    })
    top1 <- do.call(rbind.data.frame, top1)

    top1_derivedrna <- split(current_data_markers, current_data_markers$cluster)
    top1_derivedrna <- lapply(top1_derivedrna, function(x){
      x <- x[order(x[,wt]),]
      if(nrow(x) > 0){
        x <- x[1:ifelse(nrow(x)>1, 1, nrow(x)),]
      }
    })

    top1_derivedrna <- do.call(rbind.data.frame, top1_derivedrna)

    print(paste("Generating regional statistics for ", pheno_data[i,"FILE"], "..", sep = ""))
    DefaultAssay(current_seurat) <- "peaks"
    current_seurat <- RegionStats(current_seurat, genome = cref)

    print(paste("Linking peaks to genes for ", pheno_data[i,"FILE"], "..", sep = ""))
    current_seurat <- LinkPeaks(
      object = current_seurat,
      peak.assay = "peaks",
      expression.assay = "ACTIVITY",
      genes.use = unlist(top1_derivedrna$gene))

    clusters <- sort(as.numeric(as.character(unique(top1$cluster))),decreasing = F)
    DefaultAssay(current_seurat) <- "peaks"
    for(k in 1:length(clusters)){
      if(length(which(unlist(top1$cluster) == clusters[k])) > 0){
        results[['p10plots']][[length(results[['p10plots']])+1]] <- CoveragePlot(
          object = current_seurat,
          region = unlist(top1[which(top1$cluster == clusters[k]),"gene"]),
          features = unlist(top1_derivedrna[which(top1_derivedrna$cluster == clusters[k]),"gene"]),
          expression.assay = "ACTIVITY",
          extend.upstream = 10000,
          extend.downstream = 10000)+
          plot_annotation(title = paste("TOP PEAK AND GENE ACTIVITY\n",pheno_data[i,"SID"], ": CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
        names(results[['p10plots']])[length(results[['p10plots']])] <- paste(pheno_data[i,"SID"],"|","CLUSTER_",clusters[k],sep = "")
        somePNGPath <- paste(cdir,"10URSA_PLOT_scATAC_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP1_PEAK_TOP1_GENE_ACTIVITY_IN_CLUSTER_",clusters[k],"_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
        png(somePNGPath, width = 4000, height = 4000*5/6, units = "px", res = 300)
        print(results[['p10plots']][[length(results[['p10plots']])]])
        dev.off()
      }
    }

    n <- 6
    topn_genes <- split(current_data_markers, current_data_markers$cluster)
    topn_genes <- lapply(topn_genes, function(x){
      x <- x[order(x[,wt]),]
      x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
    })
    topn_genes <- do.call(rbind.data.frame, topn_genes)
    results[['p11data']][[i]] <- topn_genes
    names(results[['p11data']])[i] <- pheno_data[i,"SID"]
    write.csv(results[['p11data']][[i]], paste(cdir, "11URSA_TABLE_scATAC_TOP_",n,"_GENES_IN_CLUSTERS_EASY_TABLE_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)


    print(paste("Add motif information for ", pheno_data[i,"FILE"], "..", sep = ""))
    current_assay <- "peaks"
    used_type <- "ACTIVITY"
    DefaultAssay(current_seurat) <- "peaks"
    chosen_genes <- NULL
    position_freq <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- AddMotifs(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = position_freq)
      motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
      if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
        chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
        chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
        current_seurat <- Footprint(object = current_seurat, motif.name = chosen_genes,
                                    genome = BSgenome.Hsapiens.UCSC.hg19)
      }
    }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- AddMotifs(current_seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = position_freq)
      motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
      if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
        chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
        chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
        current_seurat <- Footprint(object = current_seurat, motif.name = chosen_genes,
                                    genome = BSgenome.Hsapiens.UCSC.hg38)
      }
    }

    p <- NULL
    p <- PlotFootprint(current_seurat, features = chosen_genes, group.by	= "seurat_clusters")
    p <- p + plot_layout(ncol = 2)+
      plot_annotation(title = paste(pheno_data[i,"SID"], ": TOP ",used_type," GENES \n(Top group with highest accessibility in motif flanking region are labeled)", sep = ""),
                      theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
    results[['p12plots']][[i]] <- p
    names(results[['p12plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"12URSA_PLOT_scATAC_",toupper(current_assay),"_TOP_",used_type,"_GENES_FOOTPRINTING_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 4000, units = "px", res = 250)
    print(results[['p12plots']][[i]])
    dev.off()

    p_thresh <- 0.05
    DefaultAssay(current_seurat) <- "peaks"

    enriched.motifs <- FindMotifs(object = current_seurat,features = da_peaks[which(da_peaks$p_val_adj < p_thresh),"gene"])
    results[['p13data']][[i]] <- enriched.motifs
    names(results[['p13data']])[i] <- pheno_data[i,"SID"]
    write.csv(results[['p13data']][[i]], paste(cdir, "13URSA_TABLE_scATAC_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

    n <- 10
    results[['p14plots']][[i]] <- MotifPlot(object = current_seurat, motifs = rownames(enriched.motifs)[1:n])+
      plot_annotation(title = paste(pheno_data[i,"SID"], "\nPOSITION WEIGHT MATRICES OF TOP ENRICHED MOTIFS", sep = ""),
                      theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
    results[['p14plots']][[i]] <- adjust_theme(results[['p14plots']][[i]], xsize = 12, title_size = 20, strip_size = 15)
    names(results[['p14plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"14URSA_PLOT_scATAC_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 2800, units = "px", res = 300)
    print(results[['p14plots']][[i]])
    dev.off()


    print(paste("Computing peak set variability for ", pheno_data[i,"FILE"], "..", sep = ""))
    # Per-cell motif activity score
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RunChromVAR(object = current_seurat,genome = BSgenome.Hsapiens.UCSC.hg19)
    }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RunChromVAR(object = current_seurat,genome = BSgenome.Hsapiens.UCSC.hg38)
    }

    n <- 6
    DefaultAssay(current_seurat) <- 'chromvar'
    p <- NULL
    p <- FeaturePlot(
      object = current_seurat,
      features = enriched.motifs$motif[1:n],
      min.cutoff = 'q10',
      max.cutoff = 'q90',
      pt.size = 0.1, label = T,
      cols = c("green", "blue")) +
      plot_annotation(title = "Motif Activities",
                      theme = theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold")))

    results[['p15plots']][[i]] <- results[['p5plots']][[i]] + p +
      plot_annotation(title = paste(pheno_data[i,"SID"], "\nPeaks VS Top Enriched Motifs", sep = ""),
                      theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
    names(results[['p15plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"15URSA_PLOT_scATAC_UMAP_VS_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 2000, units = "px", res = 200)
    print(results[['p15plots']][[i]])
    dev.off()

    DefaultAssay(current_seurat) <- "peaks"
    n <- 1000
    topn <- split(da_peaks, da_peaks$cluster)
    topn <- lapply(topn, function(x){
      x <- x[order(x[,wt]),]
      x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
    })
    topn <- do.call(rbind.data.frame, topn)
    n <- 5000
    topn <- unique(unlist(topn$gene)[1:ifelse(nrow(topn) < n, nrow(topn), n)])
    # current_seurat <- subset(current_seurat, features = topn)
    results[['data']][[i]] <- current_seurat
    names(results[['data']])[i] <- pheno_data[i,"SID"]
    rm(current_seurat)
    rm(closest_genes)
    rm(position_freq)
    # rm(top_FC)
    # rm(top.da.peak)

  }

  data_ref <- unique(toupper(pheno_data[,"REF_GENOME"]))
  if((length(results[['data']]) > 1) & (length(data_ref) == 1)){
    print(paste("Merging samples..", sep = ""))
    combined_peaks <- NULL
    for(i in 1:(length(results[['data']]) - 1)){
      combined_peaks <- reduce(x = c(results[['data']][[i]]@assays$peaks@ranges, results[['data']][[i+1]]@assays$peaks@ranges))
    }

    peakwidths <- width(combined_peaks)
    combined_peaks <- combined_peaks[peakwidths  < 10000 & peakwidths > 20]

    for(i in 1:length(results[['data']])){

      current <- FeatureMatrix(
        fragments = results[['data']][[i]]@assays$peaks@fragments,
        features = combined_peaks,
        cells = colnames(results[['data']][[i]]))

      current <- CreateChromatinAssay(current, fragments = results[['data']][[i]]@assays$peaks@fragments)
      results[['data']][[i]] <- CreateSeuratObject(current, assay = "peaks")
      results[['data']][[i]]$Sample <- sample_names[i]

    }

    rm(current)

    data <- merge(x = results[['data']][[1]], y = results[['data']][2:length(results[['data']])],
                  add.cell.ids = sample_names)

    print(paste("Running dimension reduction and clustering for merged samples..", sep = ""))

    n <- ifelse((length(data@assays$peaks@var.features)-1) < 50, (length(data@assays$peaks@var.features) - 1), 50)
    data <- RunTFIDF(data)
    data <- FindTopFeatures(data, min.cutoff = 20)
    data <- RunSVD(data) # , n = n
    data <- RunUMAP(data, reduction = 'lsi', dims = 2:n)
    # data <- RunTSNE(data, reduction = 'lsi', check_duplicates = FALSE) # dims = 2:n,
    data <- FindNeighbors(data, reduction = 'lsi') # dims = 2:n,
    data <- FindClusters(data, resolution = 0.8, n.iter = 1000)

    print(paste("Finding top peaks for merged samples..", sep = ""))
    data_markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25)

    p_thresh <- 0.05
    data_markers <- data_markers[data_markers$p_val_adj < p_thresh,]
    data_markers <- data_markers[order(data_markers$p_val_adj, decreasing = F),]
    results[['data_markers']] <- data_markers
    write.csv(results[['data_markers']], paste(cdir, "16URSA_TABLE_scATAC_TOP_PEAKS_MERGED_SAMPLES.csv",sep = ""), row.names = F, quote = F)

    plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                        Sample = data$Sample, CLUSTER = data$seurat_clusters)
    plotx$Sample <- factor(plotx$Sample)
    plotx$CLUSTER <- factor(plotx$CLUSTER, levels = sort(unique(as.numeric(as.character(plotx$CLUSTER)))))

    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(plotx$CLUSTER)))
    names(cluster_colors) <- sort(unique(plotx$CLUSTER), decreasing = F)

    p1_umap <- NULL
    p1_umap <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "Sample", plot_title = "All Samples Integrated scATAC UMAP\nBy Samples", col = sample_colors, annot = TRUE, legend_position = "right", point_size = 1, numeric = F)

    p2_umap <- NULL
    p2_umap <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "CLUSTER", plot_title = "All Samples Integrated scATAC UMAP\nBy Clusters", col = cluster_colors, annot = TRUE, legend_position = "right", point_size = 1, numeric = T)

    results[['p16plots']] <- p1_umap+p2_umap

    somePNGPath <- paste(cdir,"16URSA_PLOT_scATAC_UMAP_INTEGRATED_ATAC_SAMPLES_AUTOCLUST_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 1800, units = "px", res = 200)
    print(results[['p16plots']])
    dev.off()

    wt <- colnames(results[['data_markers']])[grep("log.*FC", colnames(results[['data_markers']]), ignore.case = T)]
    top1 <- results[['data_markers']] %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1 <- top1[order(top1$cluster, decreasing = F),]

    plots <- NULL
    for(k in 1:length(top1$gene)){

      plots[[k]] <- FeaturePlot(data, features = top1$gene[k], cols = c("green", "blue"),
                                label = T, pt.size = 0.1, max.cutoff = 'q95')+
        ggtitle(paste("CLUSTER: ", top1$cluster[k], "\n",top1$gene[k],sep = ""))
    }

    results[['p17plots']] <- plots
    somePNGPath <- paste(cdir,"17URSA_PLOT_scATAC_INTEGRATED_ATAC_SAMPLES_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = ceiling(length(top1$gene)/4)*1000, units = "px", res = 300)
    print(do.call("grid.arrange", c(results[['p17plots']], ncol=4)))
    dev.off()

    n <- 10
    wt <- colnames(results[['data_markers']])[grep("log.*FC", colnames(results[['data_markers']]), ignore.case = T)]
    top_n <- results[['data_markers']] %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
    current_clusters <- as.numeric(as.character(sort(unique(results[['data_markers']]$cluster))))
    for(k in 1:length(current_clusters)){
      if(length(which(unlist(top_n$cluster) == current_clusters[k])) > 0){
        current <- top_n[which(top_n$cluster == current_clusters[k]),]
        current <- current[order(current$p_val_adj, decreasing = F),]
        results[['p18plots']][[length(results[['p18plots']])+1]] <- VlnPlot(data, features = current$gene, pt.size = 0,
                                                                            # log = TRUE,
                                                                            cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) +
          plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", current_clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
        names(results[['p18plots']])[length(results[['p18plots']])] <- paste("CLUSTER_",current_clusters[k], sep = "")
        somePNGPath <- paste(cdir,"18URSA_PLOT_scATAC_INTEGRATED_ATAC_SAMPLES_TOP_",n,"_MARKERS_IN_CLUSTER_",current_clusters[k],"_",project_name,".png", sep = "")
        png(somePNGPath, width = 4000, height = ceiling(length(current_clusters)/4)*800, units = "px", res = 300)
        print(results[['p18plots']][[length(results[['p18plots']])]])
        dev.off()
      }
    }

    options(download.file.method = "curl")
    if(length(grep("hg19",data_ref, ignore.case = T)) > 0){
      ref_genome <- "hg19"
      ref_annot <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
      seqlevelsStyle(ref_annot) <- "UCSC"
      genome(ref_annot) <- "hg19"
    }else if(length(grep("hg38|grch38",data_ref, ignore.case = T)) > 0){
      ref_genome <- "GRCh38.p12"
      ref_annot <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      seqlevelsStyle(ref_annot) <- "UCSC"
      genome(ref_annot) <- "GRCh38.p12"
    }

    Annotation(data) <- ref_annot
    rm(ref_annot)

    print(paste("Calculating gene activities for merged samples..", sep = ""))
    gene_activities <- GeneActivity(data)
    data[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
    rm(gene_activities)

    DefaultAssay(data) <- "ACTIVITY"

    print(paste("Finding top genes for merged samples..", sep = ""))
    data_activity_markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25)
    p_thresh <- 0.05
    data_activity_markers <- data_activity_markers[data_activity_markers$p_val_adj < p_thresh,]
    data_activity_markers <- data_activity_markers[order(data_activity_markers$p_val_adj, decreasing = F),]
    results[['data_activity_markers']] <- data_activity_markers
    write.csv(results[['data_activity_markers']], paste(cdir, "18URSA_TABLE_scATAC_TOP_GENES_MERGED_SAMPLES.csv",sep = ""), row.names = F, quote = F)

    wt <- colnames(data_activity_markers)[grep("log.*FC", colnames(data_activity_markers), ignore.case = T)]
    top1_activities <- data_activity_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1_activities <- top1_activities[order(top1_activities$cluster, decreasing = F),]

    top1 <- data_activity_markers %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
    top1 <- top1[order(top1$cluster, decreasing = F),]

    clusters <- unique(sort(as.numeric(as.character(top1$cluster))))

    DefaultAssay(data) <- "peaks"
    for(k in 1:length(clusters)){
      if(length(which(unlist(top1$cluster) == clusters[k])) > 0){
        results[['p19plots']][[length(results[['p19plots']])+1]] <- CoveragePlot(
          object = data,
          expression.assay = "ACTIVITY",
          region = unlist(top1[which(top1$cluster == clusters[k]),"gene"]),
          features = unlist(top1_activities[which(top1_activities$cluster == clusters[k]),"gene"]),
          extend.upstream = 10000,
          extend.downstream = 10000)+
          # scale_color_manual(values = gen_colors(color_conditions$tenx, length(unique(data$Sample))))+
          plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
        names(results[['p19plots']])[length(results[['p19plots']])] <- paste("CLUSTER_",clusters[k], sep = "")
        somePNGPath <- paste(cdir,"19URSA_PLOT_scATAC_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP_PEAK_GENE_ACTIVITY_IN_CLUSTER_",current_clusters[k],"_",project_name,".png", sep = "")
        png(somePNGPath, width = 3000, height = 2000, units = "px", res = 300)
        print(results[['p19plots']][[length(results[['p19plots']])]])
        dev.off()
      }
    }
  }
  saveRDS(results, paste(cdir,"20URSA_DATA_scATAC_",project_name,".RDS", sep = ""))
  print("Completed!")
}

