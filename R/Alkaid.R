############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Alkaid: scATAC
# Version: V1.2.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last Update Date: 2023-01-29
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
#' @param peaks.findallmarkers.multipletestingcorrection Correction method for
#' multiple testing for the 'peaks' assay. Default to "benferroni".
#' Methods include: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". (refer to the p.adjust() function in stats package for details).
#' @param peaks.findallmarkers.padjust_threshold Threshold for retaining entries
#' in the FindAllMarkers differentially accessible peaks for the 'peaks' assay
#' with corrected p-values. Default to 0.1.
#' @param geneactivity.findallmarkers.multipletestingcorrection Correction method for
#' multiple testing for the 'ACTIVITY' assay if gene expression is quantified from the
#' scATAC-seq data. Default to "benferroni". Methods include: "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", "none". (refer to the p.adjust()
#' function in stats package for details).
#' @param geneactivity.findallmarkers.padjust_threshold Threshold for retaining entries
#' in the FindAllMarkers differentially expressed genes for the 'ACTIVITY' assay
#' with corrected p-values. Default to 0.1.
#' @export
#'
scATACPip <- function(project_name = "Ursa_scATAC",
                       input_dir = "./",
                       output_dir = "./",
                       pheno_file,
                      
                      # CreateChromatinAssay
                      signac.createchromatinassay.min.cells = 10,
                      signac.createchromatinassay.min.features = 200,
                      signac.createchromatinassay.max.cells = NULL,
                      signac.createchromatinassay.ranges = NULL,
                      signac.createchromatinassay.motifs = NULL,
                      signac.createchromatinassay.annotation = NULL,
                      signac.createchromatinassay.bias = NULL,
                      signac.createchromatinassay.positionEnrichment = NULL,
                      signac.createchromatinassay.validate.fragments = TRUE,
                      
                      # TSSEnrichment
                      signac.tssenrichment.tss.positions = NULL,
                      signac.tssenrichment.n = NULL,
                      signac.tssenrichment.fast = FALSE,
                      signac.tssenrichment.assay = NULL,
                      signac.tssenrichment.cells = NULL,
                      signac.tssenrichment.process_n = 2000,
                      
                      # subset
                      # Peak_Region_Fragments
                      signac.min.peak_region_fragments = 3000,
                      signac.max.peak_region_fragments = 20000,
                      
                      # nCount_peaks
                      signac.min.ncount_peaks = 1000,
                      signac.max.ncount_peaks = 100000,
                      
                      # PCT_Reads_in_Peaks
                      signac.min.pct_reads_in_peaks = 15,
                      
                      # Blacklist_Ratio
                      signac.max.blacklist_ratio = 0.05,
                      
                      # Nucleosome_Signal
                      signac.max.nucleosome_signal = 4,
                      
                      # TSS_Enrichment
                      signac.min.tss_enrichment = 2,
                      
                      # nCount_RNA (if scRNA-seq assay is present)
                      signac.min.ncount_rna = 1000,
                      signac.max.ncount_rna = 25000,
                      
                      # nFeature_RNA (if scRNA-seq assay is present)
                      signac.min.nfeature_rna = 200,
                      signac.max.nfeature_rna = 2500,
                      
                      # Percent_Mito (if scRNA-seq assay is present)
                      signac.max.percent_mito = 10,
                      
                      # RunTFIDF
                      signac.runtfidf.assay = NULL,
                      signac.runtfidf.method = 1,
                      signac.runtfidf.scale.factor = 10000,
                      signac.runtfidf.idf = NULL,
                      
                      # FindTopFeatures
                      signac.findtopfeatures.assay = NULL,
                      signac.findtopfeatures.min.cutoff = "q0",
                      
                      # RunSVD
                      runsvd.assay = NULL,
                      runsvd.features = NULL,
                      runsvd.n = 50,
                      runsvd.reduction.key = "LSI_",
                      runsvd.reduction.name = "lsi",
                      runsvd.scale.max = NULL,
                      
                      # DepthCor
                      signac.depthcor.assay = NULL,
                      signac.depthcor.reduction = "lsi",
                      signac.depthcor.n = 10,
                      
                      # RunUMAP
                      signac.runumap.lsi.dimension = 30,
                      seurat.runumap.features = NULL,
                      seurat.runumap.graph = NULL,
                      seurat.runumap.assay = NULL,
                      seurat.runumap.nn.name = NULL,
                      seurat.runumap.slot = "data",
                      seurat.runumap.umap.method = "uwot",
                      seurat.runumap.reduction.model = NULL,
                      seurat.runumap.return.model = FALSE,
                      seurat.runumap.n.neighbors = 30L,
                      seurat.runumap.n.components = 2L,
                      seurat.runumap.metric = "cosine",
                      seurat.runumap.n.epochs = NULL,
                      seurat.runumap.learning.rate = 1,
                      seurat.runumap.min.dist = 0.3,
                      seurat.runumap.spread = 1,
                      seurat.runumap.set.op.mix.ratio = 1,
                      seurat.runumap.local.connectivity = 1L,
                      seurat.runumap.repulsion.strength = 1,
                      seurat.runumap.negative.sample.rate = 5L,
                      seurat.runumap.a = NULL,
                      seurat.runumap.b = NULL,
                      seurat.runumap.uwot.sgd = FALSE,
                      seurat.runumap.seed.use = 42L,
                      seurat.runumap.metric.kwds = NULL,
                      seurat.runumap.angular.rp.forest = FALSE,
                      seurat.runumap.densmap = FALSE,
                      seurat.runumap.dens.lambda = 2,
                      seurat.runumap.dens.frac = 0.3,
                      seurat.runumap.dens.var.shift = 0.1,
                      seurat.runumap.verbose = TRUE,
                      seurat.runumap.reduction.name = "umap",
                      seurat.runumap.reduction.key = "UMAP_",
                      
                      # FindNeighbors
                      seurat.findneighbors.assay = NULL,
                      seurat.findneighbors.features = NULL,
                      seurat.findneighbors.k.param = 20,
                      seurat.findneighbors.return.neighbor = FALSE,
                      seurat.findneighbors.prune.SNN = 1/15,
                      seurat.findneighbors.nn.method = "annoy",
                      seurat.findneighbors.n.trees = 50,
                      seurat.findneighbors.annoy.metric = "euclidean",
                      seurat.findneighbors.nn.eps = 0,
                      seurat.findneighbors.verbose = TRUE,
                      seurat.findneighbors.force.recalc = FALSE,
                      seurat.findneighbors.do.plot = FALSE,
                      seurat.findneighbors.graph.name = NULL,
                      seurat.findneighbors.l2.norm = FALSE,
                      seurat.findneighbors.cache.index = FALSE,
                      
                      # FindClusters
                      seurat.findclusters.graph.name = NULL,
                      seurat.findclusters.modularity.fxn = 1,
                      seurat.findclusters.initial.membership = NULL,
                      seurat.findclusters.node.sizes = NULL,
                      seurat.findclusters.resolution = 0.8,
                      seurat.findclusters.method = "matrix",
                      seurat.findclusters.algorithm = 3,
                      seurat.findclusters.n.start = 10,
                      seurat.findclusters.n.iter = 1000,
                      seurat.findclusters.random.seed = 0,
                      seurat.findclusters.group.singletons = TRUE,
                      seurat.findclusters.temp.file.location = NULL,
                      seurat.findclusters.edge.file.name = NULL,
                      seurat.findclusters.verbose = TRUE,
                      
                      # Peaks - FindAllMarkers
                      seurat.peaks.findallmarkers.assay = NULL,
                      seurat.peaks.findallmarkers.features = NULL,
                      seurat.peaks.findallmarkers.logfc.threshold = 0.25,
                      seurat.peaks.findallmarkers.test.use = 'LR',
                      seurat.peaks.findallmarkers.slot = "data",
                      seurat.peaks.findallmarkers.min.pct = 0.2,
                      seurat.peaks.findallmarkers.min.diff.pct = -Inf,
                      seurat.peaks.findallmarkers.node = NULL,
                      seurat.peaks.findallmarkers.verbose = TRUE,
                      seurat.peaks.findallmarkers.only.pos = FALSE,
                      seurat.peaks.findallmarkers.max.cells.per.ident = Inf,
                      seurat.peaks.findallmarkers.random.seed = 1,
                      seurat.peaks.findallmarkers.latent.vars = NULL,
                      seurat.peaks.findallmarkers.min.cells.feature = 3,
                      seurat.peaks.findallmarkers.min.cells.group = 3,
                      seurat.peaks.findallmarkers.mean.fxn = NULL,
                      seurat.peaks.findallmarkers.fc.name = NULL,
                      seurat.peaks.findallmarkers.base = 2,
                      seurat.peaks.findallmarkers.return.thresh = 0.05,
                      seurat.peaks.findallmarkers.densify = FALSE,
                      
                      peaks.findallmarkers.padjust_threshold = 0.1,
                      peaks.findallmarkers.multipletestingcorrection = "bonferroni",
                      
                      # GeneActivity
                      seurat.geneactivity.assay = NULL,
                      seurat.geneactivity.features = NULL,
                      seurat.geneactivity.extend.upstream = 2000,
                      seurat.geneactivity.extend.downstream = 0,
                      seurat.geneactivity.biotypes = "protein_coding",
                      seurat.geneactivity.max.width = 5e+05,
                      seurat.geneactivity.process_n = 2000,
                      seurat.geneactivity.gene.id = FALSE,
                      
                      # CreateAssayObject
                      seurat.createassayobject.min.cells = 0,
                      seurat.createassayobject.min.features = 0,
                      seurat.createassayobject.check.matrix = FALSE,
                      
                      # NormalizeData
                      seurat.normalizedata.assay = NULL,
                      seurat.normalizedata.normalization.method = "LogNormalize",
                      seurat.normalizedata.scale.factor = 10000,
                      seurat.normalizedata.margin = 1,
                      
                      # ScaleData
                      seurat.scaledata.features = NULL,
                      seurat.scaledata.vars.to.regress = NULL,
                      seurat.scaledata.split.by = NULL,
                      seurat.scaledata.model.use = "linear",
                      seurat.scaledata.use.umi = FALSE,
                      seurat.scaledata.do.scale = TRUE,
                      seurat.scaledata.do.center = TRUE,
                      seurat.scaledata.scale.max = 10,
                      seurat.scaledata.block.size = 1000,
                      seurat.scaledata.min.cells.to.block = 3000,
                      
                      # Gene Activity - FindAllMarkers
                      seurat.geneactivity.findallmarkers.assay = NULL,
                      seurat.geneactivity.findallmarkers.features = NULL,
                      seurat.geneactivity.findallmarkers.logfc.threshold = 0.25,
                      seurat.geneactivity.findallmarkers.test.use = 'wilcox',
                      seurat.geneactivity.findallmarkers.slot = "data",
                      seurat.geneactivity.findallmarkers.min.pct = 0.25,
                      seurat.geneactivity.findallmarkers.min.diff.pct = -Inf,
                      seurat.geneactivity.findallmarkers.node = NULL,
                      seurat.geneactivity.findallmarkers.verbose = TRUE,
                      seurat.geneactivity.findallmarkers.only.pos = FALSE,
                      seurat.geneactivity.findallmarkers.max.cells.per.ident = Inf,
                      seurat.geneactivity.findallmarkers.random.seed = 1,
                      seurat.geneactivity.findallmarkers.latent.vars = NULL,
                      seurat.geneactivity.findallmarkers.min.cells.feature = 3,
                      seurat.geneactivity.findallmarkers.min.cells.group = 3,
                      seurat.geneactivity.findallmarkers.mean.fxn = NULL,
                      seurat.geneactivity.findallmarkers.fc.name = NULL,
                      seurat.geneactivity.findallmarkers.base = 2,
                      seurat.geneactivity.findallmarkers.return.thresh = 0.05,
                      seurat.geneactivity.findallmarkers.densify = FALSE,
                      
                      geneactivity.findallmarkers.padjust_threshold = 0.1,
                      geneactivity.findallmarkers.multipletestingcorrection = "bonferroni",
                      
                      # RegionStats
                      signac.regionstats.assay = NULL,
                      
                      # LinkPeaks
                      signac.linkpeaks.peak.slot = "counts",
                      signac.linkpeaks.expression.slot = "data",
                      signac.linkpeaks.method = "pearson",
                      signac.linkpeaks.gene.coords = NULL,
                      signac.linkpeaks.distance = 5e+05,
                      signac.linkpeaks.min.distance = NULL,
                      signac.linkpeaks.min.cells = 10,
                      signac.linkpeaks.genes.use = NULL,
                      signac.linkpeaks.n_sample = 200,
                      signac.linkpeaks.pvalue_cutoff = 0.05,
                      signac.linkpeaks.score_cutoff = 0.05,
                      signac.linkpeaks.gene.id = FALSE,
                      
                      # AddMotifs
                      signac.addmotifs.assay = NULL,
                      
                      # Footprint
                      signac.footprint.regions = NULL,
                      signac.footprint.assay = NULL,
                      
                      # FindMotifs
                      signac.findmotifs.background = 40000,
                      signac.findmotifs.assay = NULL,
                      signac.findmotifs.p.adjust.method = "BH",
                      
                      # RunChromVAR
                      signac.runchromvar.motif.matrix = NULL,
                      signac.runchromvar.assay = NULL,
                      
                      # combined_peaks subset
                      signac.post_combined_peaks.max.peaks = 10000,
                      signac.post_combined_peaks.min.peaks = 20){
  
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
  dir.create(cdir, recursive = T)

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
      min.cells = signac.createchromatinassay.min.cells,
      min.features = signac.createchromatinassay.min.features,
      max.cells = signac.createchromatinassay.max.cells,
      ranges = signac.createchromatinassay.ranges,
      motifs = signac.createchromatinassay.motifs,
      annotation = signac.createchromatinassay.annotation,
      bias = signac.createchromatinassay.bias,
      positionEnrichment = signac.createchromatinassay.positionEnrichment,
      validate.fragments = signac.createchromatinassay.validate.fragments)

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
    current_seurat <- TSSEnrichment(object = current_seurat,
                                    tss.positions = signac.tssenrichment.tss.positions,
                                    n = signac.tssenrichment.n,
                                    fast = signac.tssenrichment.fast,
                                    assay = signac.tssenrichment.assay,
                                    cells = signac.tssenrichment.cells,
                                    process_n = signac.tssenrichment.process_n)

    print(paste("Calculating fragments for ", pheno_data[i,"FILE"], "..", sep = ""))
    if(length(grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat$pct_reads_in_peaks <- current_seurat@meta.data[,grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)] / current_seurat@meta.data[,grep("passed_filters", colnames(current_seurat@meta.data),ignore.case = T)] * 100
      if(length(grep("blacklist_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
        current_seurat$blacklist_ratio <- current_seurat@meta.data[,grep("blacklist_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)] / current_seurat@meta.data[,grep("peak_region_fragments", colnames(current_seurat@meta.data),ignore.case = T)]
      }
    }
    current_seurat$high.tss <- ifelse(current_seurat@meta.data[,grep("TSS.enrichment",colnames(current_seurat@meta.data), ignore.case = T)] > 3, 'High', 'Low')

    p <- NULL
    p <- DensityScatter(current_seurat, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
    
    # results[['p0plots']][[i]] <- p
    # names(results[['p0plots']])[i] <- pheno_data[i,"SID"]
    
    somePNGPath <- paste(cdir,"0URSA_PLOT_scATAC_QUALITY_CONTROL_DENSITY_PLOT_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
    print(p+ggtitle(paste(pheno_data[i,"SID"], sep = "")))
    dev.off()

    if(signac.tssenrichment.fast == TRUE){
      p <- NULL
      p <- TSSPlot(current_seurat, group.by = 'high.tss') + NoLegend()+ scale_color_manual(values = color_conditions$alternate, length(current_seurat$high.tss))
      
      # results[['p1plots']][[i]] <- p
      # names(results[['p1plots']])[i] <- pheno_data[i,"SID"]
      
      somePNGPath <- paste(cdir,"1URSA_PLOT_scATAC_NORM_TSS_ENRICHMENT_SCORES_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
      png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
      print(p+ggtitle(paste("TSS Enrichment: ", pheno_data[i,"SID"], sep = "")))
      dev.off()
    }

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

    # results[['p2plots']][[i]] <- p
    # names(results[['p2plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"2URSA_PLOT_scATAC_FRAGMENT_LENGTH_PERIODICITY_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
    print(p+ggtitle(paste("Fragment Length Periodicity: ", pheno_data[i,"SID"], sep = "")))
    dev.off()

    colnames(current_seurat@meta.data)[grep("pct_reads_in_peaks", colnames(current_seurat@meta.data), ignore.case = T)] <- "PCT_Reads_in_Peaks"
    colnames(current_seurat@meta.data)[grep("blacklist_ratio", colnames(current_seurat@meta.data), ignore.case = T)] <- "Blacklist_Ratio"
    colnames(current_seurat@meta.data)[grep("peak_region_fragments", colnames(current_seurat@meta.data), ignore.case = T)] <- "Peak_Region_Fragments"
    colnames(current_seurat@meta.data)[grep("TSS.enrichment", colnames(current_seurat@meta.data), ignore.case = T)] <- "TSS_Enrichment"
    colnames(current_seurat@meta.data)[grep("nucleosome_signal", colnames(current_seurat@meta.data), ignore.case = T)] <- "Nucleosome_Signal"

    current_features <- colnames(current_seurat@meta.data)[grep("PCT_Reads_in_Peaks|Peak_Region_Fragments$|TSS_Enrichment$|Blacklist_Ratio$|Nucleosome_Signal$", colnames(current_seurat@meta.data), ignore.case = T)]

    p <- NULL

    p <- violin_plot(current_seurat, features = current_features, ncol = length(current_features),
                                             col = sample_colors, x_lab = "SID", log_status = TRUE)
    # results[['p3plots']][[i]] <- p
    # names(results[['p3plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"3URSA_PLOT_scATAC_QUALITY_CONTROL_PEAKS_LOGSCALE_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 150)
    print(p)
    dev.off()

    if(length(grep("Peak_Region_Fragments", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Peak_Region_Fragments > signac.min.peak_region_fragments & Peak_Region_Fragments < signac.max.peak_region_fragments)
    }

    if(length(grep("^nCount_peaks$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = nCount_peaks > signac.min.ncount_peaks & nCount_peaks < signac.max.ncount_peaks)
    }

    if(length(grep("^PCT_Reads_in_Peaks$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = PCT_Reads_in_Peaks > signac.min.pct_reads_in_peaks)
    }

    if(length(grep("^Blacklist_Ratio$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Blacklist_Ratio < signac.max.blacklist_ratio)
    }

    if(length(grep("^Nucleosome_Signal$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Nucleosome_Signal < signac.max.nucleosome_signal)
    }

    if(length(grep("^TSS_Enrichment$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = TSS_Enrichment > signac.min.tss_enrichment)
    }

    if(length(grep("^nCount_RNA$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = nCount_RNA > signac.min.ncount_rna & nCount_RNA < signac.max.ncount_rna)
    }

    if(length(grep("^nFeature_RNA$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = nFeature_RNA > signac.min.nfeature_rna & nFeature_RNA < signac.max.nfeature_rna)
    }

    if(length(grep("^Percent_Mito$", colnames(current_seurat@meta.data),ignore.case = T)) > 0){
      current_seurat <- subset(x = current_seurat, subset = Percent_Mito < signac.max.percent_mito)
    }

    print(paste("Running dimension reduction and clustering for ", pheno_data[i,"FILE"], "..", sep = ""))
    DefaultAssay(current_seurat) <- "peaks"
    current_seurat <- RunTFIDF(current_seurat,
                               assay = signac.runtfidf.assay,
                               method = signac.runtfidf.method,
                               scale.factor = signac.runtfidf.scale.factor,
                               idf = signac.runtfidf.idf)
    current_seurat <- FindTopFeatures(current_seurat,
                                      assay = signac.findtopfeatures.assay,
                                      min.cutoff = signac.findtopfeatures.min.cutoff)
    current_seurat <- RunSVD(current_seurat,
                             assay = runsvd.assay,
                             features = runsvd.features,
                             n = runsvd.n,
                             reduction.key = runsvd.reduction.key,
                             reduction.name = runsvd.reduction.name,
                             scale.max = runsvd.scale.max)
    p <- NULL
    p <- DepthCor(current_seurat,
                                          assay = signac.depthcor.assay,
                                          reduction = signac.depthcor.reduction,
                                          n = signac.depthcor.n)
    # results[['p4plots']][[i]]
    # names(results[['p4plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"4URSA_PLOT_scATAC_CORRELATION_SEQUENCING_DEPTH_LSI_COMPONENTS_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 2000, height = 1000, units = "px", res = 300)
    print(p)
    dev.off()

    dim_i <- 2 # Skip first component with high correlation
    m <- ifelse((length(current_seurat@assays$peaks@var.features)-1) < 30, (length(current_seurat@assays$peaks@var.features) - 1), 30)
    if(is.null(signac.runumap.lsi.dimension)){
      signac.runumap.lsi.dimension <- m
    }
    
    current_seurat <- RunUMAP(object = current_seurat,
                              reduction = 'lsi',
                              dims = dim_i:signac.runumap.lsi.dimension,
                              features = seurat.runumap.features,
                              graph = seurat.runumap.graph,
                              assay = seurat.runumap.assay,
    nn.name = seurat.runumap.nn.name,
    slot = seurat.runumap.slot,
    umap.method = seurat.runumap.umap.method,
    reduction.model = seurat.runumap.reduction.model,
    return.model = seurat.runumap.return.model,
    n.neighbors = seurat.runumap.n.neighbors,
    n.components = seurat.runumap.n.components,
    metric = seurat.runumap.metric,
    n.epochs = seurat.runumap.n.epochs,
    learning.rate = seurat.runumap.learning.rate,
    min.dist = seurat.runumap.min.dist,
    spread = seurat.runumap.spread,
    set.op.mix.ratio = seurat.runumap.set.op.mix.ratio,
    local.connectivity = seurat.runumap.local.connectivity,
    repulsion.strength = seurat.runumap.repulsion.strength,
    negative.sample.rate = seurat.runumap.negative.sample.rate,
    a = seurat.runumap.a,
    b = seurat.runumap.b,
    uwot.sgd = seurat.runumap.uwot.sgd,
    seed.use = seurat.runumap.seed.use,
    metric.kwds = seurat.runumap.metric.kwds,
    angular.rp.forest = seurat.runumap.angular.rp.forest,
    densmap = seurat.runumap.densmap,
    dens.lambda = seurat.runumap.dens.lambda,
    dens.frac = seurat.runumap.dens.frac,
    dens.var.shift = seurat.runumap.dens.var.shift,
    verbose = seurat.runumap.verbose,
    reduction.name = seurat.runumap.reduction.name,
    reduction.key = seurat.runumap.reduction.key)

    current_seurat <- FindNeighbors(object = current_seurat,
                                    reduction = 'lsi',
                                    dims = dim_i:signac.runumap.lsi.dimension,
                                    assay = seurat.findneighbors.assay,
                                    features = seurat.findneighbors.features,
                                    k.param = seurat.findneighbors.k.param,
                                    return.neighbor = seurat.findneighbors.return.neighbor,
                                    prune.SNN = seurat.findneighbors.prune.SNN,
                                    nn.method = seurat.findneighbors.nn.method,
                                    n.trees = seurat.findneighbors.n.trees,
                                    annoy.metric = seurat.findneighbors.annoy.metric,
                                    nn.eps = seurat.findneighbors.nn.eps,
                                    verbose = seurat.findneighbors.verbose,
                                    force.recalc = seurat.findneighbors.force.recalc,
                                    do.plot = seurat.findneighbors.do.plot,
                                    graph.name = seurat.findneighbors.graph.name,
                                    l2.norm = seurat.findneighbors.l2.norm,
                                    cache.index = seurat.findneighbors.cache.index)
    
    current_seurat <- FindClusters(object = current_seurat,
                                   graph.name = seurat.findclusters.graph.name,
                                   modularity.fxn = seurat.findclusters.modularity.fxn,
                                   initial.membership = seurat.findclusters.initial.membership,
                                   node.sizes = seurat.findclusters.node.sizes,
                                   resolution = seurat.findclusters.resolution,
                                   method = seurat.findclusters.method,
                                   algorithm = seurat.findclusters.algorithm,
                                   n.start = seurat.findclusters.n.start,
                                   n.iter = seurat.findclusters.n.iter,
                                   random.seed = seurat.findclusters.random.seed,
                                   group.singletons = seurat.findclusters.group.singletons,
                                   temp.file.location = seurat.findclusters.temp.file.location,
                                   edge.file.name = seurat.findclusters.edge.file.name,
                                   verbose = seurat.findclusters.verbose)

    plotx <- NULL
    plotx <- gen10x_plotx(current_seurat, selected = c("UMAP"), include_meta = T)
    plotx$CLUSTER <- plotx$seurat_clusters
    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(plotx$CLUSTER)))
    names(cluster_colors) <- sort(unique(plotx$CLUSTER), decreasing = F)

    p <- NULL
    p <- plot_bygroup(plotx, "UMAP_1", "UMAP_2", group = "CLUSTER", plot_title = project_name, col = cluster_colors, annot = TRUE,
                      numeric = T, legend_position = "right", point_size = 1, label_size = 8)
    # results[['p5plots']][[i]] <- p
    # names(results[['p5plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"5SCA_scATAC_UMAP_scATAC_AUTOCLUST_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 3000, units = "px", res = 300)
    print(p)
    dev.off()

    if(is.null(seurat.peaks.findallmarkers.latent.vars)){
      seurat.peaks.findallmarkers.latent.vars <- if(is.null(current_seurat@meta.data$Peak_Region_Fragments)) newobj <- NULL else newobj <- "Peak_Region_Fragments"
    }
    
    print(paste("Finding top peaks for ", pheno_data[i,"FILE"], "..", sep = ""))
    p_thresh <- 0.05
    DefaultAssay(current_seurat) <- "peaks"
    da_peaks <- NULL
    da_peaks <- FindAllMarkers(
      object = current_seurat,
      assay = seurat.peaks.findallmarkers.assay,
      features = seurat.peaks.findallmarkers.features,
      logfc.threshold = seurat.peaks.findallmarkers.logfc.threshold,
      test.use = seurat.peaks.findallmarkers.test.use,
      slot = seurat.peaks.findallmarkers.slot,
      min.pct = seurat.peaks.findallmarkers.min.pct,
      min.diff.pct = seurat.peaks.findallmarkers.min.diff.pct,
      node = seurat.peaks.findallmarkers.node,
      verbose = seurat.peaks.findallmarkers.verbose,
      only.pos = seurat.peaks.findallmarkers.only.pos,
      max.cells.per.ident = seurat.peaks.findallmarkers.max.cells.per.ident,
      random.seed = seurat.peaks.findallmarkers.random.seed,
      latent.vars = seurat.peaks.findallmarkers.latent.vars,
      min.cells.feature = seurat.peaks.findallmarkers.min.cells.feature,
      min.cells.group = seurat.peaks.findallmarkers.min.cells.group,
      mean.fxn = seurat.peaks.findallmarkers.mean.fxn,
      fc.name = seurat.peaks.findallmarkers.fc.name,
      base = seurat.peaks.findallmarkers.base,
      return.thresh = seurat.peaks.findallmarkers.return.thresh,
      densify = seurat.peaks.findallmarkers.densify)

    if(toupper(peaks.findallmarkers.multipletestingcorrection) != toupper("bonferroni")){
      da_peaks$p_val_adj <- p.adjust(da_peaks$p_val, method = peaks.findallmarkers.multipletestingcorrection)
    }
    
    if(nrow(da_peaks) > 0){
      da_peaks <- da_peaks[da_peaks$p_val_adj < peaks.findallmarkers.padjust_threshold,]
      wt <- colnames(da_peaks)[grep("log.*FC", colnames(da_peaks), ignore.case = T)]

      da_peaks <- da_peaks[order(da_peaks[,wt], decreasing = TRUE),]
      # results[['p6data']][[i]] <- da_peaks
      # names(results[['p6data']])[i] <- pheno_data[i,"SID"]

      write.csv(da_peaks, paste(cdir, "6URSA_TABLE_scATAC_TOP_PEAKS_IN_CLUSTERS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

      n <- 6
      topn <- split(da_peaks, da_peaks$cluster)
      topn <- lapply(topn, function(x){
        x <- x[order(x[,wt]),]
        x <- x[1:ifelse(nrow(x)>n, n, nrow(x)),]
      })
      topn <- do.call(rbind.data.frame, topn)
      # results[['p7data']][[i]] <- topn
      # names(results[['p7data']])[i] <- pheno_data[i,"SID"]
      write.csv(topn, paste(cdir, "7URSA_TABLE_scATAC_TOP_",n,"_PEAKS_IN_CLUSTERS_EASY_TABLE_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

      closest_genes <- NULL
      cclusters <- unique(topn$cluster)
      for(k in 1:length(cclusters)){
        current <- unique(topn[which(topn$cluster == cclusters[k]),"gene"])
        current <- data.frame(cluster = cclusters[k], ClosestFeature(current_seurat, regions = current))
        closest_genes <- rbind(closest_genes, current)
      }
      # results[['p8data']][[i]] <- closest_genes
      # names(results[['p8data']])[i] <- pheno_data[i,"SID"]
      write.csv(closest_genes, paste(cdir, "8URSA_TABLE_scATAC_CLOSEST_GENE_NEAR_TOP_PEAKS_IN_CLUSTERS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)
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
    gene_activities <- GeneActivity(current_seurat,
                                    assay = seurat.geneactivity.assay,
                                    features = seurat.geneactivity.features,
                                    extend.upstream = seurat.geneactivity.extend.upstream,
                                    extend.downstream = seurat.geneactivity.extend.downstream,
                                    biotypes = seurat.geneactivity.biotypes,
                                    max.width = seurat.geneactivity.max.width,
                                    process_n = seurat.geneactivity.process_n,
                                    gene.id = seurat.geneactivity.gene.id)
    
    current_seurat[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities,
                                                      min.cells = seurat.createassayobject.min.cells,
                                                      min.features = seurat.createassayobject.min.features,
                                                      check.matrix = seurat.createassayobject.check.matrix)
    DefaultAssay(current_seurat) <- "ACTIVITY"
    current_seurat <- NormalizeData(current_seurat,
                                    assay = seurat.normalizedata.assay,
                                    normalization.method = seurat.normalizedata.normalization.method,
                                    scale.factor = seurat.normalizedata.scale.factor,
                                    margin = seurat.normalizedata.margin)
    
    current_seurat <- ScaleData(current_seurat,
                                features = seurat.scaledata.features,
                                vars.to.regress = seurat.scaledata.vars.to.regress,
                                split.by = seurat.scaledata.split.by,
                                model.use = seurat.scaledata.model.use,
                                use.umi = seurat.scaledata.use.umi,
                                do.scale = seurat.scaledata.do.scale,
                                do.center = seurat.scaledata.do.center,
                                scale.max = seurat.scaledata.scale.max,
                                block.size = seurat.scaledata.block.size,
                                min.cells.to.block = seurat.scaledata.min.cells.to.block)
    
    Idents(current_seurat) <- "seurat_clusters"

    print(paste("Finding top genes for ", pheno_data[i,"FILE"], "..", sep = ""))
    current_data_markers <- FindAllMarkers(current_seurat,
                                           assay = seurat.geneactivity.findallmarkers.assay,
                                           features = seurat.geneactivity.findallmarkers.features,
                                           logfc.threshold = seurat.geneactivity.findallmarkers.logfc.threshold,
                                           test.use = seurat.geneactivity.findallmarkers.test.use,
                                           slot = seurat.geneactivity.findallmarkers.slot,
                                           min.pct = seurat.geneactivity.findallmarkers.min.pct,
                                           min.diff.pct = seurat.geneactivity.findallmarkers.min.diff.pct,
                                           node = seurat.geneactivity.findallmarkers.node,
                                           verbose = seurat.geneactivity.findallmarkers.verbose,
                                           only.pos = seurat.geneactivity.findallmarkers.only.pos,
                                           max.cells.per.ident = seurat.geneactivity.findallmarkers.max.cells.per.ident,
                                           random.seed = seurat.geneactivity.findallmarkers.random.seed,
                                           latent.vars = seurat.geneactivity.findallmarkers.latent.vars,
                                           min.cells.feature = seurat.geneactivity.findallmarkers.min.cells.feature,
                                           min.cells.group = seurat.geneactivity.findallmarkers.min.cells.group,
                                           mean.fxn = seurat.geneactivity.findallmarkers.mean.fxn,
                                           fc.name = seurat.geneactivity.findallmarkers.fc.name,
                                           base = seurat.geneactivity.findallmarkers.base,
                                           return.thresh = seurat.geneactivity.findallmarkers.return.thresh,
                                           densify = seurat.geneactivity.findallmarkers.densify)
    
    if(toupper(geneactivity.findallmarkers.multipletestingcorrection) != toupper("bonferroni")){
      current_data_markers$p_val_adj <- p.adjust(current_data_markers$p_val, method = geneactivity.findallmarkers.multipletestingcorrection)
    }

    current_data_markers <- current_data_markers[current_data_markers$p_val_adj < geneactivity.findallmarkers.padjust_threshold,]
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]
    # results[['p9data']][[i]] <- current_data_markers
    # names(results[['p9data']])[i] <- pheno_data[i,"SID"]
    write.csv(current_data_markers, paste(cdir, "9URSA_TABLE_scATAC_TOP_GENE_ACTIVITY_IN_CLUSTERS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

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
    current_seurat <- RegionStats(current_seurat,
                                  genome = cref,
                                  assay = signac.regionstats.assay)

    print(paste("Linking peaks to genes for ", pheno_data[i,"FILE"], "..", sep = ""))
    
    if(is.null(signac.linkpeaks.genes.use)){
      signac.linkpeaks.genes.use <- unlist(top1_derivedrna$gene)
    }
      
      
    current_seurat <- LinkPeaks(
      object = current_seurat,
      peak.assay = "peaks",
      expression.assay = "ACTIVITY",
      peak.slot = signac.linkpeaks.peak.slot,
      expression.slot = signac.linkpeaks.expression.slot,
      method = signac.linkpeaks.method,
      gene.coords = signac.linkpeaks.gene.coords,
      distance = signac.linkpeaks.distance,
      min.distance = signac.linkpeaks.min.distance,
      min.cells = signac.linkpeaks.min.cells,
      genes.use = signac.linkpeaks.genes.use,
      n_sample = signac.linkpeaks.n_sample,
      pvalue_cutoff = signac.linkpeaks.pvalue_cutoff,
      score_cutoff = signac.linkpeaks.score_cutoff,
      gene.id = signac.linkpeaks.gene.id)

    clusters <- sort(as.numeric(as.character(unique(top1$cluster))),decreasing = F)
    DefaultAssay(current_seurat) <- "peaks"
    for(k in 1:length(clusters)){
      if(length(which(unlist(top1$cluster) == clusters[k])) > 0){
        p <- NULL
        p <- CoveragePlot(
          object = current_seurat,
          region = unlist(top1[which(top1$cluster == clusters[k]),"gene"]),
          features = unlist(top1_derivedrna[which(top1_derivedrna$cluster == clusters[k]),"gene"]),
          expression.assay = "ACTIVITY",
          extend.upstream = 10000,
          extend.downstream = 10000)+
          plot_annotation(title = paste("TOP PEAK AND GENE ACTIVITY\n",pheno_data[i,"SID"], ": CLUSTER ", clusters[k], sep = ""),theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
        
        # results[['p10plots']][[length(results[['p10plots']])+1]] <- p
        # names(results[['p10plots']])[length(results[['p10plots']])] <- paste(pheno_data[i,"SID"],"|","CLUSTER_",clusters[k],sep = "")
        somePNGPath <- paste(cdir,"10URSA_PLOT_scATAC_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP1_PEAK_TOP1_GENE_ACTIVITY_IN_CLUSTER_",clusters[k],"_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
        png(somePNGPath, width = 4000, height = 4000*5/6, units = "px", res = 300)
        print(p)
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
    # results[['p11data']][[i]] <- topn_genes
    # names(results[['p11data']])[i] <- pheno_data[i,"SID"]
    write.csv(topn_genes, paste(cdir, "11URSA_TABLE_scATAC_TOP_",n,"_GENES_IN_CLUSTERS_EASY_TABLE_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)


    print(paste("Add motif information for ", pheno_data[i,"FILE"], "..", sep = ""))
    current_assay <- "peaks"
    used_type <- "ACTIVITY"
    DefaultAssay(current_seurat) <- "peaks"
    chosen_genes <- NULL
    position_freq <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- AddMotifs(current_seurat,
                                  genome = BSgenome.Hsapiens.UCSC.hg19,
                                  pfm = position_freq,
                                  assay = signac.addmotifs.assay)
      motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
      if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
        chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
        chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
        current_seurat <- Footprint(object = current_seurat,
                                    motif.name = chosen_genes,
                                    genome = BSgenome.Hsapiens.UCSC.hg19,
                                    regions = signac.footprint.regions,
                                    assay = signac.footprint.assay)
      }
    }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- AddMotifs(current_seurat,
                                  genome = BSgenome.Hsapiens.UCSC.hg38,
                                  pfm = position_freq,
                                  assay = signac.addmotifs.assay)
      motif.names <- as.character(unlist(GetAssayData(object = current_seurat[[current_assay]], slot = "motifs")@motif.names))
      if(length(which(unique(unlist(topn_genes$gene)) %in% motif.names)) > 0){
        chosen_genes <- unique(unlist(topn_genes$gene))[which(unique(unlist(topn_genes$gene)) %in% motif.names)]
        chosen_genes <- chosen_genes[1:ifelse(length(chosen_genes) <= n, length(chosen_genes), n)]
        current_seurat <- Footprint(object = current_seurat,
                                    motif.name = chosen_genes,
                                    genome = BSgenome.Hsapiens.UCSC.hg38,
                                    regions = signac.footprint.regions,
                                    assay = signac.footprint.assay)
      }
    }

    p <- NULL
    p <- PlotFootprint(current_seurat, features = chosen_genes, group.by	= "seurat_clusters")
    p <- p + plot_layout(ncol = 2)+
      plot_annotation(title = paste(pheno_data[i,"SID"], ": TOP ",used_type," GENES \n(Top group with highest accessibility in motif flanking region are labeled)", sep = ""),
                      theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
    # results[['p12plots']][[i]] <- p
    # names(results[['p12plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"12URSA_PLOT_scATAC_",toupper(current_assay),"_TOP_",used_type,"_GENES_FOOTPRINTING_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 4000, units = "px", res = 250)
    print(p)
    dev.off()

    p_thresh <- ifelse(is.null(seurat.peaks.findallmarkers.return.thresh), 0.05, seurat.peaks.findallmarkers.return.thresh)
    
    DefaultAssay(current_seurat) <- "peaks"

    enriched.motifs <- FindMotifs(object = current_seurat,
                                  features = da_peaks[which(da_peaks$p_val_adj < p_thresh),"gene"],
                                  background = signac.findmotifs.background,
                                  assay = signac.findmotifs.assay,
                                  p.adjust.method = signac.findmotifs.p.adjust.method)
    
    # results[['p13data']][[i]] <- enriched.motifs
    # names(results[['p13data']])[i] <- pheno_data[i,"SID"]
    write.csv(enriched.motifs, paste(cdir, "13URSA_TABLE_scATAC_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SID"],".csv",sep = ""), row.names = F, quote = F)

    n <- 10
    p <- NULL
    p <- MotifPlot(object = current_seurat, motifs = rownames(enriched.motifs)[1:n])+
      plot_annotation(title = paste(pheno_data[i,"SID"], "\nPOSITION WEIGHT MATRICES OF TOP ENRICHED MOTIFS", sep = ""),
                      theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
    p <- adjust_theme(p, xsize = 12, title_size = 20, strip_size = 15)
    
    # results[['p14plots']][[i]] <- p
    # names(results[['p14plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"14URSA_PLOT_scATAC_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 2800, units = "px", res = 300)
    print(p)
    dev.off()


    print(paste("Computing peak set variability for ", pheno_data[i,"FILE"], "..", sep = ""))
    # Per-cell motif activity score
    if(length(grep("hg19",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RunChromVAR(object = current_seurat,
                                    genome = BSgenome.Hsapiens.UCSC.hg19,
                                    motif.matrix = signac.runchromvar.motif.matrix,
                                    assay = signac.runchromvar.assay,
                                    new.assay.name = "chromvar")
    }else if(length(grep("hg38|grch38",pheno_data[i,"REF_GENOME"], ignore.case = T)) > 0){
      current_seurat <- RunChromVAR(object = current_seurat,
                                    genome = BSgenome.Hsapiens.UCSC.hg38,
                                    motif.matrix = signac.runchromvar.motif.matrix,
                                    assay = signac.runchromvar.assay,
                                    new.assay.name = "chromvar")
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

     p <- p +
      plot_annotation(title = paste(pheno_data[i,"SID"], "\nPeaks VS Top Enriched Motifs", sep = ""),
                      theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
     
    # results[['p15plots']][[i]] <- p
    # names(results[['p15plots']])[i] <- pheno_data[i,"SID"]

    somePNGPath <- paste(cdir,"15URSA_PLOT_scATAC_UMAP_VS_TOP_ENRICHED_MOTIFS_",pheno_data[i,"SID"],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 2000, units = "px", res = 200)
    print(p)
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
    combined_peaks <- combined_peaks[peakwidths  < signac.post_combined_peaks.max.peaks & peakwidths > signac.post_combined_peaks.min.peaks]

    for(i in 1:length(results[['data']])){

      current <- FeatureMatrix(
        fragments = results[['data']][[i]]@assays$peaks@fragments,
        features = combined_peaks,
        cells = colnames(results[['data']][[i]]))

      current <- CreateChromatinAssay(current, 
                                      fragments = results[['data']][[i]]@assays$peaks@fragments,
                                      min.cells = signac.createchromatinassay.min.cells,
                                      min.features = signac.createchromatinassay.min.features,
                                      max.cells = signac.createchromatinassay.max.cells,
                                      ranges = signac.createchromatinassay.ranges,
                                      motifs = signac.createchromatinassay.motifs,
                                      annotation = signac.createchromatinassay.annotation,
                                      bias = signac.createchromatinassay.bias,
                                      positionEnrichment = signac.createchromatinassay.positionEnrichment,
                                      validate.fragments = signac.createchromatinassay.validate.fragments)
      
      results[['data']][[i]] <- CreateSeuratObject(current, assay = "peaks")
      results[['data']][[i]]$Sample <- sample_names[i]

    }

    rm(current)

    data <- merge(x = results[['data']][[1]], y = results[['data']][2:length(results[['data']])],
                  add.cell.ids = sample_names)

    print(paste("Running dimension reduction and clustering for merged samples..", sep = ""))

    n <- ifelse((length(data@assays$peaks@var.features)-1) < 50, (length(data@assays$peaks@var.features) - 1), 50)
    data <- RunTFIDF(data,
                     assay = signac.runtfidf.assay,
                     method = signac.runtfidf.method,
                     scale.factor = signac.runtfidf.scale.factor,
                     idf = signac.runtfidf.idf)
    
    data <- FindTopFeatures(data,
                            assay = signac.findtopfeatures.assay,
                            min.cutoff = signac.findtopfeatures.min.cutoff)
    
    data <- RunSVD(data,
                   assay = runsvd.assay,
                   features = runsvd.features,
                   n = runsvd.n,
                   reduction.key = runsvd.reduction.key,
                   reduction.name = runsvd.reduction.name,
                   scale.max = runsvd.scale.max) # , n = n
    
    dim_i <- 2 # Skip first component with high correlation
    m <- ifelse((length(data@assays$peaks@var.features)-1) < 50, (length(data@assays$peaks@var.features) - 1), 50)
    if(is.null(signac.runumap.lsi.dimension)){
      signac.runumap.lsi.dimension <- m
    }
    
    data <- RunUMAP(data,
                    reduction = 'lsi',
                    dims = dim_i:signac.runumap.lsi.dimension,
                    features = seurat.runumap.features,
                    graph = seurat.runumap.graph,
                    assay = seurat.runumap.assay,
                    nn.name = seurat.runumap.nn.name,
                    slot = seurat.runumap.slot,
                    umap.method = seurat.runumap.umap.method,
                    reduction.model = seurat.runumap.reduction.model,
                    return.model = seurat.runumap.return.model,
                    n.neighbors = seurat.runumap.n.neighbors,
                    n.components = seurat.runumap.n.components,
                    metric = seurat.runumap.metric,
                    n.epochs = seurat.runumap.n.epochs,
                    learning.rate = seurat.runumap.learning.rate,
                    min.dist = seurat.runumap.min.dist,
                    spread = seurat.runumap.spread,
                    set.op.mix.ratio = seurat.runumap.set.op.mix.ratio,
                    local.connectivity = seurat.runumap.local.connectivity,
                    repulsion.strength = seurat.runumap.repulsion.strength,
                    negative.sample.rate = seurat.runumap.negative.sample.rate,
                    a = seurat.runumap.a,
                    b = seurat.runumap.b,
                    uwot.sgd = seurat.runumap.uwot.sgd,
                    seed.use = seurat.runumap.seed.use,
                    metric.kwds = seurat.runumap.metric.kwds,
                    angular.rp.forest = seurat.runumap.angular.rp.forest,
                    densmap = seurat.runumap.densmap,
                    dens.lambda = seurat.runumap.dens.lambda,
                    dens.frac = seurat.runumap.dens.frac,
                    dens.var.shift = seurat.runumap.dens.var.shift,
                    verbose = seurat.runumap.verbose,
                    reduction.name = seurat.runumap.reduction.name,
                    reduction.key = seurat.runumap.reduction.key)
    
    # data <- RunTSNE(data, reduction = 'lsi', check_duplicates = FALSE) # dims = 2:n,
    data <- FindNeighbors(data,
                          reduction = 'lsi',
                          dims = dim_i:signac.runumap.lsi.dimension,
                          assay = seurat.findneighbors.assay,
                          features = seurat.findneighbors.features,
                          k.param = seurat.findneighbors.k.param,
                          return.neighbor = seurat.findneighbors.return.neighbor,
                          prune.SNN = seurat.findneighbors.prune.SNN,
                          nn.method = seurat.findneighbors.nn.method,
                          n.trees = seurat.findneighbors.n.trees,
                          annoy.metric = seurat.findneighbors.annoy.metric,
                          nn.eps = seurat.findneighbors.nn.eps,
                          verbose = seurat.findneighbors.verbose,
                          force.recalc = seurat.findneighbors.force.recalc,
                          do.plot = seurat.findneighbors.do.plot,
                          graph.name = seurat.findneighbors.graph.name,
                          l2.norm = seurat.findneighbors.l2.norm,
                          cache.index = seurat.findneighbors.cache.index) # dims = 2:n,
    data <- FindClusters(data,
                         graph.name = seurat.findclusters.graph.name,
                         modularity.fxn = seurat.findclusters.modularity.fxn,
                         initial.membership = seurat.findclusters.initial.membership,
                         node.sizes = seurat.findclusters.node.sizes,
                         resolution = seurat.findclusters.resolution,
                         method = seurat.findclusters.method,
                         algorithm = seurat.findclusters.algorithm,
                         n.start = seurat.findclusters.n.start,
                         n.iter = seurat.findclusters.n.iter,
                         random.seed = seurat.findclusters.random.seed,
                         group.singletons = seurat.findclusters.group.singletons,
                         temp.file.location = seurat.findclusters.temp.file.location,
                         edge.file.name = seurat.findclusters.edge.file.name,
                         verbose = seurat.findclusters.verbose)

    print(paste("Finding top peaks for merged samples..", sep = ""))
    data_markers <- FindAllMarkers(data,
                                   latent.vars = seurat.peaks.findallmarkers.latent.vars,
                                   assay = seurat.peaks.findallmarkers.assay,
                                   features = seurat.peaks.findallmarkers.features,
                                   logfc.threshold = seurat.peaks.findallmarkers.logfc.threshold,
                                   test.use = seurat.peaks.findallmarkers.test.use,
                                   slot = seurat.peaks.findallmarkers.slot,
                                   min.pct = seurat.peaks.findallmarkers.min.pct,
                                   min.diff.pct = seurat.peaks.findallmarkers.min.diff.pct,
                                   node = seurat.peaks.findallmarkers.node,
                                   verbose = seurat.peaks.findallmarkers.verbose,
                                   only.pos = seurat.peaks.findallmarkers.only.pos,
                                   max.cells.per.ident = seurat.peaks.findallmarkers.max.cells.per.ident,
                                   random.seed = seurat.peaks.findallmarkers.random.seed,
                                   latent.vars = seurat.peaks.findallmarkers.latent.vars,
                                   min.cells.feature = seurat.peaks.findallmarkers.min.cells.feature,
                                   min.cells.group = seurat.peaks.findallmarkers.min.cells.group,
                                   mean.fxn = seurat.peaks.findallmarkers.mean.fxn,
                                   fc.name = seurat.peaks.findallmarkers.fc.name,
                                   base = seurat.peaks.findallmarkers.base,
                                   return.thresh = seurat.peaks.findallmarkers.return.thresh,
                                   densify = seurat.peaks.findallmarkers.densify)

    if(toupper(peaks.findallmarkers.multipletestingcorrection) != toupper("bonferroni")){
      data_markers$p_val_adj <- p.adjust(data_markers$p_val, method = peaks.findallmarkers.multipletestingcorrection)
    }
    
    data_markers <- data_markers[data_markers$p_val_adj < peaks.findallmarkers.padjust_threshold,]
    data_markers <- data_markers[order(data_markers$p_val_adj, decreasing = F),]
    results[['data_markers']] <- data_markers
    write.csv(data_markers, paste(cdir, "16URSA_TABLE_scATAC_TOP_PEAKS_MERGED_SAMPLES.csv",sep = ""), row.names = F, quote = F)

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

    # results[['p16plots']] <- p1_umap+p2_umap

    somePNGPath <- paste(cdir,"16URSA_PLOT_scATAC_UMAP_INTEGRATED_ATAC_SAMPLES_AUTOCLUST_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = 1800, units = "px", res = 200)
    print(p1_umap+p2_umap)
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

    # results[['p17plots']] <- plots
    somePNGPath <- paste(cdir,"17URSA_PLOT_scATAC_INTEGRATED_ATAC_SAMPLES_UMAP_DENSITY_TOP_1_MARKER_IN_CLUSTERS_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height = ceiling(length(top1$gene)/4)*1000, units = "px", res = 300)
    print(do.call("grid.arrange", c(plots, ncol=4)))
    dev.off()

    n <- 10
    wt <- colnames(results[['data_markers']])[grep("log.*FC", colnames(results[['data_markers']]), ignore.case = T)]
    top_n <- results[['data_markers']] %>% group_by(cluster) %>% top_n(n = n, wt = eval(parse(text = wt)))
    current_clusters <- as.numeric(as.character(sort(unique(results[['data_markers']]$cluster))))
    for(k in 1:length(current_clusters)){
      if(length(which(unlist(top_n$cluster) == current_clusters[k])) > 0){
        current <- top_n[which(top_n$cluster == current_clusters[k]),]
        current <- current[order(current$p_val_adj, decreasing = F),]
        p <- NULL
        p <- VlnPlot(data, features = current$gene, pt.size = 0,
                                                                            # log = TRUE,
                                                                            cols = gen_colors(color_conditions$tenx,length(unique(current_clusters)))) +
          plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", current_clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
        # results[['p18plots']][[length(results[['p18plots']])+1]]
        # names(results[['p18plots']])[length(results[['p18plots']])] <- paste("CLUSTER_",current_clusters[k], sep = "")
        somePNGPath <- paste(cdir,"18URSA_PLOT_scATAC_INTEGRATED_ATAC_SAMPLES_TOP_",n,"_MARKERS_IN_CLUSTER_",current_clusters[k],"_",project_name,".png", sep = "")
        png(somePNGPath, width = 4000, height = ceiling(length(current_clusters)/4)*800, units = "px", res = 300)
        print(p)
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
    gene_activities <- GeneActivity(data,
                                    assay = seurat.geneactivity.assay,
                                    features = seurat.geneactivity.features,
                                    extend.upstream = seurat.geneactivity.extend.upstream,
                                    extend.downstream = seurat.geneactivity.extend.downstream,
                                    biotypes = seurat.geneactivity.biotypes,
                                    max.width = seurat.geneactivity.max.width,
                                    process_n = seurat.geneactivity.process_n,
                                    gene.id = seurat.geneactivity.gene.id)
    data[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities,
                                            min.cells = seurat.createassayobject.min.cells,
                                            min.features = seurat.createassayobject.min.features,
                                            check.matrix = seurat.createassayobject.check.matrix)
    rm(gene_activities)

    DefaultAssay(data) <- "ACTIVITY"

    print(paste("Finding top genes for merged samples..", sep = ""))
    data_activity_markers <- FindAllMarkers(data,
                                            assay = seurat.geneactivity.findallmarkers.assay,
                                            features = seurat.geneactivity.findallmarkers.features,
                                            logfc.threshold = seurat.geneactivity.findallmarkers.logfc.threshold,
                                            test.use = seurat.geneactivity.findallmarkers.test.use,
                                            slot = seurat.geneactivity.findallmarkers.slot,
                                            min.pct = seurat.geneactivity.findallmarkers.min.pct,
                                            min.diff.pct = seurat.geneactivity.findallmarkers.min.diff.pct,
                                            node = seurat.geneactivity.findallmarkers.node,
                                            verbose = seurat.geneactivity.findallmarkers.verbose,
                                            only.pos = seurat.geneactivity.findallmarkers.only.pos,
                                            max.cells.per.ident = seurat.geneactivity.findallmarkers.max.cells.per.ident,
                                            random.seed = seurat.geneactivity.findallmarkers.random.seed,
                                            latent.vars = seurat.geneactivity.findallmarkers.latent.vars,
                                            min.cells.feature = seurat.geneactivity.findallmarkers.min.cells.feature,
                                            min.cells.group = seurat.geneactivity.findallmarkers.min.cells.group,
                                            mean.fxn = seurat.geneactivity.findallmarkers.mean.fxn,
                                            fc.name = seurat.geneactivity.findallmarkers.fc.name,
                                            base = seurat.geneactivity.findallmarkers.base,
                                            return.thresh = seurat.geneactivity.findallmarkers.return.thresh,
                                            densify = seurat.geneactivity.findallmarkers.densify)
    if(toupper(geneactivity.findallmarkers.multipletestingcorrection) != toupper("bonferroni")){
      data_activity_markers$p_val_adj <- p.adjust(data_activity_markers$p_val, method = geneactivity.findallmarkers.multipletestingcorrection)
    }

    data_activity_markers <- data_activity_markers[data_activity_markers$p_val_adj < geneactivity.findallmarkers.padjust_threshold,]
    data_activity_markers <- data_activity_markers[order(data_activity_markers$p_val_adj, decreasing = F),]
    # results[['data_activity_markers']] <- data_activity_markers
    write.csv(data_activity_markers, paste(cdir, "18URSA_TABLE_scATAC_TOP_GENES_MERGED_SAMPLES.csv",sep = ""), row.names = F, quote = F)

    wt <- colnames(data_activity_markers)[grep("log.*FC", colnames(data_activity_markers), ignore.case = T)]
    top1_activities <- data_activity_markers %>% group_by(cluster) %>% top_n(n = 1, wt = eval(parse(text = wt)))
    top1_activities <- top1_activities[order(top1_activities$cluster, decreasing = F),]

    top1 <- data_activity_markers %>% group_by(cluster) %>% dplyr::filter((eval(parse(text = wt)) > log2(1.5))) %>% dplyr::slice(1:1) # %>% filter(avg_log2FC > 1)
    top1 <- top1[order(top1$cluster, decreasing = F),]

    clusters <- unique(sort(as.numeric(as.character(top1$cluster))))

    DefaultAssay(data) <- "peaks"
    for(k in 1:length(clusters)){
      if(length(which(unlist(top1$cluster) == clusters[k])) > 0){
        p <- NULL
        p <- CoveragePlot(
          object = data,
          expression.assay = "ACTIVITY",
          region = unlist(top1[which(top1$cluster == clusters[k]),"gene"]),
          features = unlist(top1_activities[which(top1_activities$cluster == clusters[k]),"gene"]),
          extend.upstream = 10000,
          extend.downstream = 10000)+
          # scale_color_manual(values = gen_colors(color_conditions$tenx, length(unique(data$Sample))))+
          plot_annotation(title = paste("INTEGRATED DATA: CLUSTER ", clusters[k], sep = ""),
                          theme = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")))
        
        # results[['p19plots']][length(results[['p19plots']])+1] <- p
        # names(results[['p19plots']])[length(results[['p19plots']])] <- paste("CLUSTER_",clusters[k], sep = "")
        somePNGPath <- paste(cdir,"19URSA_PLOT_scATAC_INTEGRATED_SAMPLES_Tn5_INSERTION_COVERAGE_TOP_PEAK_GENE_ACTIVITY_IN_CLUSTER_",current_clusters[k],"_",project_name,".png", sep = "")
        png(somePNGPath, width = 3000, height = 2000, units = "px", res = 300)
        print(p)
        dev.off()
      }
    }
  }
  saveRDS(data, paste(cdir,"20URSA_INTEGRATED_DATA_scATAC_",project_name,".RDS", sep = ""))
  saveRDS(results[['data']], paste(cdir,"20URSA_SAMPLE_DATA_scATAC_",project_name,".RDS", sep = ""))
  
  print("Completed!")
}

