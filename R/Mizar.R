############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Mizar: Spatial
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
#' @param multipletestingcorrection Correction method for multiple testing.
#' Default to "bonferroni". Methods include: "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none". (refer to the p.adjust() function
#' in stats package for details).
#' @param padjust_threshold Threshold for retaining entries with corrected
#' p-values. Default to 0.1.
#' @param rnaseq_dir scRNA-seq data files directory. We recommend putting the scRNA-seq
#' data in a separate directory separating them from spatial input files.
#' scRNA-seq files should be filtered .h5 format files from 10X Genomics and their file
#' names should have the same sample ID as file prefix to their correspond spatial files.
#' @param run_rnaseq If scRNA-seq integration with spatial data should be run. Default to FALSE.
#' There is no need to prepare scRNA-seq folder if run_rnaseq is set to FALSE.
#' @param seurat.scrnaseq.createseuratobject.min.cells If run_rnaseq is set to TRUE, the quality
#' control filter for minimum number of cells with features detected. Default to 3 (refer to the
#' subset() function in Seurat package for details).
#' @param seurat.scrnaseq.subset.min.nFeature_RNA If run_rnaseq is set to TRUE, the quality
#' control filter for cells with at least this number of features detected.
#' Default to 200 (refer to the subset() function in Seurat package for details).
#' @param seurat.scrnaseq.subset.max.nFeature_RNA If run_rnaseq is set to TRUE, the quality
#' control filter for cells with at most this number of features detected.
#' Default to 6000 (refer to the subset() function in Seurat package for details).
#' @param seurat.scrnaseq.subset.max.percent.mt If run_rnaseq is set to TRUE, the quality
#' control filter for cells high proportion of mitochondrial reads.
#' Default to 10.
#' @param ... Arguments passed to other parameters in the dependency pages.
#' Parameters with the long format: xxx.xxx.xxx usually indicates in lowercases
#' the parameter origin: <dependent package name>.<function name>.<parameter name>(),
#' for example: seurat.load10x_spatial.assay() indicates this parameter 
#' originates from the package Seurat under the function Load10X_Spatial() and its
#' parameter 'assay'. Users could subsequently refer to the dependency 
#' package for detailed information on parameters and their usage.
#' @export
#'

SpatialPip <- function(project_name = "Ursa_Spatial",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file,
                     rnaseq_dir = "./",
                     run_rnaseq = FALSE,
                     multipletestingcorrection = "bonferroni",
                     padjust_threshold = 0.1,
                     
                     # Load10X_Spatial
                     seurat.load10x_spatial.assay = "Spatial",
                     seurat.load10x_spatial.slice = "slice1",
                     seurat.load10x_spatial.filter.matrix = TRUE,
                     seurat.load10x_spatial.to.upper = FALSE,
                     seurat.load10x_spatial.image = NULL,
                     
                     # SCTransform
                     seurat.sctransform.new.assay.name = "SCT",
                     seurat.sctransform.reference.SCT.model = NULL,
                     seurat.sctransform.do.correct.umi = TRUE,
                     seurat.sctransform.ncells = 5000,
                     seurat.sctransform.residual.features = NULL,
                     seurat.sctransform.variable.features.n = 3000,
                     seurat.sctransform.variable.features.rv.th = 1.3,
                     seurat.sctransform.vars.to.regress = NULL,
                     seurat.sctransform.do.scale = FALSE,
                     seurat.sctransform.do.center = TRUE,
                     # seurat.sctransform.clip.range = c(-sqrt(x = ncol(x = object[[assay]])/30), sqrt(x = ncol(x = object[[assay]])/30)),
                     seurat.sctransform.conserve.memory = FALSE,
                     seurat.sctransform.return.only.var.genes = TRUE,
                     seurat.sctransform.seed.use = 1448145,
                     
                     # FindVariableFeatures
                     seurat.findvariablefeatures.selection.method = "vst",
                     seurat.findvariablefeatures.loess.span = 0.3,
                     seurat.findvariablefeatures.clip.max = "auto",
                     # seurat.findvariablefeatures.mean.function = FastExpMean,
                     # seurat.findvariablefeatures.dispersion.function = FastLogVMR,
                     seurat.findvariablefeatures.num.bin = 20,
                     seurat.findvariablefeatures.binning.method = "equal_width",

                     # RunPCA
                     seurat.runpca.assay = NULL,
                     seurat.runpca.features = NULL,
                     seurat.runpca.npcs = 50,
                     seurat.runpca.rev.pca = FALSE,
                     seurat.runpca.weight.by.var = TRUE,
                     seurat.runpca.ndims.print = 1:5,
                     seurat.runpca.nfeatures.print = 30,
                     seurat.runpca.reduction.name = "pca",
                     seurat.runpca.reduction.key = "PC_",
                     seurat.runpca.seed.use = 42,
                     
                     # RunUMAP
                     seurat.runumap.pcs = 30,
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
                     seurat.runumap.reduction.name = "umap",
                     seurat.runumap.reduction.key = "UMAP_",
                     
                     # RunTSNE
                     seurat.runtsne.cells = NULL,
                     seurat.runtsne.features = NULL,
                     seurat.runtsne.seed.use = 1,
                     seurat.runtsne.tsne.method = "Rtsne",
                     seurat.runtsne.dim.embed = 2,
                     seurat.runtsne.distance.matrix = NULL,
                     seurat.runtsne.reduction.name = "tsne",
                     seurat.runtsne.reduction.key = "tSNE_",
                     seurat.runtsne.check_duplicates = FALSE,
                     
                     # FindNeighbors
                     seurat.findneighbors.assay = NULL,
                     seurat.findneighbors.features = NULL,
                     seurat.findneighbors.k.param = 20,
                     seurat.findneighbors.return.neighbor = FALSE,
                     # seurat.findneighbors.compute.SNN = !return.neighbor,
                     seurat.findneighbors.prune.SNN = 1/15,
                     seurat.findneighbors.nn.method = "annoy",
                     seurat.findneighbors.n.trees = 50,
                     seurat.findneighbors.annoy.metric = "euclidean",
                     seurat.findneighbors.nn.eps = 0,
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
                     seurat.findclusters.algorithm = 1,
                     seurat.findclusters.n.start = 10,
                     seurat.findclusters.n.iter = 10,
                     seurat.findclusters.random.seed = 0,
                     seurat.findclusters.group.singletons = TRUE,
                     seurat.findclusters.temp.file.location = NULL,
                     seurat.findclusters.edge.file.name = NULL,

                     # FindAllMarkers
                     seurat.findallmarkers.assay = NULL,
                     seurat.findallmarkers.features = NULL,
                     seurat.findallmarkers.logfc.threshold = 0.25,
                     seurat.findallmarkers.test.use = 'wilcox',
                     seurat.findallmarkers.slot = "data",
                     seurat.findallmarkers.min.pct = 0.25,
                     seurat.findallmarkers.min.diff.pct = -Inf,
                     seurat.findallmarkers.node = NULL,
                     seurat.findallmarkers.only.pos = FALSE,
                     seurat.findallmarkers.max.cells.per.ident = Inf,
                     seurat.findallmarkers.random.seed = 1,
                     seurat.findallmarkers.latent.vars = NULL,
                     seurat.findallmarkers.min.cells.feature = 3,
                     seurat.findallmarkers.min.cells.group = 3,
                     seurat.findallmarkers.mean.fxn = NULL,
                     seurat.findallmarkers.fc.name = NULL,
                     seurat.findallmarkers.base = 2,
                     seurat.findallmarkers.return.thresh = 0.05,
                     seurat.findallmarkers.densify = FALSE,
                     
                     # SingleR
                     singler.ref = "HumanPrimaryCellAtlas",
                     singler.labels = "HumanPrimaryCellAtlasLevel2",
                     singler.method = NULL,
                     singler.genes = "de",
                     singler.sd.thresh = 1,
                     singler.de.method = "classic",
                     singler.de.n = NULL,
                     singler.de.args = list(),
                     singler.aggr.ref = FALSE,
                     singler.aggr.args = list(),
                     singler.recompute = TRUE,
                     singler.restrict = NULL,
                     singler.quantile = 0.8,
                     singler.fine.tune = TRUE,
                     singler.tune.thresh = 0.05,
                     singler.prune = TRUE,
                     singler.assay.type.test = "logcounts",
                     singler.assay.type.ref = "logcounts",
                     singler.check.missing = TRUE,
                     # singler.num.threads = bpnworkers(BPPARAM),
                     # singler.BNPARAM = NULL,
                     # singler.BPPARAM = SerialParam(),
                     
                     # CreateSeuratObject
                     seurat.scrnaseq.createseuratobject.min.cells = 3,

                     # subset - scRNA-seq
                     seurat.scrnaseq.subset.min.nFeature_RNA = 200,
                     seurat.scrnaseq.subset.max.nFeature_RNA = 6000,
                     seurat.scrnaseq.subset.max.percent.mt = 10,
                     
                     # SCTransform
                     seurat.scrnaseq.sctransform.assay = "RNA",
                     seurat.scrnaseq.sctransform.new.assay.name = "SCT",
                     seurat.scrnaseq.sctransform.reference.SCT.model = NULL,
                     seurat.scrnaseq.sctransform.do.correct.umi = TRUE,
                     seurat.scrnaseq.sctransform.ncells = 3000,
                     seurat.scrnaseq.sctransform.residual.features = NULL,
                     seurat.scrnaseq.sctransform.variable.features.n = 3000,
                     seurat.scrnaseq.sctransform.variable.features.rv.th = 1.3,
                     seurat.scrnaseq.sctransform.vars.to.regress = NULL,
                     seurat.scrnaseq.sctransform.do.scale = FALSE,
                     seurat.scrnaseq.sctransform.do.center = TRUE,
                     # seurat.scrnaseq.sctransform.clip.range = c(-sqrt(x = ncol(x = object[[assay]])/30), sqrt(x = ncol(x = object[[assay]])/30)),
                     seurat.scrnaseq.sctransform.conserve.memory = FALSE,
                     seurat.scrnaseq.sctransform.return.only.var.genes = TRUE,
                     seurat.scrnaseq.sctransform.seed.use = 1448145,
                     
                     # RunPCA
                     seurat.scrnaseq.runpca.assay = NULL,
                     seurat.scrnaseq.runpca.features = NULL,
                     seurat.scrnaseq.runpca.npcs = 50,
                     seurat.scrnaseq.runpca.rev.pca = FALSE,
                     seurat.scrnaseq.runpca.weight.by.var = TRUE,
                     seurat.scrnaseq.runpca.ndims.print = 1:5,
                     seurat.scrnaseq.runpca.nfeatures.print = 30,
                     seurat.scrnaseq.runpca.reduction.name = "pca",
                     seurat.scrnaseq.runpca.reduction.key = "PC_",
                     seurat.scrnaseq.runpca.seed.use = 42,
                     
                     # RunUMAP
                     seurat.scrnaseq.runumap.pcs = 30,
                     seurat.scrnaseq.runumap.features = NULL,
                     seurat.scrnaseq.runumap.graph = NULL,
                     seurat.scrnaseq.runumap.assay = NULL,
                     seurat.scrnaseq.runumap.nn.name = NULL,
                     seurat.scrnaseq.runumap.slot = "data",
                     seurat.scrnaseq.runumap.umap.method = "uwot",
                     seurat.scrnaseq.runumap.reduction.model = NULL,
                     seurat.scrnaseq.runumap.return.model = FALSE,
                     seurat.scrnaseq.runumap.n.neighbors = 30L,
                     seurat.scrnaseq.runumap.n.components = 2L,
                     seurat.scrnaseq.runumap.metric = "cosine",
                     seurat.scrnaseq.runumap.n.epochs = NULL,
                     seurat.scrnaseq.runumap.learning.rate = 1,
                     seurat.scrnaseq.runumap.min.dist = 0.3,
                     seurat.scrnaseq.runumap.spread = 1,
                     seurat.scrnaseq.runumap.set.op.mix.ratio = 1,
                     seurat.scrnaseq.runumap.local.connectivity = 1L,
                     seurat.scrnaseq.runumap.repulsion.strength = 1,
                     seurat.scrnaseq.runumap.negative.sample.rate = 5L,
                     seurat.scrnaseq.runumap.a = NULL,
                     seurat.scrnaseq.runumap.b = NULL,
                     seurat.scrnaseq.runumap.uwot.sgd = FALSE,
                     seurat.scrnaseq.runumap.seed.use = 42L,
                     seurat.scrnaseq.runumap.metric.kwds = NULL,
                     seurat.scrnaseq.runumap.angular.rp.forest = FALSE,
                     seurat.scrnaseq.runumap.densmap = FALSE,
                     seurat.scrnaseq.runumap.dens.lambda = 2,
                     seurat.scrnaseq.runumap.dens.frac = 0.3,
                     seurat.scrnaseq.runumap.dens.var.shift = 0.1,
                     seurat.scrnaseq.runumap.reduction.name = "umap",
                     seurat.scrnaseq.runumap.reduction.key = "UMAP_",
                     
                     # RunTSNE
                     seurat.scrnaseq.runtsne.cells = NULL,
                     seurat.scrnaseq.runtsne.features = NULL,
                     seurat.scrnaseq.runtsne.seed.use = 1,
                     seurat.scrnaseq.runtsne.tsne.method = "Rtsne",
                     seurat.scrnaseq.runtsne.dim.embed = 2,
                     seurat.scrnaseq.runtsne.distance.matrix = NULL,
                     seurat.scrnaseq.runtsne.reduction.name = "tsne",
                     seurat.scrnaseq.runtsne.reduction.key = "tSNE_",
                     seurat.scrnaseq.runtsne.check_duplicates = FALSE,
                     
                     # Seurat - FindNeighbors
                     seurat.scrnaseq.findneighbors.assay = NULL,
                     seurat.scrnaseq.findneighbors.features = NULL,
                     seurat.scrnaseq.findneighbors.k.param = 20,
                     seurat.scrnaseq.findneighbors.return.neighbor = FALSE,
                     # seurat.scrnaseq.findneighbors.compute.SNN = !return.neighbor,
                     seurat.scrnaseq.findneighbors.prune.SNN = 1/15,
                     seurat.scrnaseq.findneighbors.nn.method = "annoy",
                     seurat.scrnaseq.findneighbors.n.trees = 50,
                     seurat.scrnaseq.findneighbors.annoy.metric = "euclidean",
                     seurat.scrnaseq.findneighbors.nn.eps = 0,
                     seurat.scrnaseq.findneighbors.force.recalc = FALSE,
                     seurat.scrnaseq.findneighbors.do.plot = FALSE,
                     seurat.scrnaseq.findneighbors.graph.name = NULL,
                     seurat.scrnaseq.findneighbors.l2.norm = FALSE,
                     seurat.scrnaseq.findneighbors.cache.index = FALSE,
                     
                     # Seurat - FindClusters
                     seurat.scrnaseq.findclusters.graph.name = NULL,
                     seurat.scrnaseq.findclusters.modularity.fxn = 1,
                     seurat.scrnaseq.findclusters.initial.membership = NULL,
                     seurat.scrnaseq.findclusters.node.sizes = NULL,
                     seurat.scrnaseq.findclusters.resolution = 0.8,
                     seurat.scrnaseq.findclusters.method = "matrix",
                     seurat.scrnaseq.findclusters.algorithm = 1,
                     seurat.scrnaseq.findclusters.n.start = 10,
                     seurat.scrnaseq.findclusters.n.iter = 10,
                     seurat.scrnaseq.findclusters.random.seed = 0,
                     seurat.scrnaseq.findclusters.group.singletons = TRUE,
                     seurat.scrnaseq.findclusters.temp.file.location = NULL,
                     seurat.scrnaseq.findclusters.edge.file.name = NULL,
                     
                     # SingleR
                     singler.scrnaseq.ref = "HumanPrimaryCellAtlas",
                     singler.scrnaseq.labels = "HumanPrimaryCellAtlasLevel1",
                     singler.scrnaseq.method = NULL,
                     singler.scrnaseq.genes = "de",
                     singler.scrnaseq.sd.thresh = 1,
                     singler.scrnaseq.de.method = "classic",
                     singler.scrnaseq.de.n = NULL,
                     singler.scrnaseq.de.args = list(),
                     singler.scrnaseq.aggr.ref = FALSE,
                     singler.scrnaseq.aggr.args = list(),
                     singler.scrnaseq.recompute = TRUE,
                     singler.scrnaseq.restrict = NULL,
                     singler.scrnaseq.quantile = 0.8,
                     singler.scrnaseq.fine.tune = TRUE,
                     singler.scrnaseq.tune.thresh = 0.05,
                     singler.scrnaseq.prune = TRUE,
                     singler.scrnaseq.assay.type.test = "logcounts",
                     singler.scrnaseq.assay.type.ref = "logcounts",
                     singler.scrnaseq.check.missing = TRUE,
                     # singler.scrnaseq.num.threads = bpnworkers(BPPARAM),
                     # singler.scrnaseq.BNPARAM = NULL,
                     # singler.scrnaseq.BPPARAM = SerialParam(),
                     
                     # FindTransferAnchors
                     seurat.multimodalintegration.findfransferanchors.normalization.method = "SCT",
                     seurat.multimodalintegration.findfransferanchors.recompute.residuals = TRUE,
                     seurat.multimodalintegration.findfransferanchors.reference.assay = NULL,
                     seurat.multimodalintegration.findfransferanchors.reference.neighbors = NULL,
                     seurat.multimodalintegration.findfransferanchors.query.assay = NULL,
                     seurat.multimodalintegration.findfransferanchors.reduction = "pcaproject",
                     seurat.multimodalintegration.findfransferanchors.reference.reduction = NULL,
                     seurat.multimodalintegration.findfransferanchors.project.query = FALSE,
                     seurat.multimodalintegration.findfransferanchors.features = NULL,
                     seurat.multimodalintegration.findfransferanchors.scale = TRUE,
                     seurat.multimodalintegration.findfransferanchors.npcs = 30,
                     seurat.multimodalintegration.findfransferanchors.l2.norm = TRUE,
                     seurat.multimodalintegration.findfransferanchors.k.anchor = 5,
                     seurat.multimodalintegration.findfransferanchors.k.filter = 200,
                     seurat.multimodalintegration.findfransferanchors.k.score = 30,
                     seurat.multimodalintegration.findfransferanchors.max.features = 200,
                     seurat.multimodalintegration.findfransferanchors.nn.method = "annoy",
                     seurat.multimodalintegration.findfransferanchors.n.trees = 50,
                     seurat.multimodalintegration.findfransferanchors.eps = 0,
                     seurat.multimodalintegration.findfransferanchors.approx.pca = TRUE,
                     seurat.multimodalintegration.findfransferanchors.mapping.score.k = NULL,
                     
                     # TransferData
                     seurat.multimodalintegration.transferdata.reference = NULL,
                     seurat.multimodalintegration.transferdata.query = NULL,
                     seurat.multimodalintegration.transferdata.weight.reduction = "pcaproject",
                     seurat.multimodalintegration.transferdata.l2.norm = FALSE,
                     seurat.multimodalintegration.transferdata.k.weight = 50,
                     seurat.multimodalintegration.transferdata.sd.weight = 1,
                     seurat.multimodalintegration.transferdata.eps = 0,
                     seurat.multimodalintegration.transferdata.n.trees = 50,
                     seurat.multimodalintegration.transferdata.slot = "data",
                     seurat.multimodalintegration.transferdata.prediction.assay = TRUE,
                     seurat.multimodalintegration.transferdata.store.weights = TRUE){
  
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "SPATIAL", isDir = T)
  color_conditions <- color_ini()
  ctime <- time_ini()
  hpca.se <- HumanPrimaryCellAtlasData()
  if(toupper(singler.ref) == toupper("HumanPrimaryCellAtlas")){
    singler.ref <- hpca.se
    if(toupper(singler.labels) == toupper("HumanPrimaryCellAtlasLevel1")){
      singler.labels <- hpca.se$label.main
    }else if(toupper(singler.labels) == toupper("HumanPrimaryCellAtlasLevel2")){
      singler.labels <- hpca.se$label.fine
    }
  }
  
  if(toupper(singler.scrnaseq.ref) == toupper("HumanPrimaryCellAtlas")){
    singler.scrnaseq.ref <- hpca.se
    if(toupper(singler.scrnaseq.labels) == toupper("HumanPrimaryCellAtlasLevel1")){
      singler.scrnaseq.labels <- hpca.se$label.main
    }else if(toupper(singler.scrnaseq.labels) == toupper("HumanPrimaryCellAtlasLevel2")){
      singler.scrnaseq.labels <- hpca.se$label.fine
    }
  }
  
  hs <- org.Hs.eg.db
  hgnc.table <- data("hgnc.table", package="HGNChelper")
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir, recursive = T)

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
    current <- Load10X_Spatial(data.dir = gsub("(.*\\/).*","\\1",sample_files[i], ignore.case = T),
                               filename = current_sample,
                               assay = seurat.load10x_spatial.assay,
                               slice = seurat.load10x_spatial.slice,
                               filter.matrix = seurat.load10x_spatial.filter.matrix,
                               to.upper = seurat.load10x_spatial.to.upper,
                               image = seurat.load10x_spatial.image)
    current <- SCTransform(current,
                           assay = "Spatial",
                           new.assay.name = seurat.sctransform.new.assay.name,
                           reference.SCT.model = seurat.sctransform.reference.SCT.model,
                           do.correct.umi = seurat.sctransform.do.correct.umi,
                           ncells = seurat.sctransform.ncells,
                           residual.features = seurat.sctransform.residual.features,
                           variable.features.n = seurat.sctransform.variable.features.n,
                           variable.features.rv.th = seurat.sctransform.variable.features.rv.th,
                           vars.to.regress = seurat.sctransform.vars.to.regress,
                           do.scale = seurat.sctransform.do.scale,
                           do.center = seurat.sctransform.do.center,
                           # clip.range = seurat.sctransform.clip.range,
                           conserve.memory = seurat.sctransform.conserve.memory,
                           return.only.var.genes = seurat.sctransform.return.only.var.genes,
                           seed.use = seurat.sctransform.seed.use,
                           verbose = FALSE)
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
    current <- FindVariableFeatures(current,
                                    selection.method = seurat.findvariablefeatures.selection.method,
                                    loess.span = seurat.findvariablefeatures.loess.span,
                                    clip.max = seurat.findvariablefeatures.clip.max,
                                    # mean.function = seurat.findvariablefeatures.mean.function,
                                    # dispersion.function = seurat.findvariablefeatures.dispersion.function,
                                    num.bin = seurat.findvariablefeatures.num.bin,
                                    binning.method = seurat.findvariablefeatures.binning.method)
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
    current <- RunPCA(current,
                      assay = seurat.runpca.assay,
                      features = seurat.runpca.features,
                      npcs = seurat.runpca.npcs,
                      rev.pca = seurat.runpca.rev.pca,
                      weight.by.var = seurat.runpca.weight.by.var,
                      ndims.print = seurat.runpca.ndims.print,
                      nfeatures.print = seurat.runpca.nfeatures.print,
                      reduction.name = seurat.runpca.reduction.name,
                      reduction.key = seurat.runpca.reduction.key,
                      seed.use = seurat.runpca.seed.use,
                      verbose = FALSE)
    
    selected_pcs <- NULL
    selected_pcs <- find_selected_pcs(current,
                                      pcs = seurat.runumap.pcs,
                                      pca = "pca",
                                      cname = cname)
    current <- RunUMAP(current,
                       reduction = "pca",
                       dims = 1:selected_pcs,
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
    reduction.name = seurat.runumap.reduction.name,
    reduction.key = seurat.runumap.reduction.key)
    
    current <- RunTSNE(current,
                       reduction = "pca",
                       dims = 1:selected_pcs,
                       cells = seurat.runtsne.cells,
                       features = seurat.runtsne.features,
                       seed.use = seurat.runtsne.seed.use,
                       tsne.method = seurat.runtsne.tsne.method,
                       dim.embed = seurat.runtsne.dim.embed,
                       distance.matrix = seurat.runtsne.distance.matrix,
                       reduction.name = seurat.runtsne.reduction.name,
                       reduction.key = seurat.runtsne.reduction.key,
                       check_duplicates = seurat.runtsne.check_duplicates)
    current <- FindNeighbors(current,
                             reduction = "pca",
                             dims = 1:selected_pcs,
                             assay = seurat.findneighbors.assay,
                             features = seurat.findneighbors.features,
                             k.param = seurat.findneighbors.k.param,
                             return.neighbor = seurat.findneighbors.return.neighbor,
                             # compute.SNN = seurat.findneighbors.compute.SNN,
                             prune.SNN = seurat.findneighbors.prune.SNN,
                             nn.method = seurat.findneighbors.nn.method,
                             n.trees = seurat.findneighbors.n.trees,
                             annoy.metric = seurat.findneighbors.annoy.metric,
                             nn.eps = seurat.findneighbors.nn.eps,
                             force.recalc = seurat.findneighbors.force.recalc,
                             do.plot = seurat.findneighbors.do.plot,
                             graph.name = seurat.findneighbors.graph.name,
                             l2.norm = seurat.findneighbors.l2.norm,
                             cache.index = seurat.findneighbors.cache.index)
    current <- FindClusters(current,
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
                            verbose = FALSE)

    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "SCT"
    current_de_markers <- NULL
    current_de_markers <- FindAllMarkers(current,
                                         assay = seurat.findallmarkers.assay,
                                         features = seurat.findallmarkers.features,
                                         logfc.threshold = seurat.findallmarkers.logfc.threshold,
                                         test.use = seurat.findallmarkers.test.use,
                                         slot = seurat.findallmarkers.slot,
                                         min.pct = seurat.findallmarkers.min.pct,
                                         min.diff.pct = seurat.findallmarkers.min.diff.pct,
                                         node = seurat.findallmarkers.node,
                                         only.pos = seurat.findallmarkers.only.pos,
                                         max.cells.per.ident = seurat.findallmarkers.max.cells.per.ident,
                                         random.seed = seurat.findallmarkers.random.seed,
                                         latent.vars = seurat.findallmarkers.latent.vars,
                                         min.cells.feature = seurat.findallmarkers.min.cells.feature,
                                         min.cells.group = seurat.findallmarkers.min.cells.group,
                                         mean.fxn = seurat.findallmarkers.mean.fxn,
                                         fc.name = seurat.findallmarkers.fc.name,
                                         base = seurat.findallmarkers.base,
                                         return.thresh = seurat.findallmarkers.return.thresh,
                                         densify = seurat.findallmarkers.densify)
    current_de_markers <- data.frame(SAMPLE = annot_names[i], current_de_markers)
    if(toupper(multipletestingcorrection) != toupper("bonferroni")){
      current_de_markers$p_val_adj <- p.adjust(current_de_markers$p_val, method = multipletestingcorrection)
    }
    current_de_markers <- current_de_markers[which(current_de_markers$p_val_adj < padjust_threshold),]
    current_de_markers <- current_de_markers[order(current_de_markers$avg_log2FC, decreasing = T),]

    data_markers[[i]] <- current_de_markers
    names(data_markers)[i] <- annot_names[i]

    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "SCT"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(current)),
                       clusters = current$seurat_clusters,
                       ref = singler.ref,
                       labels = singler.labels,
                       method = singler.method,
                       genes = singler.genes,
                       sd.thresh = singler.sd.thresh,
                       de.method = singler.de.method,
                       de.n = singler.de.n,
                       de.args = singler.de.args,
                       aggr.ref = singler.aggr.ref,
                       aggr.args = singler.aggr.args,
                       recompute = singler.recompute,
                       restrict = singler.restrict,
                       quantile = singler.quantile,
                       fine.tune = singler.fine.tune,
                       tune.thresh = singler.tune.thresh,
                       prune = singler.prune,
                       assay.type.test = singler.assay.type.test,
                       assay.type.ref = singler.assay.type.ref,
                       check.missing = singler.check.missing)
                       # num.threads = singler.num.threads
                       # BNPARAM = singler.BNPARAM
                       # BPPARAM = singler.BPPARAM
    
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
      current <- NULL
      current <- data_current[[which(names(data_current) == pheno_data[i,"SAMPLE_ID"])]]
      if((length(cfile) > 0)){
      scrna_data <- Read10X_h5(filename = cfile)
      if((length(scrna_data) == 1) | (length(scrna_data) > 10)){
        scrna_data <- CreateSeuratObject(counts = scrna_data, project = pheno_data[i,"SAMPLE_ID"], min.cells = seurat.scrnaseq.createseuratobject.min.cells)
      }else{
        scrna_data <- CreateSeuratObject(counts = scrna_data$`Gene Expression`, project = pheno_data[i,"SAMPLE_ID"], min.cells = seurat.createseuratobject.min.cells)
      }
      scrna_data[["Percent_Mito"]] <- PercentageFeatureSet(scrna_data, pattern = "^MT-")
      scrna_data <- subset(scrna_data, subset = nFeature_RNA > seurat.scrnaseq.subset.min.nFeature_RNA & nFeature_RNA < seurat.scrnaseq.subset.max.nFeature_RNA & Percent_Mito < seurat.scrnaseq.subset.max.percent.mt)
    scrna_data <- SCTransform(scrna_data,
                              assay = seurat.scrnaseq.sctransform.assay,
                              new.assay.name = seurat.scrnaseq.sctransform.new.assay.name,
                              reference.SCT.model = seurat.scrnaseq.sctransform.reference.SCT.model,
                              do.correct.umi = seurat.scrnaseq.sctransform.do.correct.umi,
                              ncells = seurat.scrnaseq.sctransform.ncells,
                              residual.features = seurat.scrnaseq.sctransform.residual.features,
                              variable.features.n = seurat.scrnaseq.sctransform.variable.features.n,
                              variable.features.rv.th = seurat.scrnaseq.sctransform.variable.features.rv.th,
                              vars.to.regress = seurat.scrnaseq.sctransform.vars.to.regress,
                              do.scale = seurat.scrnaseq.sctransform.do.scale,
                              do.center = seurat.scrnaseq.sctransform.do.center,
                              # clip.range = seurat.scrnaseq.sctransform.clip.range,
                              conserve.memory = seurat.scrnaseq.sctransform.conserve.memory,
                              return.only.var.genes = seurat.scrnaseq.sctransform.return.only.var.genes,
                              seed.use = seurat.scrnaseq.sctransform.seed.use,
                              verbose = FALSE) %>% RunPCA(assay = seurat.runpca.assay,
                                                          features = seurat.runpca.features,
                                                          npcs = seurat.runpca.npcs,
                                                          rev.pca = seurat.runpca.rev.pca,
                                                          weight.by.var = seurat.runpca.weight.by.var,
                                                          ndims.print = seurat.runpca.ndims.print,
                                                          nfeatures.print = seurat.runpca.nfeatures.print,
                                                          reduction.name = seurat.runpca.reduction.name,
                                                          reduction.key = seurat.runpca.reduction.key,
                                                          seed.use = seurat.runpca.seed.use,
                                                          verbose = FALSE)
    
    selected_pcs <- NULL
    selected_pcs <- find_selected_pcs(scrna_data,
                                      pcs = seurat.scrnaseq.runumap.pcs,
                                      pca = "pca",
                                      cname = pheno_data[i,"SAMPLE_ID"])
    scrna_data <- RunUMAP(scrna_data,
                          reduction = "pca",
                          dims = 1:selected_pcs,
                          features = seurat.scrnaseq.runumap.features,
                          graph = seurat.scrnaseq.runumap.graph,
                          assay = seurat.scrnaseq.runumap.assay,
    nn.name = seurat.scrnaseq.runumap.nn.name,
    slot = seurat.scrnaseq.runumap.slot,
    umap.method = seurat.scrnaseq.runumap.umap.method,
    reduction.model = seurat.scrnaseq.runumap.reduction.model,
    return.model = seurat.scrnaseq.runumap.return.model,
    n.neighbors = seurat.scrnaseq.runumap.n.neighbors,
    n.components = seurat.scrnaseq.runumap.n.components,
    metric = seurat.scrnaseq.runumap.metric,
    n.epochs = seurat.scrnaseq.runumap.n.epochs,
    learning.rate = seurat.scrnaseq.runumap.learning.rate,
    min.dist = seurat.scrnaseq.runumap.min.dist,
    spread = seurat.scrnaseq.runumap.spread,
    set.op.mix.ratio = seurat.scrnaseq.runumap.set.op.mix.ratio,
    local.connectivity = seurat.scrnaseq.runumap.local.connectivity,
    repulsion.strength = seurat.scrnaseq.runumap.repulsion.strength,
    negative.sample.rate = seurat.scrnaseq.runumap.negative.sample.rate,
    a = seurat.scrnaseq.runumap.a,
    b = seurat.scrnaseq.runumap.b,
    uwot.sgd = seurat.scrnaseq.runumap.uwot.sgd,
    seed.use = seurat.scrnaseq.runumap.seed.use,
    metric.kwds = seurat.scrnaseq.runumap.metric.kwds,
    angular.rp.forest = seurat.scrnaseq.runumap.angular.rp.forest,
    densmap = seurat.scrnaseq.runumap.densmap,
    dens.lambda = seurat.scrnaseq.runumap.dens.lambda,
    dens.frac = seurat.scrnaseq.runumap.dens.frac,
    dens.var.shift = seurat.scrnaseq.runumap.dens.var.shift,
    reduction.name = seurat.scrnaseq.runumap.reduction.name,
    reduction.key = seurat.scrnaseq.runumap.reduction.key)
    
    scrna_data <- RunTSNE(scrna_data,
                          reduction = "pca",
                          dims = 1:selected_pcs,
                          cells = seurat.scrnaseq.runtsne.cells,
                          features = seurat.scrnaseq.runtsne.features,
                          seed.use = seurat.scrnaseq.runtsne.seed.use,
                          tsne.method = seurat.scrnaseq.runtsne.tsne.method,
                          dim.embed = seurat.scrnaseq.runtsne.dim.embed,
                          distance.matrix = seurat.scrnaseq.runtsne.distance.matrix,
                          reduction.name = seurat.scrnaseq.runtsne.reduction.name,
                          reduction.key = seurat.scrnaseq.runtsne.reduction.key,
                          check_duplicates = seurat.scrnaseq.runtsne.check_duplicates)
    scrna_data <- FindNeighbors(scrna_data,
                                reduction = "pca",
                                dims = 1:selected_pcs,
                                assay = seurat.scrnaseq.findneighbors.assay,
                                features = seurat.scrnaseq.findneighbors.features,
                                k.param = seurat.scrnaseq.findneighbors.k.param,
                                return.neighbor = seurat.scrnaseq.findneighbors.return.neighbor,
                                # compute.SNN = seurat.scrnaseq.findneighbors.compute.SNN,
                                prune.SNN = seurat.scrnaseq.findneighbors.prune.SNN,
                                nn.method = seurat.scrnaseq.findneighbors.nn.method,
                                n.trees = seurat.scrnaseq.findneighbors.n.trees,
                                annoy.metric = seurat.scrnaseq.findneighbors.annoy.metric,
                                nn.eps = seurat.scrnaseq.findneighbors.nn.eps,
                                force.recalc = seurat.scrnaseq.findneighbors.force.recalc,
                                do.plot = seurat.scrnaseq.findneighbors.do.plot,
                                graph.name = seurat.scrnaseq.findneighbors.graph.name,
                                l2.norm = seurat.scrnaseq.findneighbors.l2.norm,
                                cache.index = seurat.scrnaseq.findneighbors.cache.index)
    scrna_data <- FindClusters(scrna_data,
                               graph.name = seurat.scrnaseq.findclusters.graph.name,
                               modularity.fxn = seurat.scrnaseq.findclusters.modularity.fxn,
                               initial.membership = seurat.scrnaseq.findclusters.initial.membership,
                               node.sizes = seurat.scrnaseq.findclusters.node.sizes,
                               resolution = seurat.scrnaseq.findclusters.resolution,
                               method = seurat.scrnaseq.findclusters.method,
                               algorithm = seurat.scrnaseq.findclusters.algorithm,
                               n.start = seurat.scrnaseq.findclusters.n.start,
                               n.iter = seurat.scrnaseq.findclusters.n.iter,
                               random.seed = seurat.scrnaseq.findclusters.random.seed,
                               group.singletons = seurat.scrnaseq.findclusters.group.singletons,
                               temp.file.location = seurat.scrnaseq.findclusters.temp.file.location,
                               edge.file.name = seurat.scrnaseq.findclusters.edge.file.name,
                               verbose = FALSE)

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
                       ref = singler.scrnaseq.ref,
                       labels = singler.scrnaseq.labels,
                       method = singler.scrnaseq.method,
                       genes = singler.scrnaseq.genes,
                       sd.thresh = singler.scrnaseq.sd.thresh,
                       de.method = singler.scrnaseq.de.method,
                       de.n = singler.scrnaseq.de.n,
                       de.args = singler.scrnaseq.de.args,
                       aggr.ref = singler.scrnaseq.aggr.ref,
                       aggr.args = singler.scrnaseq.aggr.args,
                       recompute = singler.scrnaseq.recompute,
                       restrict = singler.scrnaseq.restrict,
                       quantile = singler.scrnaseq.quantile,
                       fine.tune = singler.scrnaseq.fine.tune,
                       tune.thresh = singler.scrnaseq.tune.thresh,
                       prune = singler.scrnaseq.prune,
                       assay.type.test = singler.scrnaseq.assay.type.test,
                       assay.type.ref = singler.scrnaseq.assay.type.ref,
                       check.missing = singler.scrnaseq.check.missing)
                       # num.threads = singler.scrnaseq.num.threads
                       # BNPARAM = singler.scrnaseq.BNPARAM
                       # BPPARAM = singler.scrnaseq.BPPARAM

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

    data_anchors <- FindTransferAnchors(reference = scrna_data,
                                        query = current,
                                        normalization.method = seurat.multimodalintegration.findfransferanchors.normalization.method, # SCT
                                        recompute.residuals = seurat.multimodalintegration.findfransferanchors.recompute.residuals,
                                        reference.assay = seurat.multimodalintegration.findfransferanchors.reference.assay,
                                        reference.neighbors = seurat.multimodalintegration.findfransferanchors.reference.neighbors,
                                        query.assay = seurat.multimodalintegration.findfransferanchors.query.assay,
                                        reduction = seurat.multimodalintegration.findfransferanchors.reduction,
                                        reference.reduction = seurat.multimodalintegration.findfransferanchors.reference.reduction,
                                        project.query = seurat.multimodalintegration.findfransferanchors.project.query,
                                        features = seurat.multimodalintegration.findfransferanchors.features,
                                        scale = seurat.multimodalintegration.findfransferanchors.scale,
                                        npcs = seurat.multimodalintegration.findfransferanchors.npcs,
                                        l2.norm = seurat.multimodalintegration.findfransferanchors.l2.norm,
                                        dims = 1:selected_pcs,
                                        k.anchor = seurat.multimodalintegration.findfransferanchors.k.anchor,
                                        k.filter = seurat.multimodalintegration.findfransferanchors.k.filter,
                                        k.score = seurat.multimodalintegration.findfransferanchors.k.score,
                                        max.features = seurat.multimodalintegration.findfransferanchors.max.features,
                                        nn.method = seurat.multimodalintegration.findfransferanchors.nn.method,
                                        n.trees = seurat.multimodalintegration.findfransferanchors.n.trees,
                                        eps = seurat.multimodalintegration.findfransferanchors.eps,
                                        approx.pca = seurat.multimodalintegration.findfransferanchors.approx.pca,
                                        mapping.score.k = seurat.multimodalintegration.findfransferanchors.mapping.score.k)
    
    predictions_assay <- TransferData(anchorset = data_anchors,
                                      refdata = scrna_data$Cell_Type,
                                      dims = 1:selected_pcs,
                                      reference = seurat.multimodalintegration.transferdata.reference,
                                      query = seurat.multimodalintegration.transferdata.query,
                                      weight.reduction = seurat.multimodalintegration.transferdata.weight.reduction,
                                      l2.norm = seurat.multimodalintegration.transferdata.l2.norm,
                                      k.weight = seurat.multimodalintegration.transferdata.k.weight,
                                      sd.weight = seurat.multimodalintegration.transferdata.sd.weight,
                                      eps = seurat.multimodalintegration.transferdata.eps,
                                      n.trees = seurat.multimodalintegration.transferdata.n.trees,
                                      slot = seurat.multimodalintegration.transferdata.slot,
                                      prediction.assay = seurat.multimodalintegration.transferdata.prediction.assay,
                                      store.weights = seurat.multimodalintegration.transferdata.store.weights)
    
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
      saveRDS(scrna_data, paste(cdir,"18URSA_SPATIAL_SCRNASEQ_DATA_",gsub("\\/","",pheno_data[i,"SAMPLE_ID"]),"_",project_name,".RDS", sep = ""))
      saveRDS(current, paste(cdir,"19URSA_SPATIAL_SAMPLE_DATA_WITH_SCRNASEQ_TRANSFERED_LABELS_",gsub("\\/","",pheno_data[i,"SAMPLE_ID"]),"_",project_name,".RDS", sep = ""))
    }
    }
saveRDS(data_current, paste(cdir,"20URSA_SPATIAL_SAMPLE_DATA_",project_name,".RDS", sep = ""))

print("Completed!")
}
