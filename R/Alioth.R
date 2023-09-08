############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Alioth: scRNASEQ
# Version: V1.2.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Last Update Date: 2023-01-29
############################################################################################################
#' @include ini.R
#' @include common.R
#' @import ALRA
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
#' @import metap
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
#' @param seurat.subset.min.nFeature_RNA Lower bound to select for cells with
#' number of genes more than this value. Default is 200.
#' @param seurat.subset.max.nFeature_RNA Upper bound to select for cells with
#' number of genes less than or equals to this value Default is 25000.
#' @param seurat.max.mito.percent Threshold percentage to select for cells with
#' percentage of mitochondria genes less than this value.
#' This is remove cells which are debris or artefacts and doublets.
#' Default is 10.
#' @param integration_method Integration method.
#' Accepts 'Seurat' or "Harmony' integration. Default is Seurat.
#' @param cc_regression Cell cycle regression method.
#' Accepts 0 for no regression,
#' 1 for regression with two-phases, 2 for regression with phase
#' difference (See Seurat). Default is 0.
#' @param seurat.runumap.pcs (1). Number of principal components (PCs) or Harmonys
#' (if harmony was chosen for integration) to be used for the analysis. Default to 30.
#' Please state a reasonable number which is no more than the total number
#' of cells submitted to avoid running into errors. Similarly for integrated data,
#' the corresponding parameter is: seurat.integration.runumap.pcs;
#' (2). alternatively, provide a method name to run an automatic selection of the
#' number of PCs or Harmonys. Methods include 'none', 'all', 'piecewise linear model',
#' 'first derivative', 'second derivative', 'preceding residual', 'perpendicular line',
#' and 'k-means clustering'. If 'all' is used, all methods will be assessed and the
#' final PC number will be determined via the mean of the output of all
#' methods. For more information, please visit:
#' https://github.com/haotian-zhuang/findPC
#' (Zhuang et al., Bioinformatics, 2022).
#' Default value for seurat.runumap.pcs and seurat.integration.runumap.pcs
#' is set to 30.
#' @param to_impute Pass 'YES' to perform imputation on individual
#' samples or 'NO' to skip imputation process. Default is set to 'NO'.
#' @param find_doublet Pass 'YES' to perform doublets removal using
#'  DoubletFinder (Christopher S. M., Cell Systems, 2019). Default
#' is set to 'YES'. If set to 'NO', potential doublets will not be removed.
#' @param run_unbias_vis Plot additional unbias visualizations for top
#' features using SCUBI (Wenpin H. & Zhicheng J., Cell Reports Methods,
#' 2021) for dimensionally reduced plots such as UMAP. Pass 'YES' to
#' plot, default to 'NO'.
#' @param ... Arguments passed to other parameters in the dependency pages.
#' Parameters with the long format: xxx.xxx.xxx usually indicates in lowercases
#' the parameter origin: <dependent package name>.<function name>.<parameter name>(),
#' for example: seurat.load10x_spatial.assay() indicates this parameter 
#' originates from the package Seurat under the function Load10X_Spatial() and its
#' parameter 'assay'. Users could subsequently refer to the dependency 
#' package for detailed information on parameters and their usage.
#' @export

scRNASEQPip <- function(project_name = "Ursa_scRNASEQ",
                      input_dir = "./",
                      output_dir = "./",
                      pheno_file = "study_meta.csv",
                      integration_method = "Seurat",
                      cc_regression = 0,
                      to_impute = "NO",
                      find_doublet = "YES",
                      run_unbias_vis = "NO",
                      
                      # CreateSeuratObject
                      seurat.subset.min.cells = 3,
                      seurat.subset.min.features = 200,
                      
                      # Subset
                      seurat.subset.min.nFeature_RNA = 200,
                      seurat.subset.max.nFeature_RNA = 25000,
                      seurat.max.mito.percent = 10,
                      
                      # NormalizeData
                      seurat.normalizedata.normalization.method = "LogNormalize",
                      seurat.normalizedata.norm.scale.factor = 10000,
                      seurat.normalizedata.norm.margin = 1,
                      
                      # FindVariableFeatures
                      seurat.findvariablefeatures.nfeatures = 2000,
                      seurat.findvariablefeatures.selection.method = "vst",
                      seurat.findvariablefeatures.loess.span = 0.3,
                      seurat.findvariablefeatures.clip.max = "auto",
                      seurat.findvariablefeatures.mean.function = "FastExpMean",
                      seurat.findvariablefeatures.dispersion.function = "FastLogVMR",
                      seurat.findvariablefeatures.num.bin = 20,
                      seurat.findvariablefeatures.binning.method = "equal_width",
                      
                      # ScaleData
                      seurat.scaledata.features = NULL,
                      seurat.scaledata.vars.to.regress = NULL,
                      # seurat.scaledata.latent.data = NULL,
                      seurat.scaledata.split.by = NULL,
                      seurat.scaledata.model.use = "linear",
                      seurat.scaledata.use.umi = FALSE,
                      seurat.scaledata.do.scale = TRUE,
                      seurat.scaledata.do.center = TRUE,
                      seurat.scaledata.scale.max = 10,
                      seurat.scaledata.block.size = 1000,
                      seurat.scaledata.min.cells.to.block = 3000,
                      
                      # CellCycleScoring
                      seurat.cellcyclescoring.ctrl = NULL,
                      seurat.cellcyclescoring.set.ident = FALSE,
                      
                      # CellCycleScoring - RunPCA
                      seurat.cellcyclescoring.seurat.runpca.npcs = 50,
                      seurat.cellcyclescoring.seurat.runpca.rev.pca = FALSE,
                      seurat.cellcyclescoring.seurat.runpca.weight.by.var = TRUE,
                      seurat.cellcyclescoring.seurat.runpca.ndims.print = 5,
                      seurat.cellcyclescoring.seurat.runpca.nfeatures.print = 30,
                      seurat.cellcyclescoring.seurat.runpca.reduction.key = "PC_",
                      seurat.cellcyclescoring.seurat.runpca.seed.use = 42,
                      seurat.cellcyclescoring.seurat.runpca.approx = TRUE,
                      
                      # No CellCycleScoring - RunPCA
                      seurat.runpca.features = NULL,
                      seurat.runpca.npcs = 50,
                      seurat.runpca.rev.pca = FALSE,
                      seurat.runpca.weight.by.var = TRUE,
                      seurat.runpca.ndims.print = 5,
                      seurat.runpca.nfeatures.print = 30,
                      seurat.runpca.reduction.key = "PC_",
                      seurat.runpca.seed.use = 42,
                      seurat.runpca.approx = TRUE,
                      
                      # RunUMAP
                      seurat.runumap.pcs = 30,
                      seurat.runumap.reduction = "pca",
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
                      
                      # Choose k for imputation
                      alra.choosek.k = 100,
                      alra.choosek.thresh = 6,
                      alra.choosek.noise_start = 80,
                      alra.choosek.q = 2,
                      alra.choosek.use.mkl = F,
                      alra.choosek.mkl.seed = -1,
                      
                      # Imputation - alra
                      alra.alra.q = 10,
                      alra.alra.quantile.prob = 0.001,
                      alra.alra.use.mkl = F,
                      alra.alra.mkl.seed = -1,
                      
                      # DoubletFinder - paramSweep_v3
                      doubletfinder.paramSweep_v3.PCs = 30,
                      doubletfinder.paramSweep_v3.sct = FALSE,
                      
                      # DoubletFinder - summarizeSweep
                      doubletfinder.summarizeSweep.gt = FALSE,
                      doubletfinder.summarizeSweep.classifications = NULL,
                      
                      # DoubletFinder - doubletFinder_v3
                      doubletfinder.homotypic.prop = 0.075,
                      doubletfinder.doubletFinder_v3.PCs = 30,
                      doubletfinder.doubletFinder_v3.pN = 0.25,
                      doubletfinder.doubletFinder_v3.reuse.pANN = FALSE,
                      doubletfinder.doubletFinder_v3.sct = FALSE,
                      doubletfinder.doubletFinder_v3.annotations = NULL,
                      
                      # FindNeighbors
                      seurat.findneighbors.reduction = "pca",
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
                      
                      # RunTSNE
                      seurat.runtsne.reduction = "pca",
                      seurat.runtsne.cells = NULL,
                      seurat.runtsne.features = NULL,
                      seurat.runtsne.seed.use = 1,
                      seurat.runtsne.tsne.method = "Rtsne",
                      seurat.runtsne.dim.embed = 2,
                      seurat.runtsne.distance.matrix = NULL,
                      seurat.runtsne.reduction.name = "tsne",
                      seurat.runtsne.reduction.key = "tSNE_",
                      seurat.runtsne.check_duplicates = FALSE,
                      
                      # In deanalysis function - FindAllMarkers
                      deanalysis.p_threshold = 0.05,
                      seurat.findallmarkers.min.pct = 0.5,
                      seurat.findallmarkers.logfc.threshold = 0.25,
                      seurat.findallmarkers.assay = NULL,
                      seurat.findallmarkers.features = NULL,
                      seurat.findallmarkers.test.use = "wilcox",
                      seurat.findallmarkers.slot = "data",
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
                      seurat.findallmarkers.return.thresh = 0.01,
                      seurat.findallmarkers.densify = FALSE,
                      
                      # In deanalysis function - FindConservedMarkers
                      seurat.findconservedmarkers.ident.2 = NULL,
                      seurat.findconservedmarkers.assay = "RNA",
                      seurat.findconservedmarkers.slot = "data",
                      seurat.findconservedmarkers.min.cells.group = 3,
                      seurat.findconservedmarkers.meta.method = metap::minimump,
                      
                      # In deanalysis function - FindMarkers
                      seurat.findmarkers.group.by = NULL,
                      seurat.findmarkers.subset.ident = NULL,
                      seurat.findmarkers.assay = NULL,
                      seurat.findmarkers.slot = "data",
                      seurat.findmarkers.reduction = NULL,
                      seurat.findmarkers.features = NULL,
                      seurat.findmarkers.logfc.threshold = 0.25,
                      seurat.findmarkers.test.use = "wilcox",
                      seurat.findmarkers.min.pct = 0.1,
                      seurat.findmarkers.min.diff.pct = -Inf,
                      seurat.findmarkers.only.pos = FALSE,
                      seurat.findmarkers.max.cells.per.ident = Inf,
                      seurat.findmarkers.random.seed = 1,
                      seurat.findmarkers.latent.vars = NULL,
                      seurat.findmarkers.min.cells.feature = 3,
                      seurat.findmarkers.min.cells.group = 3,
                      seurat.findmarkers.mean.fxn = NULL,
                      seurat.findmarkers.fc.name = NULL,
                      seurat.findmarkers.base = 2,
                      seurat.findmarkers.densify = FALSE,
                      
                      # scubi
                      # scubi.resolution = 10 * 20/max(diff(range(dim1)), diff(range(dim2))),
                      scubi.palette = rainbow(15)[c(11:1,15)],
                      
                      # SingleR
                      singler.ref = "HumanPrimaryCellAtlas",
                      singler.labels.main = "HumanPrimaryCellAtlasLevel1",
                      singler.labels.fine = "HumanPrimaryCellAtlasLevel2",
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
                      
                      # DOSE - enrichDGN
                      dose.enrichdgn.pvalueCutoff = 0.05,
                      dose.enrichdgn.pAdjustMethod = "BH",
                      dose.enrichdgn.universe = NULL,
                      dose.enrichdgn.minGSSize = 10,
                      dose.enrichdgn.maxGSSize = 500,
                      dose.enrichdgn.qvalueCutoff = 0.05,
                      dose.enrichdgn.readable = TRUE,
                      
                      # enrichplot - pairwise_termsim
                      enrichplot.pairwise_termsim.method = "JC",
                      enrichplot.pairwise_termsim.semData = NULL,
                      
                      # clusterProfiler - compareCluster
                      clusterprofiler.comparecluster.pvaluecutoff = 0.05,
                      
                      # DOSE - gseNCG
                      dose.gsencg.exponent = 1,
                      dose.gsencg.minGSSize = 10,
                      dose.gsencg.maxGSSize = 500,
                      dose.gsencg.pvalueCutoff = 0.05,
                      dose.gsencg.pAdjustMethod = "BH",
                      dose.gsencg.seed = FALSE,
                      dose.gsencg.by = "fgsea",
                      
                      # Seurat - SelectIntegrationFeatures
                      seurat.selectintegrationfeatures.nfeatures = 2000,
                      seurat.selectintegrationfeatures.assay = NULL,
                      seurat.selectintegrationfeatures.fvf.nfeatures = 2000,
                      
                      # Seurat - FindIntegrationAnchors
                      seurat.findintegrationanchors.assay = NULL,
                      seurat.findintegrationanchors.reference = NULL,
                      seurat.findintegrationanchors.anchor.features = 2000,
                      seurat.findintegrationanchors.scale = TRUE,
                      seurat.findintegrationanchors.normalization.method = c("LogNormalize", "SCT"),
                      seurat.findintegrationanchors.sct.clip.range = NULL,
                      seurat.findintegrationanchors.reduction = "rpca",
                      seurat.findintegrationanchors.l2.norm = TRUE,
                      seurat.findintegrationanchors.dims = 1:30,
                      seurat.findintegrationanchors.k.anchor = 5,
                      seurat.findintegrationanchors.k.filter = 200,
                      seurat.findintegrationanchors.k.score = 30,
                      seurat.findintegrationanchors.max.features = 200,
                      seurat.findintegrationanchors.nn.method = "annoy",
                      seurat.findintegrationanchors.n.trees = 50,
                      seurat.findintegrationanchors.eps = 0,
                      
                      # Seurat - IntegrateData
                      seurat.integratedata.dims = 30,
                      seurat.integratedata.new.assay.name = "integrated",
                      seurat.integratedata.normalization.method = c("LogNormalize", "SCT"),
                      seurat.integratedata.features = NULL,
                      seurat.integratedata.features.to.integrate = NULL,
                      seurat.integratedata.k.weight = 100,
                      seurat.integratedata.weight.reduction = NULL,
                      seurat.integratedata.sd.weight = 1,
                      seurat.integratedata.sample.tree = NULL,
                      seurat.integratedata.preserve.order = FALSE,
                      seurat.integratedata.eps = 0,
                      
                      # Seurat - Integration ScaleData
                      seurat.integration.scaledata.features = NULL,
                      seurat.integration.scaledata.vars.to.regress = NULL,
                      # seurat.integration.scaledata.latent.data = NULL,
                      seurat.integration.scaledata.split.by = NULL,
                      seurat.integration.scaledata.model.use = "linear",
                      seurat.integration.scaledata.use.umi = FALSE,
                      seurat.integration.scaledata.do.scale = TRUE,
                      seurat.integration.scaledata.do.center = TRUE,
                      seurat.integration.scaledata.scale.max = 10,
                      seurat.integration.scaledata.block.size = 1000,
                      seurat.integration.scaledata.min.cells.to.block = 3000,
                      
                      # Seurat - Integration FindVariableFeatures
                      seurat.integration.findvariablefeatures.nfeatures = 2000,
                      seurat.integration.findvariablefeatures.selection.method = "vst",
                      seurat.integration.findvariablefeatures.loess.span = 0.3,
                      seurat.integration.findvariablefeatures.clip.max = "auto",
                      seurat.integration.findvariablefeatures.mean.function = "FastExpMean",
                      seurat.integration.findvariablefeatures.dispersion.function = "FastLogVMR",
                      seurat.integration.findvariablefeatures.num.bin = 20,
                      seurat.integration.findvariablefeatures.binning.method = "equal_width",
                      
                      # Seurat - Integration RunPCA
                      seurat.integration.runpca.features = NULL,
                      seurat.integration.runpca.npcs = 50,
                      seurat.integration.runpca.rev.pca = FALSE,
                      seurat.integration.runpca.weight.by.var = TRUE,
                      seurat.integration.runpca.ndims.print = 1:5,
                      seurat.integration.runpca.nfeatures.print = 30,
                      seurat.integration.runpca.reduction.key = "PC_",
                      seurat.integration.runpca.seed.use = 42,
                      seurat.integration.runpca.approx = TRUE,
                      
                      # Seurat - Integration RunUMAP
                      seurat.integration.runumap.pcs = 30,
                      seurat.integration.runumap.features = NULL,
                      seurat.integration.runumap.graph = NULL,
                      seurat.integration.runumap.assay = NULL,
                      seurat.integration.runumap.nn.name = NULL,
                      seurat.integration.runumap.slot = "data",
                      seurat.integration.runumap.umap.method = "uwot",
                      seurat.integration.runumap.reduction.model = NULL,
                      seurat.integration.runumap.return.model = FALSE,
                      seurat.integration.runumap.n.neighbors = 30L,
                      seurat.integration.runumap.n.components = 2L,
                      seurat.integration.runumap.metric = "cosine",
                      seurat.integration.runumap.n.epochs = NULL,
                      seurat.integration.runumap.learning.rate = 1,
                      seurat.integration.runumap.min.dist = 0.3,
                      seurat.integration.runumap.spread = 1,
                      seurat.integration.runumap.set.op.mix.ratio = 1,
                      seurat.integration.runumap.local.connectivity = 1L,
                      seurat.integration.runumap.repulsion.strength = 1,
                      seurat.integration.runumap.negative.sample.rate = 5L,
                      seurat.integration.runumap.a = NULL,
                      seurat.integration.runumap.b = NULL,
                      seurat.integration.runumap.uwot.sgd = FALSE,
                      seurat.integration.runumap.seed.use = 42L,
                      seurat.integration.runumap.metric.kwds = NULL,
                      seurat.integration.runumap.angular.rp.forest = FALSE,
                      seurat.integration.runumap.densmap = FALSE,
                      seurat.integration.runumap.dens.lambda = 2,
                      seurat.integration.runumap.dens.frac = 0.3,
                      seurat.integration.runumap.dens.var.shift = 0.1,
                      seurat.integration.runumap.reduction.name = "umap",
                      seurat.integration.runumap.reduction.key = "UMAP_",
                      
                      # Seurat - Integration RunTSNE
                      seurat.integration.runtsne.cells = NULL,
                      seurat.integration.runtsne.features = NULL,
                      seurat.integration.runtsne.seed.use = 1,
                      seurat.integration.runtsne.tsne.method = "Rtsne",
                      seurat.integration.runtsne.dim.embed = 2,
                      seurat.integration.runtsne.distance.matrix = NULL,
                      seurat.integration.runtsne.reduction.name = "tsne",
                      seurat.integration.runtsne.reduction.key = "tSNE_",
                      seurat.integration.runtsne.check_duplicates = FALSE,
                      
                      # Harmony - RunHarmony
                      harmony.runharmony.reduction = "pca",
                      harmony.runharmony.group.by.vars = "BATCH",
                      harmony.runharmony.theta = NULL,
                      harmony.runharmony.lambda = NULL,
                      harmony.runharmony.sigma = 0.1,
                      harmony.runharmony.nclust = NULL,
                      harmony.runharmony.tau = 0,
                      harmony.runharmony.block.size = 0.05,
                      harmony.runharmony.max.iter.harmony = 10,
                      harmony.runharmony.max.iter.cluster = 20,
                      harmony.runharmony.epsilon.cluster = 1e-05,
                      harmony.runharmony.epsilon.harmony = 1e-04,
                      harmony.runharmony.plot_convergence = FALSE,
                      harmony.runharmony.reference_values = NULL,
                      harmony.runharmony.reduction.save = "harmony",
                      harmony.runharmony.assay.use = NULL,
                      harmony.runharmony.project.dim = TRUE,

                      # Seurat - Inetgration FindNeighbors
                      seurat.integration.findneighbors.assay = NULL,
                      seurat.integration.findneighbors.features = NULL,
                      seurat.integration.findneighbors.k.param = 20,
                      seurat.integration.findneighbors.return.neighbor = FALSE,
                      # seurat.integration.findneighbors.compute.SNN = !return.neighbor,
                      seurat.integration.findneighbors.prune.SNN = 1/15,
                      seurat.integration.findneighbors.nn.method = "annoy",
                      seurat.integration.findneighbors.n.trees = 50,
                      seurat.integration.findneighbors.annoy.metric = "euclidean",
                      seurat.integration.findneighbors.nn.eps = 0,
                      seurat.integration.findneighbors.force.recalc = FALSE,
                      seurat.integration.findneighbors.do.plot = FALSE,
                      seurat.integration.findneighbors.graph.name = NULL,
                      seurat.integration.findneighbors.l2.norm = FALSE,
                      seurat.integration.findneighbors.cache.index = FALSE,
                      
                      # Seurat - Integration FindClusters
                      seurat.integration.findclusters.graph.name = NULL,
                      seurat.integration.findclusters.modularity.fxn = 1,
                      seurat.integration.findclusters.initial.membership = NULL,
                      seurat.integration.findclusters.node.sizes = NULL,
                      seurat.integration.findclusters.resolution = 0.8,
                      seurat.integration.findclusters.method = "matrix",
                      seurat.integration.findclusters.algorithm = 1,
                      seurat.integration.findclusters.n.start = 10,
                      seurat.integration.findclusters.n.iter = 10,
                      seurat.integration.findclusters.random.seed = 0,
                      seurat.integration.findclusters.group.singletons = TRUE,
                      seurat.integration.findclusters.temp.file.location = NULL,
                      seurat.integration.findclusters.edge.file.name = NULL,
                      
                      # Celltalker - celltalk
                      metadata_grouping = "INTEGRATED_CELL_TYPE",
                      celltalker.celltalk.ligand_receptor_pairs = ramilowski_pairs,
                      celltalker.celltalk.number_cells_required = 10,
                      celltalker.celltalk.min_expression = 10,
                      celltalker.celltalk.max_expression = 20000,
                      celltalker.celltalk.scramble_times = 100){
  
  print("Initialising pipeline environment..")
  hs <- org.Hs.eg.db
  data("hgnc.table", package="HGNChelper")
  head(hgnc.table)
  hpca.se <- HumanPrimaryCellAtlasData()
  if(toupper(singler.ref) == toupper("HumanPrimaryCellAtlas")){
    singler.ref <- hpca.se
    if(toupper(singler.labels.main) == toupper("HumanPrimaryCellAtlasLevel1")){
      singler.labels.main <-hpca.se$label.main
    }
    if(toupper(singler.labels.fine) == toupper("HumanPrimaryCellAtlasLevel2")){
      singler.labels.fine <- hpca.se$label.fine
    }
  }
  cc.s.genes <- cc.genes$s.genes
  cc.g2m.genes <- cc.genes$g2m.genes
  pheno_data <- pheno_ini(pheno_file, pipeline = "scRNASEQ", isDir = T)
  pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, sep = "_")
  color_conditions <- color_ini()
  ctime <- time_ini()
  sample_colors <- gen_colors(n = length(unique(pheno_data$SID)))
  names(sample_colors) <- unique(pheno_data$SID)
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir, recursive = T)

  sample_files <- list.files(input_dir, pattern = ".*h5", full.names = T, ignore.case = T, recursive = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]

  # Imputation
  # source("https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")

  #######################################################################################################################################
  print("Preparing..")
  data <- NULL
  data_current <- NULL
  annot_names <- NULL
  results <- NULL

  for(j in 1:nrow(pheno_data)){
    print(paste("Running ", pheno_data[j,"FILE"], "..", sep = ""))
    annot_names <- c(annot_names, pheno_data[j,"SID"])
    data_current[[j]] <- NULL
    data_current[[j]] <- sample_files[grep(pheno_data[j,"FILE"], sample_files, ignore.case = T)]
    data_current[[j]] <- Read10X_h5(data_current[[j]])
    if((length(data_current[[j]]) == 1) | (length(data_current[[j]]) > 10)){
      data_current[[j]] <- CreateSeuratObject(counts = data_current[[j]], project = annot_names[j], min.cells = seurat.subset.min.cells, min.features = seurat.subset.min.features)
    }else{
      data_current[[j]] <- CreateSeuratObject(counts = data_current[[j]]$`Gene Expression`, project = annot_names[j], min.cells = seurat.subset.min.cells, min.features = seurat.subset.min.features)
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

    data_current[[j]] <- NormalizeData(data_current[[j]],
                                       normalization.method = seurat.normalizedata.normalization.method,
                                       scale.factor = seurat.normalizedata.norm.scale.factor,
                                       margin = seurat.normalizedata.norm.margin,
                                       verbose = TRUE)
    
    data_current[[j]]@assays$RNA@data@x[is.na(data_current[[j]]@assays$RNA@data@x)] <- 0
    print("FindVariableFeatures..")
    data_current[[j]] <- FindVariableFeatures(data_current[[j]],
                                              nfeatures = seurat.findvariablefeatures.nfeatures,
                                              selection.method = seurat.findvariablefeatures.selection.method,
                                              loess.span = seurat.findvariablefeatures.loess.span,
                                              clip.max = seurat.findvariablefeatures.clip.max,
                                              mean.function = seurat.findvariablefeatures.mean.function,
                                              dispersion.function = seurat.findvariablefeatures.dispersion.function,
                                              num.bin = seurat.findvariablefeatures.num.bin,
                                              binning.method = seurat.findvariablefeatures.binning.method,
                                              verbose = TRUE)
    data_current[[j]] <- ScaleData(data_current[[j]],
                                   features = seurat.scaledata.features,
                                   vars.to.regress = seurat.scaledata.vars.to.regress,
                                   # latent.data = seurat.scaledata.latent.data,
                                   split.by = seurat.scaledata.split.by,
                                   model.use = seurat.scaledata.model.use,
                                   use.umi = seurat.scaledata.use.umi,
                                   do.scale = seurat.scaledata.do.scale,
                                   do.center = seurat.scaledata.do.center,
                                   scale.max = seurat.scaledata.scale.max,
                                   block.size = seurat.scaledata.block.size,
                                   min.cells.to.block = seurat.scaledata.min.cells.to.block,
                                   verbose = TRUE)

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

      data_current[[j]] <- CellCycleScoring(data_current[[j]],
                                            s.features = cc.s.genes,
                                            g2m.features = cc.g2m.genes,
                                            ctrl = seurat.cellcyclescoring.ctrl,
                                            set.ident = seurat.cellcyclescoring.set.ident)
      data_current[[j]]@meta.data$Phase <- factor(data_current[[j]]@meta.data$Phase, levels = c("G1","S","G2M"))
      data_current[[j]] <- RunPCA(data_current[[j]],
                                  features = c(cc.s.genes, cc.g2m.genes),
                                  npcs = seurat.cellcyclescoring.seurat.runpca.npcs,
                                  rev.pca = seurat.cellcyclescoring.seurat.runpca.rev.pca,
                                  weight.by.var = seurat.cellcyclescoring.seurat.runpca.weight.by.var,
                                  ndims.print = 1:seurat.cellcyclescoring.seurat.runpca.ndims.print,
                                  nfeatures.print = seurat.cellcyclescoring.seurat.runpca.nfeatures.print,
                                  reduction.key = seurat.cellcyclescoring.seurat.runpca.reduction.key,
                                  seed.use = seurat.cellcyclescoring.seurat.runpca.seed.use,
                                  approx = seurat.cellcyclescoring.seurat.runpca.approx,
                                  verbose = TRUE)
      data_current[[j]]@reductions$pca_selected <- data_current[[j]]@reductions$pca

      if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
        row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
        row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
        row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
      }

      data_current[[j]] <- RunPCA(data_current[[j]],
                                  features = seurat.runpca.features,
                                  npcs = seurat.runpca.npcs,
                                  rev.pca = seurat.runpca.rev.pca,
                                  weight.by.var = seurat.runpca.weight.by.var,
                                  ndims.print = 1:seurat.runpca.ndims.print,
                                  nfeatures.print = seurat.runpca.nfeatures.print,
                                  reduction.key = seurat.runpca.reduction.key,
                                  seed.use = seurat.runpca.seed.use,
                                  approx = seurat.runpca.approx,
                                  verbose = TRUE)

      # G2/M and S Phase Markers: Tirosh et al, 2015
      p1 <- own_2d_scatter(data_current[[j]], "pca", "Phase", "PCA Based on Variable Features")
      p2 <- own_2d_scatter(data_current[[j]], "pca_selected", "Phase", "PCA Based on G2/M and S Phase Markers")

      somePDFPath = paste(cdir,"3URSA_PLOT_scRNASEQ_CELL_CYCLE_PHASES_PCA_BEFORE_CC_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=7,pointsize=12)
      print(p1 / p2 + plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

      selected_pcs <- NULL
      selected_pcs <- find_selected_pcs(data_current[[j]],
                                        pcs = seurat.runumap.pcs,
                                        pca = "pca_selected",
                                        cname = annot_names[j])
      
      data_current[[j]] <- RunUMAP(data_current[[j]],
                                   reduction = "pca_selected",
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
                                   reduction.key = seurat.runumap.reduction.key,
                                   verbose = TRUE)
      data_current[[j]]@reductions$umap_selected <- data_current[[j]]@reductions$umap
      data_current[[j]] <- RunUMAP(data_current[[j]],
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
                                   reduction.key = seurat.runumap.reduction.key,
                                   verbose = TRUE)

      p1 <- own_2d_scatter(data_current[[j]], "umap", "Phase", "UMAP Based on Variable Features")
      p2 <- own_2d_scatter(data_current[[j]], "umap_selected", "Phase", "UMAP Based on G2/M and S Phase Markers")

      somePDFPath = paste(cdir,"4URSA_PLOT_scRNASEQ_CELLCYCLE_UMAP_BEFORE_REGRESSION_",annot_names[j],"_",project_name,".pdf", sep = "")
      pdf(file=somePDFPath, width=10, height=7,pointsize=12)
      print(p1 / p2 + plot_annotation(title = paste("Before Cell Cycle Regression - ",annot_names[j], sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
      dev.off()

    }else{
      data_current[[j]] <- RunPCA(data_current[[j]],
                                  features = seurat.runpca.features,
                                  npcs = seurat.runpca.npcs,
                                  rev.pca = seurat.runpca.rev.pca,
                                  weight.by.var = seurat.runpca.weight.by.var,
                                  ndims.print = 1:seurat.runpca.ndims.print,
                                  nfeatures.print = seurat.runpca.nfeatures.print,
                                  reduction.key = seurat.runpca.reduction.key,
                                  seed.use = seurat.runpca.seed.use,
                                  approx = seurat.runpca.approx,
                                  verbose = TRUE)

      selected_pcs <- NULL
      selected_pcs <- find_selected_pcs(data_current[[j]],
                                        pcs = seurat.runumap.pcs,
                                        pca = "pca",
                                        cname = annot_names[j])
      
      data_current[[j]] <- RunUMAP(data_current[[j]],
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
                                   reduction.key = seurat.runumap.reduction.key,
                                   verbose = TRUE)
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
    data_current[[j]] <- subset(data_current[[j]], subset = nFeature_RNA > seurat.subset.min.nFeature_RNA & nFeature_RNA <= seurat.subset.max.nFeature_RNA & Percent_Mito < seurat.max.mito.percent)

    if(ncol(data_current[[j]]) > 200){

      data_current[[j]] <- NormalizeData(data_current[[j]],
                                         normalization.method = seurat.normalizedata.normalization.method,
                                         scale.factor = seurat.normalizedata.norm.scale.factor,
                                         margin = seurat.normalizedata.norm.margin,
                                         verbose = TRUE)

      # Imputation
      if(toupper(to_impute) == "YES"){
        print("Performing imputation using ALRA..")
        cimpute_norm <- NULL
        cimpute_norm <- t(as.matrix(data_current[[j]]@assays$RNA@data))
        k_choice <- choose_k(cimpute_norm,
                             K = alra.choosek.k,
                             thresh = alra.choosek.thresh,
                             noise_start = alra.choosek.noise_start,
                             q = alra.choosek.q,
                             use.mkl = alra.choosek.use.mkl,
                             mkl.seed = alra.choosek.mkl.seed)
        cimpute_norm <- alra(cimpute_norm,
                             k = k_choice$k,
                             q = alra.alra.q,
                             quantile.prob = alra.alra.quantile.prob,
                             use.mkl = alra.alra.use.mkl,
                             mkl.seed = alra.alra.mkl.seed)[[3]]
        cimpute_norm <- t(cimpute_norm)
        colnames(cimpute_norm) <- colnames(data_current[[j]])
        # data_current[[j]]@assays$RNA@data_original <- data_current[[j]]@assays$RNA@data
        data_current[[j]]@assays$RNA@data <- cimpute_norm
      }

      data_current[[j]]@assays$RNA@data[is.na(data_current[[j]]@assays$RNA@data)] <- 0
      data_current[[j]] <- FindVariableFeatures(data_current[[j]],
                                                nfeatures = seurat.findvariablefeatures.nfeatures,
                                                selection.method = seurat.findvariablefeatures.selection.method,
                                                loess.span = seurat.findvariablefeatures.loess.span,
                                                clip.max = seurat.findvariablefeatures.clip.max,
                                                mean.function = seurat.findvariablefeatures.mean.function,
                                                dispersion.function = seurat.findvariablefeatures.dispersion.function,
                                                num.bin = seurat.findvariablefeatures.num.bin,
                                                binning.method = seurat.findvariablefeatures.binning.method,
                                                verbose = TRUE)
      data_current[[j]] <- ScaleData(data_current[[j]],
                                     features = seurat.scaledata.features,
                                     vars.to.regress = seurat.scaledata.vars.to.regress,
                                     # latent.data = seurat.scaledata.latent.data,
                                     split.by = seurat.scaledata.split.by,
                                     model.use = seurat.scaledata.model.use,
                                     use.umi = seurat.scaledata.use.umi,
                                     do.scale = seurat.scaledata.do.scale,
                                     do.center = seurat.scaledata.do.center,
                                     scale.max = seurat.scaledata.scale.max,
                                     block.size = seurat.scaledata.block.size,
                                     min.cells.to.block = seurat.scaledata.min.cells.to.block,
                                     verbose = TRUE)

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

        data_current[[j]] <- CellCycleScoring(data_current[[j]],
                                              s.features = cc.s.genes,
                                              g2m.features = cc.g2m.genes,
                                              ctrl = seurat.cellcyclescoring.ctrl,
                                              set.ident = seurat.cellcyclescoring.set.ident)

        data_current[[j]]@meta.data$Phase <- factor(data_current[[j]]@meta.data$Phase, levels = c("G1","S","G2M"))

        if((cc_regression == 1)){
          cc_method <- c("S.Score", "G2M.Score")
          cc_type <- "Phase Scores"

        }else if (cc_regression == 2){
          data_current[[j]]$CC.Difference <- data_current[[j]]$S.Score - data_current[[j]]$G2M.Score
          cc_method <- "CC.Difference"
          cc_type <- "Difference"
        }

        data_current[[j]] <- ScaleData(data_current[[j]],
                                       vars.to.regress = cc_method,
                                       features = seurat.scaledata.features,
                                       # latent.data = seurat.scaledata.latent.data,
                                       split.by = seurat.scaledata.split.by,
                                       model.use = seurat.scaledata.model.use,
                                       use.umi = seurat.scaledata.use.umi,
                                       do.scale = seurat.scaledata.do.scale,
                                       do.center = seurat.scaledata.do.center,
                                       scale.max = seurat.scaledata.scale.max,
                                       block.size = seurat.scaledata.block.size,
                                       min.cells.to.block = seurat.scaledata.min.cells.to.block,
                                       verbose = TRUE)
        data_current[[j]] <- RunPCA(data_current[[j]],
                                    features = c(cc.s.genes, cc.g2m.genes),
                                    npcs = seurat.cellcyclescoring.seurat.runpca.npcs,
                                    rev.pca = seurat.cellcyclescoring.seurat.runpca.rev.pca,
                                    weight.by.var = seurat.cellcyclescoring.seurat.runpca.weight.by.var,
                                    ndims.print = 1:seurat.cellcyclescoring.seurat.runpca.ndims.print,
                                    nfeatures.print = seurat.cellcyclescoring.seurat.runpca.nfeatures.print,
                                    reduction.key = seurat.cellcyclescoring.seurat.runpca.reduction.key,
                                    seed.use = seurat.cellcyclescoring.seurat.runpca.seed.use,
                                    approx = seurat.cellcyclescoring.seurat.runpca.approx,
                                    verbose = TRUE)
        data_current[[j]]@reductions$pca_selected <- data_current[[j]]@reductions$pca

        if(length(grep("ENSG[0-9]+",row.names(data_current)[j], ignore.case = T)) > nrow(data_current[[j]])/2){
          row.names(data_current[[j]]@assays$RNA@counts) <- orig_gene_names
          row.names(data_current[[j]]@assays$RNA@data) <- orig_gene_names
          row.names(data_current[[j]]@assays$RNA@scale.data) <- scale_orig_gene_names
        }

        data_current[[j]] <- RunPCA(data_current[[j]],
                                    features = seurat.runpca.features,
                                    npcs = seurat.runpca.npcs,
                                    rev.pca = seurat.runpca.rev.pca,
                                    weight.by.var = seurat.runpca.weight.by.var,
                                    ndims.print = 1:seurat.runpca.ndims.print,
                                    nfeatures.print = seurat.runpca.nfeatures.print,
                                    reduction.key = seurat.runpca.reduction.key,
                                    seed.use = seurat.runpca.seed.use,
                                    approx = seurat.runpca.approx,
                                    verbose = TRUE)
        
        selected_pcs <- NULL
        selected_pcs <- find_selected_pcs(data_current[[j]],
                                          pcs = seurat.runumap.pcs,
                                          pca = "pca_selected",
                                          cname = annot_names[j])
        
        data_current[[j]] <- RunUMAP(data_current[[j]],
                                     reduction = "pca_selected",
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
                                     reduction.key = seurat.runumap.reduction.key,
                                     verbose = TRUE)
        data_current[[j]]@reductions$umap_selected <- data_current[[j]]@reductions$umap
        data_current[[j]] <- RunUMAP(data_current[[j]],
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
                                     reduction.key = seurat.runumap.reduction.key,
                                     verbose = TRUE)

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
        
        selected_pcs <- NULL
        selected_pcs <- find_selected_pcs(data_current[[j]],
                                          pcs = seurat.runumap.pcs,
                                          pca = "pca",
                                          cname = annot_names[j])
        
        DefaultAssay(data_current[[j]]) <- 'RNA'
        data_current[[j]] <- RunPCA(data_current[[j]],
                                    features = seurat.runpca.features,
                                    npcs = seurat.runpca.npcs,
                                    rev.pca = seurat.runpca.rev.pca,
                                    weight.by.var = seurat.runpca.weight.by.var,
                                    ndims.print = 1:seurat.runpca.ndims.print,
                                    nfeatures.print = seurat.runpca.nfeatures.print,
                                    reduction.key = seurat.runpca.reduction.key,
                                    seed.use = seurat.runpca.seed.use,
                                    approx = seurat.runpca.approx,
                                    verbose = TRUE)
        data_current[[j]] <- RunUMAP(data_current[[j]],
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
                                     reduction.key = seurat.runumap.reduction.key,
                                     verbose = TRUE)

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
        temp <- NormalizeData(data_current[[j]],
                              normalization.method = seurat.normalizedata.normalization.method,
                              scale.factor = seurat.normalizedata.norm.scale.factor,
                              margin = seurat.normalizedata.norm.margin,
                              verbose = TRUE)
        data_current[[j]]@commands$NormalizeData.RNA <- temp@commands$NormalizeData.RNA
      }
      
      if(find_doublet == "YES"){
      temp <- NULL
      temp <- data_current[[j]]
      sweep.res.list <- paramSweep_v3(temp,
                                      PCs = 1:doubletfinder.paramSweep_v3.PCs,
                                      sct = doubletfinder.paramSweep_v3.sct)
      sweep.stats <- summarizeSweep(sweep.res.list,
                                    GT = doubletfinder.summarizeSweep.gt,
                                    GT.calls = doubletfinder.summarizeSweep.classifications)
      bcmvn <- find.pK(sweep.stats)
      bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
      cpK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),"pK"]))
      nExp_poi <- round(doubletfinder.homotypic.prop*nrow(temp@meta.data))
      temp <- doubletFinder_v3(temp,
                               PCs = 1:doubletfinder.doubletFinder_v3.PCs,
                               pK = cpK,
                               nExp = nExp_poi,
                               pN = doubletfinder.doubletFinder_v3.pN,
                               reuse.pANN = doubletfinder.doubletFinder_v3.reuse.pANN,
                               sct = doubletfinder.doubletFinder_v3.sct,
                               annotations = doubletfinder.doubletFinder_v3.annotations)
      data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, temp@meta.data[,grep("^pANN_|^DF\\.classifications_", colnames(temp@meta.data), ignore.case = T)])
      cdoublets <- NULL
      cclassname <- NULL
      cclassname <- colnames(data_current[[j]]@meta.data)[grep("^DF\\.classifications_", colnames(data_current[[j]]@meta.data), ignore.case = T)]
      cdoublets <- table(data_current[[j]]@meta.data[,cclassname])
      cdoublets <- cdoublets[grep("Doublet", names(cdoublets), ignore.case = T)]
      print(paste(annot_names[j],": found ",cdoublets, " doublets out of ", ncol(data_current[[j]])," cells..removing..", sep = ""))
      Idents(data_current[[j]]) <- cclassname
      data_current[[j]] <- subset(data_current[[j]], cells = row.names(data_current[[j]]@meta.data[which(data_current[[j]]@meta.data[,cclassname] == "Singlet"),]))
    }
    
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
      data_current[[j]] <- RunPCA(data_current[[j]],
                                  features = seurat.runpca.features,
                                  npcs = seurat.runpca.npcs,
                                  rev.pca = seurat.runpca.rev.pca,
                                  weight.by.var = seurat.runpca.weight.by.var,
                                  ndims.print = 1:seurat.runpca.ndims.print,
                                  nfeatures.print = seurat.runpca.nfeatures.print,
                                  reduction.key = seurat.runpca.reduction.key,
                                  seed.use = seurat.runpca.seed.use,
                                  approx = seurat.runpca.approx,
                                  verbose = TRUE)

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
      
      selected_pcs <- NULL
      selected_pcs <- find_selected_pcs(data_current[[j]],
                                        pcs = seurat.runumap.pcs,
                                        pca = "pca",
                                        cname = annot_names[j])

      data_current[[j]] <- FindNeighbors(data_current[[j]],
                                         dims = 1:selected_pcs,
                                         reduction = seurat.findneighbors.reduction,
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
      data_current[[j]] <- FindClusters(data_current[[j]],
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
                                        edge.file.name = seurat.findclusters.edge.file.name)
      data_current[[j]] <- RunTSNE(data_current[[j]],
                                   reduction = seurat.runtsne.reduction,
                                   cells = seurat.runtsne.cells,
                                   dims = 1:selected_pcs,
                                   features = seurat.runtsne.features,
                                   seed.use = seurat.runtsne.seed.use,
                                   tsne.method = seurat.runtsne.tsne.method,
                                   dim.embed = seurat.runtsne.dim.embed,
                                   distance.matrix = seurat.runtsne.distance.matrix,
                                   reduction.name = seurat.runtsne.reduction.name,
                                   reduction.key = seurat.runtsne.reduction.key,
                                   check_duplicates = seurat.runtsne.check_duplicates)
      
      data_current[[j]] <- RunUMAP(data_current[[j]],
                                   reduction = seurat.runumap.reduction,
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
                                   reduction.key = seurat.runumap.reduction.key,
                                   verbose = TRUE)

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
      current_out <- deanalysis(data_current[[j]],
                                current_clusters,
                                plot_title = annot_names[j],
                                group=NULL,
                                de_analysis = "findallmarkers",
                                deanalysis.p_threshold = 0.05,
                                
                                # FindAllMarkers
                                seurat.findallmarkers.min.pct = seurat.findallmarkers.min.pct,
                                seurat.findallmarkers.logfc.threshold = seurat.findallmarkers.logfc.threshold,
                                seurat.findallmarkers.assay = seurat.findallmarkers.assay,
                                seurat.findallmarkers.features = seurat.findallmarkers.features,
                                seurat.findallmarkers.test.use = seurat.findallmarkers.test.use, #
                                seurat.findallmarkers.slot = seurat.findallmarkers.slot,
                                seurat.findallmarkers.min.diff.pct = seurat.findallmarkers.min.diff.pct,
                                seurat.findallmarkers.node = seurat.findallmarkers.node,
                                seurat.findallmarkers.only.pos = seurat.findallmarkers.only.pos,
                                seurat.findallmarkers.max.cells.per.ident = seurat.findallmarkers.max.cells.per.ident,
                                seurat.findallmarkers.random.seed = seurat.findallmarkers.random.seed,
                                seurat.findallmarkers.latent.vars = seurat.findallmarkers.latent.vars,
                                seurat.findallmarkers.min.cells.feature = seurat.findallmarkers.min.cells.feature,
                                seurat.findallmarkers.min.cells.group = seurat.findallmarkers.min.cells.group,
                                seurat.findallmarkers.mean.fxn = seurat.findallmarkers.mean.fxn,
                                seurat.findallmarkers.fc.name = seurat.findallmarkers.fc.name,
                                seurat.findallmarkers.base = seurat.findallmarkers.base,
                                seurat.findallmarkers.return.thresh = seurat.findallmarkers.return.thresh,
                                seurat.findallmarkers.densify = seurat.findallmarkers.densify,
                                
                                # FindConservedMarkers
                                seurat.findconservedmarkers.ident.2 = seurat.findconservedmarkers.ident.2,
                                seurat.findconservedmarkers.assay = seurat.findconservedmarkers.assay,
                                seurat.findconservedmarkers.slot = seurat.findconservedmarkers.slot,
                                seurat.findconservedmarkers.min.cells.group = seurat.findconservedmarkers.min.cells.group,
                                seurat.findconservedmarkers.meta.method = seurat.findconservedmarkers.meta.method,
                                
                                # FindMarkers
                                seurat.findmarkers.group.by = seurat.findmarkers.group.by,
                                seurat.findmarkers.subset.ident = seurat.findmarkers.subset.ident,
                                seurat.findmarkers.assay = seurat.findmarkers.assay,
                                seurat.findmarkers.slot = seurat.findmarkers.slot,
                                seurat.findmarkers.reduction = seurat.findmarkers.reduction,
                                seurat.findmarkers.features = seurat.findmarkers.features,
                                seurat.findmarkers.logfc.threshold = seurat.findmarkers.logfc.threshold,
                                seurat.findmarkers.test.use = seurat.findmarkers.test.use,
                                seurat.findmarkers.min.pct = seurat.findmarkers.min.pct,
                                seurat.findmarkers.min.diff.pct = seurat.findmarkers.min.diff.pct,
                                seurat.findmarkers.only.pos = seurat.findmarkers.only.pos,
                                seurat.findmarkers.max.cells.per.ident = seurat.findmarkers.max.cells.per.ident,
                                seurat.findmarkers.random.seed = seurat.findmarkers.random.seed,
                                seurat.findmarkers.latent.vars = seurat.findmarkers.latent.vars,
                                seurat.findmarkers.min.cells.feature = seurat.findmarkers.min.cells.feature,
                                seurat.findmarkers.min.cells.group = seurat.findmarkers.min.cells.group,
                                seurat.findmarkers.mean.fxn = seurat.findmarkers.mean.fxn,
                                seurat.findmarkers.fc.name = seurat.findmarkers.fc.name,
                                seurat.findmarkers.base = seurat.findmarkers.base,
                                seurat.findmarkers.densify = seurat.findmarkers.densify)

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
      cumapplot <- scubi_expression(dim1 = cumap[,"UMAP_1"],
                                    dim2 = cumap[,"UMAP_2"],
                                    count = data_current[[j]]@assays$RNA@counts,
                                    gene = as.list(current_out$top1$gene),
                                    palette = scubi.palette)

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
                      ncol = 4, cols = gen_colors(n = length(unique(data_current[[j]]$seurat_clusters)))) &
                xlab("CLUSTERS") &
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
                         ref = singler.ref,
                         labels = singler.labels.main,
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
      data_current[[j]]$CELL_TYPE <- clu_ann$labels[match(data_current[[j]]$seurat_clusters,row.names(clu_ann))]
      clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data_current[[j]])),
                         clusters =data_current[[j]]$seurat_clusters,
                         ref = singler.ref,
                         labels = singler.labels.fine,
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
        x <- enrichDGN(gene=unique(as.numeric(as.character(x$ENTREZ_ID))),
                       pvalueCutoff = dose.enrichdgn.pvalueCutoff,
                       pAdjustMethod = dose.enrichdgn.pAdjustMethod,
                       qvalueCutoff = dose.enrichdgn.qvalueCutoff,
                       universe = dose.enrichdgn.universe,
                       minGSSize = dose.enrichdgn.minGSSize,
                       maxGSSize = dose.enrichdgn.maxGSSize,
                       readable = dose.enrichdgn.readable)
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
            currentR <- pairwise_termsim(current,
                                         method = enrichplot.pairwise_termsim.method,
                                         semData = enrichplot.pairwise_termsim.semData)
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
      plotx <- compareCluster(current, fun="enrichKEGG", organism="hsa", pvalueCutoff=clusterprofiler.comparecluster.pvaluecutoff)
      currentR <- pairwise_termsim(plotx,
                                   method = enrichplot.pairwise_termsim.method,
                                   semData = enrichplot.pairwise_termsim.semData)

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
          GSEA_result[[i]] <- gseNCG(gene=gene_list,
                                     exponent = dose.gsencg.exponent,
                                     minGSSize = dose.gsencg.minGSSize,
                                     maxGSSize = dose.gsencg.maxGSSize,
                                     pvalueCutoff = dose.gsencg.pvalueCutoff,
                                     pAdjustMethod = dose.gsencg.pAdjustMethod,
                                     seed = dose.gsencg.seed,
                                     by = dose.gsencg.by)
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
              print(gseaplot2(current, geneSetID = which(current@result$ID != "-"), title = paste(unique(current$Description[which(current$Description != "-")]), collapse = "\n"), pvalue_table = TRUE, ES_geom = "line"))
# ggtitle(paste(annot_names[j], "\nGSEA Plot: Cluster ", names(pathway_EA_result)[i],sep = "")) +
#                 theme(plot.title = element_text(face = "bold")))
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
      x <- FindVariableFeatures(x,
                                nfeatures = seurat.findvariablefeatures.nfeatures,
                                selection.method = seurat.findvariablefeatures.selection.method,
                                loess.span = seurat.findvariablefeatures.loess.span,
                                clip.max = seurat.findvariablefeatures.clip.max,
                                mean.function = seurat.findvariablefeatures.mean.function,
                                dispersion.function = seurat.findvariablefeatures.dispersion.function,
                                num.bin = seurat.findvariablefeatures.num.bin,
                                binning.method = seurat.findvariablefeatures.binning.method,
                                verbose = TRUE)
      # x <- RunPCA(x, npcs = 30, verbose = FALSE, features = VariableFeatures(x))
    })
    integrative_features <- SelectIntegrationFeatures(data_current,
                                                      nfeatures = seurat.selectintegrationfeatures.nfeatures,
                                                      assay = seurat.selectintegrationfeatures.assay,
                                                      fvf.nfeatures = seurat.selectintegrationfeatures.fvf.nfeatures)
    data_anchors <- FindIntegrationAnchors(data_current,
                                           reduction = seurat.findintegrationanchors.reduction, # reduction = "rpca",
                                           anchor.features = integrative_features,
                                           assay = seurat.findintegrationanchors.assay,
                                           reference = seurat.findintegrationanchors.reference,
                                           scale = seurat.findintegrationanchors.scale,
                                           normalization.method = seurat.findintegrationanchors.normalization.method,
                                           sct.clip.range = seurat.findintegrationanchors.sct.clip.range,
                                           l2.norm = seurat.findintegrationanchors.l2.norm,
                                           dims = seurat.findintegrationanchors.dims,
                                           k.anchor = seurat.findintegrationanchors.k.anchor,
                                           k.filter = seurat.findintegrationanchors.k.filter,
                                           k.score = seurat.findintegrationanchors.k.score,
                                           max.features = seurat.findintegrationanchors.max.features,
                                           nn.method = seurat.findintegrationanchors.nn.method,
                                           n.trees = seurat.findintegrationanchors.n.trees,
                                           eps = seurat.findintegrationanchors.eps,
                                           verbose = TRUE)
    data <- IntegrateData(anchorset = data_anchors,
                          new.assay.name = "integrated",
                          normalization.method = seurat.integratedata.normalization.method,
                          features = seurat.integratedata.features,
                          features.to.integrate = seurat.integratedata.features.to.integrate,
                          dims = 1:seurat.integratedata.dims,
                          k.weight = seurat.integratedata.k.weight,
                          weight.reduction = seurat.integratedata.weight.reduction,
                          sd.weight = seurat.integratedata.sd.weight,
                          sample.tree = seurat.integratedata.sample.tree,
                          preserve.order = seurat.integratedata.preserve.order,
                          eps = seurat.integratedata.eps,
                          verbose = TRUE)
    rm(data_anchors)
    DefaultAssay(data) <- "integrated"
    data <- ScaleData(data,
                      features = seurat.integration.scaledata.features,
                      vars.to.regress = seurat.integration.scaledata.vars.to.regress,
                      # latent.data = seurat.integration.scaledata.latent.data,
                      split.by = seurat.integration.scaledata.split.by,
                      model.use = seurat.integration.scaledata.model.use,
                      use.umi = seurat.integration.scaledata.use.umi,
                      do.scale = seurat.integration.scaledata.do.scale,
                      do.center = seurat.integration.scaledata.do.center,
                      scale.max = seurat.integration.scaledata.scale.max,
                      block.size = seurat.integration.scaledata.block.size,
                      min.cells.to.block = seurat.integration.scaledata.min.cells.to.block,
                      verbose = TRUE)

    DefaultAssay(data) <- "RNA"
    data <- FindVariableFeatures(data,
                                 nfeatures = seurat.findvariablefeatures.nfeatures,
                                 selection.method = seurat.findvariablefeatures.selection.method,
                                 loess.span = seurat.findvariablefeatures.loess.span,
                                 clip.max = seurat.findvariablefeatures.clip.max,
                                 mean.function = seurat.findvariablefeatures.mean.function,
                                 dispersion.function = seurat.findvariablefeatures.dispersion.function,
                                 num.bin = seurat.findvariablefeatures.num.bin,
                                 binning.method = seurat.findvariablefeatures.binning.method)
    data <- ScaleData(data,
                      features = seurat.scaledata.features,
                      vars.to.regress = seurat.scaledata.vars.to.regress,
                      # latent.data = seurat.scaledata.latent.data,
                      split.by = seurat.scaledata.split.by,
                      model.use = seurat.scaledata.model.use,
                      use.umi = seurat.scaledata.use.umi,
                      do.scale = seurat.scaledata.do.scale,
                      do.center = seurat.scaledata.do.center,
                      scale.max = seurat.scaledata.scale.max,
                      block.size = seurat.scaledata.block.size,
                      min.cells.to.block = seurat.scaledata.min.cells.to.block,
                      verbose = TRUE)
    data <- RunPCA(data,
                   features = seurat.runpca.features,
                   npcs = seurat.runpca.npcs,
                   rev.pca = seurat.runpca.rev.pca,
                   weight.by.var = seurat.runpca.weight.by.var,
                   ndims.print = 1:seurat.runpca.ndims.print,
                   nfeatures.print = seurat.runpca.nfeatures.print,
                   reduction.key = seurat.runpca.reduction.key,
                   seed.use = seurat.runpca.seed.use,
                   approx = seurat.runpca.approx,
                   verbose = TRUE)
    
    selected_pcs <- NULL
    selected_pcs <- find_selected_pcs(data,
                                      pcs = seurat.runumap.pcs,
                                      pca = "pca",
                                      cname = "all samples (before integration)")

    data <- RunUMAP(data,
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
    
    data <- RunTSNE(data,
                    reduction = seurat.runtsne.reduction,
                    cells = seurat.runtsne.cells,
                    dims = 1:selected_pcs,
                    features = seurat.runtsne.features,
                    seed.use = seurat.runtsne.seed.use,
                    tsne.method = seurat.runtsne.tsne.method,
                    dim.embed = seurat.runtsne.dim.embed,
                    distance.matrix = seurat.runtsne.distance.matrix,
                    reduction.name = seurat.runtsne.reduction.name,
                    reduction.key = seurat.runtsne.reduction.key,
                    check_duplicates = seurat.runtsne.check_duplicates)
    
    data_dim <- data.frame(gen10x_plotx(data), DATA_TYPE = "BEFORE_INTEGRATION", SAMPLE_ID = data$orig.ident)

    write.table(data_dim, paste(cdir, "42URSA_TABLE_scRNASEQ_DIM_PARAMETERS_BEFORE_INTEGRATION_", project_name,".txt", sep = ""), quote = F, row.names = T, sep = "\t")

    if((toupper(integration_method) == "SEURAT") | (toupper(integration_method) == "NULL")){
      # integration_method <- "SEURAT"
      integration_name <- "SEURAT_INTEGRATED"
      DefaultAssay(data) <- "integrated"
      reduction_method <- "pca"

      data <- RunPCA(data,
                     features = seurat.integration.runpca.features,
                     npcs = seurat.integration.runpca.npcs,
                     rev.pca = seurat.integration.runpca.rev.pca,
                     weight.by.var = seurat.integration.runpca.weight.by.var,
                     ndims.print = seurat.integration.runpca.ndims.print,
                     nfeatures.print = seurat.integration.runpca.nfeatures.print,
                     reduction.key = seurat.integration.runpca.reduction.key,
                     seed.use = seurat.integration.runpca.seed.use,
                     approx = seurat.integration.runpca.approx,
                     verbose = FALSE)
      
      selected_pcs <- NULL
      selected_pcs <- find_selected_pcs(data,
                                        pcs = seurat.integration.runumap.pcs,
                                        pca = "pca",
                                        cname = "integrated data")
      
      data <- RunUMAP(data,
                      reduction = "pca",
                      dims = 1:selected_pcs,
                      features = seurat.integration.runumap.features,
                      graph = seurat.integration.runumap.graph,
                      assay = seurat.integration.runumap.assay,
      nn.name = seurat.integration.runumap.nn.name,
      slot = seurat.integration.runumap.slot,
      umap.method = seurat.integration.runumap.umap.method,
      reduction.model = seurat.integration.runumap.reduction.model,
      return.model = seurat.integration.runumap.return.model,
      n.neighbors = seurat.integration.runumap.n.neighbors,
      n.components = seurat.integration.runumap.n.components,
      metric = seurat.integration.runumap.metric,
      n.epochs = seurat.integration.runumap.n.epochs,
      learning.rate = seurat.integration.runumap.learning.rate,
      min.dist = seurat.integration.runumap.min.dist,
      spread = seurat.integration.runumap.spread,
      set.op.mix.ratio = seurat.integration.runumap.set.op.mix.ratio,
      local.connectivity = seurat.integration.runumap.local.connectivity,
      repulsion.strength = seurat.integration.runumap.repulsion.strength,
      negative.sample.rate = seurat.integration.runumap.negative.sample.rate,
      a = seurat.integration.runumap.a,
      b = seurat.integration.runumap.b,
      uwot.sgd = seurat.integration.runumap.uwot.sgd,
      seed.use = seurat.integration.runumap.seed.use,
      metric.kwds = seurat.integration.runumap.metric.kwds,
      angular.rp.forest = seurat.integration.runumap.angular.rp.forest,
      densmap = seurat.integration.runumap.densmap,
      dens.lambda = seurat.integration.runumap.dens.lambda,
      dens.frac = seurat.integration.runumap.dens.frac,
      dens.var.shift = seurat.integration.runumap.dens.var.shift,
      reduction.name = seurat.integration.runumap.reduction.name,
      reduction.key = seurat.integration.runumap.reduction.key)

      data <- RunTSNE(data,
                      reduction = "pca",
                      dims = 1:selected_pcs,
                      cells = seurat.integration.runtsne.cells,
                      features = seurat.integration.runtsne.features,
                      seed.use = seurat.integration.runtsne.seed.use,
                      tsne.method = seurat.integration.runtsne.tsne.method,
                      dim.embed = seurat.integration.runtsne.dim.embed,
                      distance.matrix = seurat.integration.runtsne.distance.matrix,
                      reduction.name = seurat.integration.runtsne.reduction.name,
                      reduction.key = seurat.integration.runtsne.reduction.key,
                      check_duplicates = seurat.integration.runtsne.check_duplicates)
      
      current <- cbind(data.frame(gen10x_plotx(data), DATA_TYPE = integration_name, SAMPLE_ID = data$orig.ident))
      
    }else if(toupper(integration_method) == "HARMONY"){
      
      integration_name <- "HARMONY_INTEGRATED"
      DefaultAssay(data) <- "RNA"
      reduction_method <- "harmony"
      data <- RunHarmony(data,
                         group.by.vars = harmony.runharmony.group.by.vars, # BATCH
                         reduction = harmony.runharmony.reduction,
                         dims.use = harmony.runharmony.dims.use,
                         theta = harmony.runharmony.theta,
                         lambda = harmony.runharmony.lambda,
                         sigma = harmony.runharmony.sigma,
                         nclust = harmony.runharmony.nclust,
                         tau = harmony.runharmony.tau,
                         block.size = harmony.runharmony.block.size,
                         max.iter.harmony = harmony.runharmony.max.iter.harmony,
                         max.iter.cluster = harmony.runharmony.max.iter.cluster,
                         epsilon.cluster = harmony.runharmony.epsilon.cluster,
                         epsilon.harmony = harmony.runharmony.epsilon.harmony,
                         plot_convergence = harmony.runharmony.plot_convergence,
                         reference_values = harmony.runharmony.reference_values,
                         reduction.save = harmony.runharmony.reduction.save,
                         assay.use = harmony.runharmony.assay.use,
                         project.dim = harmony.runharmony.project.dim,
                         verbose = FALSE)

      selected_pcs <- NULL
      selected_pcs <- find_selected_pcs(data,
                                        pcs = seurat.integration.runumap.pcs,
                                        pca = "harmony",
                                        cname = "integrated data")

      data <- RunUMAP(data,
                      reduction = "harmony",
                      dims = 1:selected_pcs,
                      features = seurat.integration.runumap.features,
                      graph = seurat.integration.runumap.graph,
                      assay = seurat.integration.runumap.assay,
                      nn.name = seurat.integration.runumap.nn.name,
                      slot = seurat.integration.runumap.slot,
                      umap.method = seurat.integration.runumap.umap.method,
                      reduction.model = seurat.integration.runumap.reduction.model,
                      return.model = seurat.integration.runumap.return.model,
                      n.neighbors = seurat.integration.runumap.n.neighbors,
                      n.components = seurat.integration.runumap.n.components,
                      metric = seurat.integration.runumap.metric,
                      n.epochs = seurat.integration.runumap.n.epochs,
                      learning.rate = seurat.integration.runumap.learning.rate,
                      min.dist = seurat.integration.runumap.min.dist,
                      spread = seurat.integration.runumap.spread,
                      set.op.mix.ratio = seurat.integration.runumap.set.op.mix.ratio,
                      local.connectivity = seurat.integration.runumap.local.connectivity,
                      repulsion.strength = seurat.integration.runumap.repulsion.strength,
                      negative.sample.rate = seurat.integration.runumap.negative.sample.rate,
                      a = seurat.integration.runumap.a,
                      b = seurat.integration.runumap.b,
                      uwot.sgd = seurat.integration.runumap.uwot.sgd,
                      seed.use = seurat.integration.runumap.seed.use,
                      metric.kwds = seurat.integration.runumap.metric.kwds,
                      angular.rp.forest = seurat.integration.runumap.angular.rp.forest,
                      densmap = seurat.integration.runumap.densmap,
                      dens.lambda = seurat.integration.runumap.dens.lambda,
                      dens.frac = seurat.integration.runumap.dens.frac,
                      dens.var.shift = seurat.integration.runumap.dens.var.shift,
                      reduction.name = seurat.integration.runumap.reduction.name,
                      reduction.key = seurat.integration.runumap.reduction.key)
      
      data <- RunTSNE(data,
                      reduction = "harmony",
                      dims = 1:selected_pcs,
                      check_duplicates = FALSE,
                      cells = seurat.integration.runtsne.cells,
                      features = seurat.integration.runtsne.features,
                      seed.use = seurat.integration.runtsne.seed.use,
                      tsne.method = seurat.integration.runtsne.tsne.method,
                      dim.embed = seurat.integration.runtsne.dim.embed,
                      distance.matrix = seurat.integration.runtsne.distance.matrix,
                      reduction.name = seurat.integration.runtsne.reduction.name,
                      reduction.key = seurat.integration.runtsne.reduction.key)
      
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
    
    selected_pcs <- NULL
    selected_pcs <- find_selected_pcs(data,
                                      pcs = seurat.integration.runumap.pcs,
                                      pca = reduction_method,
                                      cname = "integrated data")

    data <- FindNeighbors(data,
                          reduction = reduction_method,
                          dims = 1:selected_pcs,
                          assay = seurat.integration.findneighbors.assay,
                          features = seurat.integration.findneighbors.features,
                          k.param = seurat.integration.findneighbors.k.param,
                          return.neighbor = seurat.integration.findneighbors.return.neighbor,
                          # compute.SNN = seurat.integration.findneighbors.compute.SNN,
                          prune.SNN = seurat.integration.findneighbors.prune.SNN,
                          nn.method = seurat.integration.findneighbors.nn.method,
                          n.trees = seurat.integration.findneighbors.n.trees,
                          annoy.metric = seurat.integration.findneighbors.annoy.metric,
                          nn.eps = seurat.integration.findneighbors.nn.eps,
                          force.recalc = seurat.integration.findneighbors.force.recalc,
                          do.plot = seurat.integration.findneighbors.do.plot,
                          graph.name = seurat.integration.findneighbors.graph.name,
                          l2.norm = seurat.integration.findneighbors.l2.norm,
                          cache.index = seurat.integration.findneighbors.cache.index)
    data <- FindClusters(data,
                         graph.name = seurat.integration.findclusters.graph.name,
                         modularity.fxn = seurat.integration.findclusters.modularity.fxn,
                         initial.membership = seurat.integration.findclusters.initial.membership,
                         node.sizes = seurat.integration.findclusters.node.sizes,
                         resolution = seurat.integration.findclusters.resolution,
                         method = seurat.integration.findclusters.method,
                         algorithm = seurat.integration.findclusters.algorithm,
                         n.start = seurat.integration.findclusters.n.start,
                         n.iter = seurat.integration.findclusters.n.iter,
                         random.seed = seurat.integration.findclusters.random.seed,
                         group.singletons = seurat.integration.findclusters.group.singletons,
                         temp.file.location = seurat.integration.findclusters.temp.file.location,
                         edge.file.name = seurat.integration.findclusters.edge.file.name,
                         verbose = FALSE)

    if(toupper(integration_method) == "HARMONY"){
      integration_cluster <- paste("RNA_snn_res.", seurat.findclusters.resolution, sep = "")
    }else{
      integration_cluster <- paste("integrated_snn_res.", seurat.findclusters.resolution, sep = "")
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
    current_out <- deanalysis(data, current_clusters,
                              plot_title = project_name,
                              group=NULL,
                              de_analysis = "findallmarkers",
                              cluster_name = integration_cluster,
                              deanalysis.p_threshold = 0.05,
                              
                              # FindAllMarkers
                              seurat.findallmarkers.min.pct = seurat.findallmarkers.min.pct,
                              seurat.findallmarkers.logfc.threshold = seurat.findallmarkers.logfc.threshold,
                              seurat.findallmarkers.assay = seurat.findallmarkers.assay,
                              seurat.findallmarkers.features = seurat.findallmarkers.features,
                              seurat.findallmarkers.test.use = seurat.findallmarkers.test.use,
                              seurat.findallmarkers.slot = seurat.findallmarkers.slot,
                              seurat.findallmarkers.min.diff.pct = seurat.findallmarkers.min.diff.pct,
                              seurat.findallmarkers.node = seurat.findallmarkers.node,
                              seurat.findallmarkers.only.pos = seurat.findallmarkers.only.pos,
                              seurat.findallmarkers.max.cells.per.ident = seurat.findallmarkers.max.cells.per.ident,
                              seurat.findallmarkers.random.seed = seurat.findallmarkers.random.seed,
                              seurat.findallmarkers.latent.vars = seurat.findallmarkers.latent.vars,
                              seurat.findallmarkers.min.cells.feature = seurat.findallmarkers.min.cells.feature,
                              seurat.findallmarkers.min.cells.group = seurat.findallmarkers.min.cells.group,
                              seurat.findallmarkers.mean.fxn = seurat.findallmarkers.mean.fxn,
                              seurat.findallmarkers.fc.name = seurat.findallmarkers.fc.name,
                              seurat.findallmarkers.base = seurat.findallmarkers.base,
                              seurat.findallmarkers.return.thresh = seurat.findallmarkers.return.thresh,
                              seurat.findallmarkers.densify = seurat.findallmarkers.densify,
                              
                              # FindConservedMarkers
                              seurat.findconservedmarkers.ident.2 = seurat.findconservedmarkers.ident.2,
                              seurat.findconservedmarkers.assay = seurat.findconservedmarkers.assay,
                              seurat.findconservedmarkers.slot = seurat.findconservedmarkers.slot,
                              seurat.findconservedmarkers.min.cells.group = seurat.findconservedmarkers.min.cells.group,
                              seurat.findconservedmarkers.meta.method = seurat.findconservedmarkers.meta.method,
                              
                              # FindMarkers
                              seurat.findmarkers.group.by = seurat.findmarkers.group.by,
                              seurat.findmarkers.subset.ident = seurat.findmarkers.subset.ident,
                              seurat.findmarkers.assay = seurat.findmarkers.assay,
                              seurat.findmarkers.slot = seurat.findmarkers.slot,
                              seurat.findmarkers.reduction = seurat.findmarkers.reduction,
                              seurat.findmarkers.features = seurat.findmarkers.features,
                              seurat.findmarkers.logfc.threshold = seurat.findmarkers.logfc.threshold,
                              seurat.findmarkers.test.use = seurat.findmarkers.test.use,
                              seurat.findmarkers.min.pct = seurat.findmarkers.min.pct,
                              seurat.findmarkers.min.diff.pct = seurat.findmarkers.min.diff.pct,
                              seurat.findmarkers.only.pos = seurat.findmarkers.only.pos,
                              seurat.findmarkers.max.cells.per.ident = seurat.findmarkers.max.cells.per.ident,
                              seurat.findmarkers.random.seed = seurat.findmarkers.random.seed,
                              seurat.findmarkers.latent.vars = seurat.findmarkers.latent.vars,
                              seurat.findmarkers.min.cells.feature = seurat.findmarkers.min.cells.feature,
                              seurat.findmarkers.min.cells.group = seurat.findmarkers.min.cells.group,
                              seurat.findmarkers.mean.fxn = seurat.findmarkers.mean.fxn,
                              seurat.findmarkers.fc.name = seurat.findmarkers.fc.name,
                              seurat.findmarkers.base = seurat.findmarkers.base,
                              seurat.findmarkers.densify = seurat.findmarkers.densify)
    
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
                       ref = singler.ref,
                       labels = singler.labels.main,
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
    
    data$INTEGRATED_CELL_TYPE <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
    data@meta.data[which(is.na(data$INTEGRATED_CELL_TYPE)),"INTEGRATED_CELL_TYPE"] <- "Unidentifiable"

    # clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
    #                    clusters =  data$seurat_clusters,
    #                    ref = hpca.se, assay.type.test=1,
    #                    labels = hpca.se$label.main)
    # data$INTEGRATED_CELL_TYPE_LEVEL2 <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
    # data@meta.data[which(is.na(data$INTEGRATED_CELL_TYPE_LEVEL2)),"CELL_TYPE"] <- "Unidentifiable"

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
                                               color_by = "INTEGRATED_CELL_TYPE", xlabel = "UMAP_1", ylabel = "UMAP_2",
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
                                  ligand_receptor_pairs = celltalker.celltalk.ligand_receptor_pairs,
                                  number_cells_required = celltalker.celltalk.number_cells_required,
                                  min_expression = celltalker.celltalk.min_expression,
                                  max_expression = celltalker.celltalk.max_expression,
                                  scramble_times = celltalker.celltalk.scramble_times)

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
    total_counts <- data.frame(table(current$INTEGRATED_CELL_TYPE))

    current <- data.frame(table(current[,c("SID","INTEGRATED_CELL_TYPE")]))
    current <- current[current$Freq != 0,]
    colnames(current) <- c("SID","INTEGRATED_CELL_TYPE","COUNT")
    current$CELL_TYPE_COUNT <- total_counts[match(current$INTEGRATED_CELL_TYPE, total_counts$Var1),"Freq"]
    current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT

    somePDFPath = paste(cdir,"70URSA_PLOT_scRNASEQ_CELL_TYPE_SUMMARY_SAMPLE_PROPORTION_",integration_name, "_",project_name,".pdf", sep = "")
    pdf(somePDFPath, width = 12, height = 8, pointsize = 10)
    print(ggplot(current,aes(SID,PROPORTION,fill=INTEGRATED_CELL_TYPE))+
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
    current <- group_medianexpr(current_out$current_data_markers, data, group = "INTEGRATED_CELL_TYPE", cell_type = T)
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
    cell_types <- sort(as.character(unique(top_markers_cell_type$INTEGRATED_CELL_TYPE)))

    p1 <- ggplot(data=data.frame(top_markers_cell_type), aes(x=gene, y=data.frame(top_markers_cell_type)[,wt], fill=CELL_TYPE)) +
      geom_bar(stat="identity", position=position_dodge())+
      scale_fill_manual(values = ct_colors)+
      theme_classic() +
      coord_flip()+
      ylab(wt)+
      facet_wrap(~CELL_TYPE, scales = "free_y") +
      theme(axis.text.y=element_text(size = 7))

    plotx <- gen10x_plotx(data, groups = NULL)
    plotx$INTEGRATED_CELL_TYPE <- data$INTEGRATED_CELL_TYPE

    p2 <- plot_bygroup(plotx, x = "tSNE_1", y="tSNE_2", group = "INTEGRATED_CELL_TYPE",
                       plot_title = integration_name,col = ct_colors,numeric = F,
                       annot = T, legend_position = "right", point_size = 1)
    p3 <- plot_bygroup(plotx, x = "UMAP_1", y="UMAP_2", group = "INTEGRATED_CELL_TYPE",
                       plot_title = integration_name,col = ct_colors,numeric = F,
                       annot = T, legend_position = "right", point_size = 1)

    somePDFPath = paste(cdir,"72URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_TSNE_CELL_TYPES_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=24, height=ceiling(length(unique(plotx$INTEGRATED_CELL_TYPE))/3)*5,pointsize=12)
    grid.arrange(p1, p2, ncol = 2)
    dev.off()

    somePDFPath = paste(cdir,"73URSA_PLOT_scRNASEQ_AVE_LOGFC_TOP_GENES_UMAP_CELL_TYPES_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=24, height=ceiling(length(unique(plotx$INTEGRATED_CELL_TYPE))/3)*5,pointsize=12)
    grid.arrange(p1, p3, ncol = 2)
    dev.off()

    current_out$current_data_markers$INTEGRATED_CELL_TYPE <- data@meta.data[match(current_out$current_data_markers$cluster, data@meta.data[,integration_cluster]),"INTEGRATED_CELL_TYPE"]
    top_1_markers <- current_out$current_data_markers %>% group_by(INTEGRATED_CELL_TYPE) %>% top_n(n = 1, eval(parse(text = wt)))
    Idents(data) <- "INTEGRATED_CELL_TYPE"

    somePDFPath = paste(cdir,"74URSA_PLOT_scRNASEQ_RIDGE_TOP_FC_SIG_GENES_CELL_TYPES_",integration_name,"_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=16, height=length(unique(unlist(top_1_markers$gene)))/2*6,pointsize=12)
    print(RidgePlot(data, features = unique(unlist(top_1_markers$gene)), ncol = 2,
                                       cols = ct_colors)+
      plot_annotation(title = integration_name, theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
    dev.off()

    current <- data@meta.data
    total_counts <- data.frame(table(current$INTEGRATED_CELL_TYPE))

    current <- data.frame(table(current[,c("SID","INTEGRATED_CELL_TYPE")]))
    current <- current[current$Freq != 0,]
    colnames(current) <- c("SID","INTEGRATED_CELL_TYPE","COUNT")
    current$CELL_TYPE_COUNT <- total_counts[match(current$INTEGRATED_CELL_TYPE, total_counts$Var1),"Freq"]
    current$PROPORTION <- current$COUNT/current$CELL_TYPE_COUNT

    node_proportion <- current

    somePDFPath = paste(cdir,"75URSA_PLOT_scRNASEQ_CELL_TYPE_SAMPLE_PROPORTION_",integration_name, "_",project_name,".pdf", sep = "")
    pdf(somePDFPath, width = 25, height = 20, pointsize = 10)
    print(ggplot(node_proportion, aes(INTEGRATED_CELL_TYPE, PROPORTION, fill = SID))+
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
      current_out <- deanalysis(data, current_clusters,
                                plot_title = project_name,
                                group="GROUP",
                                de_analysis = "finddegroup",
                                deanalysis.p_threshold = 0.05,
                                
                                # FindAllMarkers
                                seurat.findallmarkers.min.pct = seurat.findallmarkers.min.pct,
                                seurat.findallmarkers.logfc.threshold = seurat.findallmarkers.logfc.threshold,
                                seurat.findallmarkers.assay = seurat.findallmarkers.assay,
                                seurat.findallmarkers.features = seurat.findallmarkers.features,
                                seurat.findallmarkers.test.use = seurat.findallmarkers.test.use,
                                seurat.findallmarkers.slot = seurat.findallmarkers.slot,
                                seurat.findallmarkers.min.diff.pct = seurat.findallmarkers.min.diff.pct,
                                seurat.findallmarkers.node = seurat.findallmarkers.node,
                                seurat.findallmarkers.only.pos = seurat.findallmarkers.only.pos,
                                seurat.findallmarkers.max.cells.per.ident = seurat.findallmarkers.max.cells.per.ident,
                                seurat.findallmarkers.random.seed = seurat.findallmarkers.random.seed,
                                seurat.findallmarkers.latent.vars = seurat.findallmarkers.latent.vars,
                                seurat.findallmarkers.min.cells.feature = seurat.findallmarkers.min.cells.feature,
                                seurat.findallmarkers.min.cells.group = seurat.findallmarkers.min.cells.group,
                                seurat.findallmarkers.mean.fxn = seurat.findallmarkers.mean.fxn,
                                seurat.findallmarkers.fc.name = seurat.findallmarkers.fc.name,
                                seurat.findallmarkers.base = seurat.findallmarkers.base,
                                seurat.findallmarkers.return.thresh = seurat.findallmarkers.return.thresh,
                                seurat.findallmarkers.densify = seurat.findallmarkers.densify,
                                
                                # FindConservedMarkers
                                seurat.findconservedmarkers.ident.2 = seurat.findconservedmarkers.ident.2,
                                seurat.findconservedmarkers.assay = seurat.findconservedmarkers.assay,
                                seurat.findconservedmarkers.slot = seurat.findconservedmarkers.slot,
                                seurat.findconservedmarkers.min.cells.group = seurat.findconservedmarkers.min.cells.group,
                                seurat.findconservedmarkers.meta.method = seurat.findconservedmarkers.meta.method,
                                
                                # FindMarkers
                                seurat.findmarkers.group.by = seurat.findmarkers.group.by,
                                seurat.findmarkers.subset.ident = seurat.findmarkers.subset.ident,
                                seurat.findmarkers.assay = seurat.findmarkers.assay,
                                seurat.findmarkers.slot = seurat.findmarkers.slot,
                                seurat.findmarkers.reduction = seurat.findmarkers.reduction,
                                seurat.findmarkers.features = seurat.findmarkers.features,
                                seurat.findmarkers.logfc.threshold = seurat.findmarkers.logfc.threshold,
                                seurat.findmarkers.test.use = seurat.findmarkers.test.use,
                                seurat.findmarkers.min.pct = seurat.findmarkers.min.pct,
                                seurat.findmarkers.min.diff.pct = seurat.findmarkers.min.diff.pct,
                                seurat.findmarkers.only.pos = seurat.findmarkers.only.pos,
                                seurat.findmarkers.max.cells.per.ident = seurat.findmarkers.max.cells.per.ident,
                                seurat.findmarkers.random.seed = seurat.findmarkers.random.seed,
                                seurat.findmarkers.latent.vars = seurat.findmarkers.latent.vars,
                                seurat.findmarkers.min.cells.feature = seurat.findmarkers.min.cells.feature,
                                seurat.findmarkers.min.cells.group = seurat.findmarkers.min.cells.group,
                                seurat.findmarkers.mean.fxn = seurat.findmarkers.mean.fxn,
                                seurat.findmarkers.fc.name = seurat.findmarkers.fc.name,
                                seurat.findmarkers.base = seurat.findmarkers.base,
                                seurat.findmarkers.densify = seurat.findmarkers.densify)
      
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
    selected_pcs <- find_selected_pcs(data_current[[j]],
                                      pcs = seurat.runumap.pcs,
                                      pca = "pca",
                                      cname = annot_names[j])

    data_current[[j]] <- RunUMAP(data_current[[j]],
                                 n.components = 3,
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

    data_current[[j]] <- RunTSNE(data_current[[j]],
                                 dim.embed = 3,
                                 reduction = seurat.runtsne.reduction,
                                 cells = seurat.runtsne.cells,
                                 dims = 1:selected_pcs,
                                 features = seurat.runtsne.features,
                                 seed.use = seurat.runtsne.seed.use,
                                 tsne.method = seurat.runtsne.tsne.method,
                                 distance.matrix = seurat.runtsne.distance.matrix,
                                 reduction.name = seurat.runtsne.reduction.name,
                                 reduction.key = seurat.runtsne.reduction.key,
                                 check_duplicates = seurat.runtsne.check_duplicates)

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

    data_current[[j]] <- RunTSNE(data_current[[j]],
                                 reduction = seurat.runtsne.reduction,
                                 cells = seurat.runtsne.cells,
                                 dims = 1:seurat.runumap.pcs,
                                 features = seurat.runtsne.features,
                                 seed.use = seurat.runtsne.seed.use,
                                 tsne.method = seurat.runtsne.tsne.method,
                                 dim.embed = seurat.runtsne.dim.embed,
                                 distance.matrix = seurat.runtsne.distance.matrix,
                                 reduction.name = seurat.runtsne.reduction.name,
                                 reduction.key = seurat.runtsne.reduction.key,
                                 check_duplicates = seurat.runtsne.check_duplicates)
    
    data_current[[j]] <- RunUMAP(data_current[[j]],
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
  }

  if(length(data_current) > 1){

    pseudo_para <- c("PSEUDOTIME_PC1","PSEUDOTIME_UMAP1","PSEUDOTIME_tSNE1")
    chosen_para <- c(integration_cluster, "INTEGRATED_CELL_TYPE","orig.ident")
    chosen_para_names <- c("CLUSTERS", "INTEGRATED_CELL_TYPE","SAMPLE_ID")
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
                                         group = "INTEGRATED_CELL_TYPE", label_size = 6,
                                         # cell_size = 0.5,traj_size = 0.8,
                                         paste(project_name,"\nINTEGRATED TRAJECTORY - COLOR BY CELL TYPE"),
                                         color_conditions$tenx, length(unique(mono3data$INTEGRATED_CELL_TYPE)))+
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
      selected_pcs <- find_selected_pcs(data,
                                        pcs = seurat.integration.runumap.features,
                                        pca = "pca",
                                        cname = "integrated data")

      data <- RunTSNE(data,
                      dim.embed = 3,
                      reduction = "harmony",
                      dims = 1:selected_pcs,
                      check_duplicates = FALSE,
                      cells = seurat.integration.runtsne.cells,
                      features = seurat.integration.runtsne.features,
                      seed.use = seurat.integration.runtsne.seed.use,
                      tsne.method = seurat.integration.runtsne.tsne.method,
                      distance.matrix = seurat.integration.runtsne.distance.matrix,
                      reduction.name = seurat.integration.runtsne.reduction.name,
                      reduction.key = seurat.integration.runtsne.reduction.key)
      
      data <- RunUMAP(data,
                      n.components = 3,
                      reduction = "harmony",
                      dims = 1:selected_pcs,
                      features = seurat.integration.runumap.features,
                      graph = seurat.integration.runumap.graph,
                      assay = seurat.integration.runumap.assay,
      nn.name = seurat.integration.runumap.nn.name,
      slot = seurat.integration.runumap.slot,
      umap.method = seurat.integration.runumap.umap.method,
      reduction.model = seurat.integration.runumap.reduction.model,
      return.model = seurat.integration.runumap.return.model,
      n.neighbors = seurat.integration.runumap.n.neighbors,
      metric = seurat.integration.runumap.metric,
      n.epochs = seurat.integration.runumap.n.epochs,
      learning.rate = seurat.integration.runumap.learning.rate,
      min.dist = seurat.integration.runumap.min.dist,
      spread = seurat.integration.runumap.spread,
      set.op.mix.ratio = seurat.integration.runumap.set.op.mix.ratio,
      local.connectivity = seurat.integration.runumap.local.connectivity,
      repulsion.strength = seurat.integration.runumap.repulsion.strength,
      negative.sample.rate = seurat.integration.runumap.negative.sample.rate,
      a = seurat.integration.runumap.a,
      b = seurat.integration.runumap.b,
      uwot.sgd = seurat.integration.runumap.uwot.sgd,
      seed.use = seurat.integration.runumap.seed.use,
      metric.kwds = seurat.integration.runumap.metric.kwds,
      angular.rp.forest = seurat.integration.runumap.angular.rp.forest,
      densmap = seurat.integration.runumap.densmap,
      dens.lambda = seurat.integration.runumap.dens.lambda,
      dens.frac = seurat.integration.runumap.dens.frac,
      dens.var.shift = seurat.integration.runumap.dens.var.shift,
      reduction.name = seurat.integration.runumap.reduction.name,
      reduction.key = seurat.integration.runumap.reduction.key)

    }else if(reduction_method == "pca"){
      
      selected_pcs <- NULL
      selected_pcs <- find_selected_pcs(data,
                                        pcs = seurat.integration.runumap.pcs,
                                        pca = "pca",
                                        cname = "integrated data")
      
      DefaultAssay(data) <- 'integrated'
      data <- RunTSNE(data,
                      dim.embed = 3,
                      reduction = "pca",
                      dims = 1:selected_pcs,
                      cells = seurat.integration.runtsne.cells,
                      features = seurat.integration.runtsne.features,
                      seed.use = seurat.integration.runtsne.seed.use,
                      tsne.method = seurat.integration.runtsne.tsne.method,
                      distance.matrix = seurat.integration.runtsne.distance.matrix,
                      reduction.name = seurat.integration.runtsne.reduction.name,
                      reduction.key = seurat.integration.runtsne.reduction.key,
                      check_duplicates = seurat.integration.runtsne.check_duplicates)
      data <- RunUMAP(data,
                      n.components = 3,
                      reduction = "pca",
                      dims = 1:selected_pcs,
                      features = seurat.integration.runumap.features,
                      graph = seurat.integration.runumap.graph,
                      assay = seurat.integration.runumap.assay,
                      nn.name = seurat.integration.runumap.nn.name,
                      slot = seurat.integration.runumap.slot,
                      umap.method = seurat.integration.runumap.umap.method,
                      reduction.model = seurat.integration.runumap.reduction.model,
                      return.model = seurat.integration.runumap.return.model,
                      n.neighbors = seurat.integration.runumap.n.neighbors,
                      metric = seurat.integration.runumap.metric,
                      n.epochs = seurat.integration.runumap.n.epochs,
                      learning.rate = seurat.integration.runumap.learning.rate,
                      min.dist = seurat.integration.runumap.min.dist,
                      spread = seurat.integration.runumap.spread,
                      set.op.mix.ratio = seurat.integration.runumap.set.op.mix.ratio,
                      local.connectivity = seurat.integration.runumap.local.connectivity,
                      repulsion.strength = seurat.integration.runumap.repulsion.strength,
                      negative.sample.rate = seurat.integration.runumap.negative.sample.rate,
                      a = seurat.integration.runumap.a,
                      b = seurat.integration.runumap.b,
                      uwot.sgd = seurat.integration.runumap.uwot.sgd,
                      seed.use = seurat.integration.runumap.seed.use,
                      metric.kwds = seurat.integration.runumap.metric.kwds,
                      angular.rp.forest = seurat.integration.runumap.angular.rp.forest,
                      densmap = seurat.integration.runumap.densmap,
                      dens.lambda = seurat.integration.runumap.dens.lambda,
                      dens.frac = seurat.integration.runumap.dens.frac,
                      dens.var.shift = seurat.integration.runumap.dens.var.shift,
                      reduction.name = seurat.integration.runumap.reduction.name,
                      reduction.key = seurat.integration.runumap.reduction.key)
    }

    p <- NULL
    p <- plot_3d(data@reductions$umap@cell.embeddings[,"UMAP_1"],
                                     data@reductions$umap@cell.embeddings[,"UMAP_2"],
                                     data@reductions$umap@cell.embeddings[,"UMAP_3"],
                                     color_group = data@meta.data[,"INTEGRATED_CELL_TYPE"],
                                     plot_title = names(data_current)[j],
                                     col = color_conditions$tenx,
                                     n = length(unique(data$INTEGRATED_CELL_TYPE)),
                                     lt = "CELL TYPE", t1 = "UMAP_1",t2 = "UMAP_2", t3= "UMAP_3",
                                     current_text= paste("Cluster ",data@meta.data[,integration_cluster],"\n",
                                                         "Cell: ",row.names(data@meta.data),"\n",
                                                         data@meta.data[,"INTEGRATED_CELL_TYPE"], sep = ""))

    htmlwidgets::saveWidget(p, paste(cdir, "96URSA_PLOT_scRNASEQ_UMAP_INTERACTIVE_CELL_TYPES_",integration_name,"_",project_name,".html", sep = ""))

    p <- NULL
    p <- plot_3d(data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                     data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                                     data@reductions$tsne@cell.embeddings[,"tSNE_3"],
                                     color_group = data@meta.data[,"INTEGRATED_CELL_TYPE"],
                                     plot_title = names(data_current)[j],
                                     col = color_conditions$tenx,
                                     n = length(unique(data$INTEGRATED_CELL_TYPE)),
                                     lt = "CELL TYPE", t1 = "tSNE_1",t2 = "tSNE_2", t3= "tSNE_3",
                                     current_text= paste("Cluster ",data@meta.data[,integration_cluster],"\n",
                                                         "Cell: ",row.names(data@meta.data),"\n",
                                                         data@meta.data[,"INTEGRATED_CELL_TYPE"], sep = ""))
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
