############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Dubhe: scImmune
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
#' @import ggraph
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
#' @import circlize
#' @import ggalluvial
#' @import harmony
#' @import immunarch
#' @import scRepertoire
#' @import Seurat
#' @import SingleR
#' @import vegan
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format scImmune Profiling pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification scImmune
#' Profiling pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_scImmune' by default.
#' @param input_dir Directory to all input files. Current working directory by
#' default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format
#' files.
#' @param integration_method Integration method for combining scRNASeq data.
#' Default by harmony.
#' @param ... Arguments passed to other parameters in the dependency pages.
#' Parameters with the long format: xxx.xxx.xxx usually indicates in lowercases
#' the parameter origin: <dependent package name>.<function name>.<parameter name>(),
#' for example: immunarch.repLoad.mode() indicates this parameter 
#' originates from the package immunarch under the function repLoad() and its
#' parameter 'mode'. Users could subsequently refer to the dependency 
#' package for detailed information on parameters and their usage.
#' @export
scImmunePip <- function(project_name = "Ursa_scImmune",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file,
                     integration_method = c("harmony","seurat"),
                     
                     # repLoad
                     immunarch.repLoad.mode = "paired",
                     immunarch.repLoad.coding = TRUE,
                     
                     # repOverlap - method 1
                     immunarch.repoverlap.method1.method = "morisita",
                     immunarch.repoverlap.method1.col = "aa",
                     immunarch.repoverlap.method1.a = 0.5,
                     immunarch.repoverlap.method1.b = 0.5,
                     immunarch.repoverlap.method1.step = 1000,
                     immunarch.repoverlap.method1.n.steps = 10,
                     immunarch.repoverlap.method1.downsample = FALSE,
                     immunarch.repoverlap.method1.bootstrap = NA,
                     immunarch.repoverlap.method1.force.matrix = FALSE,
                     
                     # repOverlap - method 2
                     immunarch.repoverlap.method2.method = "morisita",
                     immunarch.repoverlap.method2.col = "aa",
                     immunarch.repoverlap.method2.a = 0.5,
                     immunarch.repoverlap.method2.b = 0.5,
                     immunarch.repoverlap.method2.step = 1000,
                     immunarch.repoverlap.method2.n.steps = 10,
                     immunarch.repoverlap.method2.downsample = FALSE,
                     immunarch.repoverlap.method2.bootstrap = NA,
                     immunarch.repoverlap.method2.force.matrix = FALSE,
                     
                     # repOverlap - method 3
                     immunarch.repoverlap.method3.method = "morisita",
                     immunarch.repoverlap.method3.a = 0.5,
                     immunarch.repoverlap.method3.b = 0.5,
                     immunarch.repoverlap.method3.step = 1000,
                     immunarch.repoverlap.method3.n.steps = 10,
                     immunarch.repoverlap.method3.downsample = FALSE,
                     immunarch.repoverlap.method3.bootstrap = NA,
                     immunarch.repoverlap.method3.force.matrix = FALSE,
                     
                     # geneUsage
                     immunarch.geneusage.gene = c("hs.trbv", "HomoSapiens.TRBJ", "macmul.IGHV"),
                     immunarch.geneusage.quant = c(NA, "count"),
                     immunarch.geneusage.ambig = c("inc", "exc", "maj"),
                     immunarch.geneusage.type = c("segment", "allele", "family"),
                     immunarch.geneusage.norm = FALSE,
                     
                     # geneUsageAnalysis - method 1
                     immunarch.geneusageanalysis.method1.method = c("js+hclust", "pca+kmeans", "anova", "js+pca+kmeans"),
                     immunarch.geneusageanalysis.method1.base = 2,
                     immunarch.geneusageanalysis.method1.norm.entropy = FALSE,
                     immunarch.geneusageanalysis.method1.cor = c("pearson", "kendall", "spearman"),
                     immunarch.geneusageanalysis.method1.do.norm = TRUE,
                     immunarch.geneusageanalysis.method1.laplace = 1e-12,
                     immunarch.geneusageanalysis.method1.k = 2,
                     immunarch.geneusageanalysis.method1.eps = 0.01,
                     immunarch.geneusageanalysis.method1.perp = 1,
                     immunarch.geneusageanalysis.method1.theta = 0.1,
                     
                     # geneUsageAnalysis - method 2
                     immunarch.geneusageanalysis.method2.method = c("js+hclust", "pca+kmeans", "anova", "js+pca+kmeans"),
                     immunarch.geneusageanalysis.method2.base = 2,
                     immunarch.geneusageanalysis.method2.norm.entropy = FALSE,
                     immunarch.geneusageanalysis.method2.cor = c("pearson", "kendall", "spearman"),
                     immunarch.geneusageanalysis.method2.do.norm = TRUE,
                     immunarch.geneusageanalysis.method2.laplace = 1e-12,
                     immunarch.geneusageanalysis.method2.k = 2,
                     immunarch.geneusageanalysis.method2.eps = 0.01,
                     immunarch.geneusageanalysis.method2.perp = 1,
                     immunarch.geneusageanalysis.method2.theta = 0.1,
                     
                     # geneUsageAnalysis - method 3
                     immunarch.geneusageanalysis.method3.method = c("js+hclust", "pca+kmeans", "anova", "js+pca+kmeans"),
                     immunarch.geneusageanalysis.method3.base = 2,
                     immunarch.geneusageanalysis.method3.norm.entropy = FALSE,
                     immunarch.geneusageanalysis.method3.cor = c("pearson", "kendall", "spearman"),
                     immunarch.geneusageanalysis.method3.do.norm = TRUE,
                     immunarch.geneusageanalysis.method3.laplace = 1e-12,
                     immunarch.geneusageanalysis.method3.k = 2,
                     immunarch.geneusageanalysis.method3.eps = 0.01,
                     immunarch.geneusageanalysis.method3.perp = 1,
                     immunarch.geneusageanalysis.method3.theta = 0.1,
                     
                     # spectratype - method 1
                     immunarch.spectratype.method1.quant = "id",
                     immunarch.spectratype.method1.col = "nt",
                     
                     # spectratype - method 2
                     immunarch.spectratype.method2.quant = "count",
                     immunarch.spectratype.method2.col = "aa+v",
                     
                     # repDiversity
                     immunarch.repdiversity.col = "aa",
                     immunarch.repdiversity.max.q = 6,
                     immunarch.repdiversity.min.q = 1,
                     immunarch.repdiversity.q = 5,
                     immunarch.repdiversity.step = NA,
                     immunarch.repdiversity.quantile = c(0.025, 0.975),
                     immunarch.repdiversity.extrapolation = NA,
                     immunarch.repdiversity.perc = 50,
                     immunarch.repdiversity.norm = TRUE,
                     immunarch.repdiversity.do.norm = NA,
                     immunarch.repdiversity.laplace = 0,
                     
                     # trackClonotypes - method 1
                     immunarch.trackclonotypes.method1.col = "aa",
                     immunarch.trackclonotypes.method1.norm = TRUE,
                     
                     # trackClonotypes - method 2
                     immunarch.trackclonotypes.method2.col = "aa",
                     immunarch.trackclonotypes.method2.norm = TRUE,
                     
                     # combineBCR
                     screpertoire.combinebcr.call.related.clones = TRUE,
                     screpertoire.combinebcr.threshold = 0.85,
                     screpertoire.combinebcr.removeNA = TRUE,
                     screpertoire.combinebcr.removeMulti = FALSE,
                     screpertoire.combinebcr.filterMulti = TRUE,
                     
                     # combineTCR
                     screpertoire.combinetcr.removeNA = TRUE,
                     screpertoire.combinetcr.removeMulti = FALSE,
                     screpertoire.combinetcr.filterMulti = TRUE,
                     
                     # quantContig
                     screpertoire.quantcontig.cloneCall = "strict",
                     screpertoire.quantcontig.chain = "both",
                     screpertoire.quantcontig.scale = TRUE,
                     screpertoire.quantcontig.group.by = NULL,
                     screpertoire.quantcontig.split.by = NULL,
                     screpertoire.quantcontig.order = TRUE,
                     screpertoire.quantcontig.exportTable = TRUE,
                     
                     # abundanceContig
                     screpertoire.abundancecontig.cloneCall = "strict",
                     screpertoire.abundancecontig.chain = "both",
                     screpertoire.abundancecontig.group.by = NULL,
                     screpertoire.abundancecontig.split.by = NULL,
                     screpertoire.abundancecontig.order = TRUE,
                     
                     # lengthContig - choice 1
                     screpertoire.lengthcontig.choice1.cloneCall = "aa",
                     screpertoire.lengthcontig.choice1.chain = "both",
                     screpertoire.lengthcontig.choice1.group.by = NULL,
                     screpertoire.lengthcontig.choice1.split.by = NULL,
                     screpertoire.lengthcontig.choice1.order = TRUE,
                     screpertoire.lengthcontig.choice1.scale = FALSE,
                     screpertoire.lengthcontig.choice1.exportTable = FALSE,
                     
                     # lengthContig - choice 2
                     screpertoire.lengthcontig.choice2.cloneCall = "aa",
                     screpertoire.lengthcontig.choice2.chain = "both",
                     screpertoire.lengthcontig.choice2.group.by = NULL,
                     screpertoire.lengthcontig.choice2.split.by = NULL,
                     screpertoire.lengthcontig.choice2.order = TRUE,
                     screpertoire.lengthcontig.choice2.scale = FALSE,
                     screpertoire.lengthcontig.choice2.exportTable = FALSE,
                     
                     # compareClonotypes - choice 1
                     screpertoire.compareclonotypes.choice1.cloneCall = "strict",
                     screpertoire.compareclonotypes.choice1.chain = "both",
                     screpertoire.compareclonotypes.choice1.samples = NULL,
                     screpertoire.compareclonotypes.choice1.numbers = NULL,
                     screpertoire.compareclonotypes.choice1.split.by = NULL,
                     screpertoire.compareclonotypes.choice1.graph = "alluvial",
                     screpertoire.compareclonotypes.choice1.exportTable = FALSE,
                     
                     # compareClonotypes - choice 2
                     screpertoire.compareclonotypes.choice2.cloneCall = "aa",
                     screpertoire.compareclonotypes.choice2.chain = "both",
                     screpertoire.compareclonotypes.choice2.samples = NULL,
                     screpertoire.compareclonotypes.choice2.numbers = NULL,
                     screpertoire.compareclonotypes.choice2.split.by = NULL,
                     screpertoire.compareclonotypes.choice2.graph = "alluvial",
                     screpertoire.compareclonotypes.choice2.exportTable = FALSE,
                     
                     # compareClonotypes - choice 3
                     screpertoire.compareclonotypes.choice3.cloneCall = "nt",
                     screpertoire.compareclonotypes.choice3.chain = "both",
                     screpertoire.compareclonotypes.choice3.samples = NULL,
                     screpertoire.compareclonotypes.choice3.numbers = NULL,
                     screpertoire.compareclonotypes.choice3.split.by = NULL,
                     screpertoire.compareclonotypes.choice3.graph = "alluvial",
                     screpertoire.compareclonotypes.choice3.exportTable = FALSE,
                     
                     # compareClonotypes - choice 4
                     screpertoire.compareclonotypes.choice4.cloneCall = "gene",
                     screpertoire.compareclonotypes.choice4.chain = "both",
                     screpertoire.compareclonotypes.choice4.samples = NULL,
                     screpertoire.compareclonotypes.choice4.numbers = NULL,
                     screpertoire.compareclonotypes.choice4.split.by = NULL,
                     screpertoire.compareclonotypes.choice4.graph = "alluvial",
                     screpertoire.compareclonotypes.choice4.exportTable = FALSE,
                     
                     # vizGenes
                     screpertoire.vizgenes.gene = "V",
                     screpertoire.vizgenes.chain = "TRA",
                     screpertoire.vizgenes.plot = "heatmap",
                     screpertoire.vizgenes.y.axis = "sample",
                     screpertoire.vizgenes.order = "gene",
                     screpertoire.vizgenes.scale = TRUE,
                     screpertoire.vizgenes.group.by = NULL,
                     screpertoire.vizgenes.split.by = NULL,
                     screpertoire.vizgenes.exportTable = FALSE,
                     
                     # clonalHomeostasis
                     screpertoire.clonalhomeostasis.cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1),
                     screpertoire.clonalhomeostasis.cloneCall = "strict",
                     screpertoire.clonalhomeostasis.chain = "both",
                     screpertoire.clonalhomeostasis.group.by = NULL,
                     screpertoire.clonalhomeostasis.split.by = NULL,
                     screpertoire.clonalhomeostasis.exportTable = FALSE,
                     
                     # clonalOverlap
                     screpertoire.clonaloverlap.cloneCall = "strict",
                     screpertoire.clonaloverlap.method = c("morisita"),
                     screpertoire.clonaloverlap.chain = "both",
                     screpertoire.clonaloverlap.split.by = NULL,
                     screpertoire.clonaloverlap.exportTable = FALSE,
                     
                     # clonalProportion
                     screpertoire.clonalproportion.split = c(10, 100, 1000, 10000, 30000, 1e+05),
                     screpertoire.clonalproportion.cloneCall = "strict",
                     screpertoire.clonalproportion.chain = "both",
                     screpertoire.clonalproportion.group.by = NULL,
                     screpertoire.clonalproportion.split.by = NULL,
                     screpertoire.clonalproportion.exportTable = FALSE,
                     
                     # clonesizeDistribution
                     screpertoire.clonesizedistribution.cloneCall = "strict",
                     screpertoire.clonesizedistribution.chain = "both",
                     screpertoire.clonesizedistribution.method = "ward.D2",
                     screpertoire.clonesizedistribution.threshold = 1,
                     screpertoire.clonesizedistribution.group.by = NULL,
                     screpertoire.clonesizedistribution.split.by = NULL,
                     screpertoire.clonesizedistribution.exportTable = FALSE,
                     
                     # hclust
                     screpertoire.clonesizedistribution.hclust.method = "ward.D2",
                     
                     # clonalDiversity
                     screpertoire.clonesizediversity.chain = "both",
                     screpertoire.clonesizediversity.group.by = NULL,
                     screpertoire.clonesizediversity.x.axis = NULL,
                     screpertoire.clonesizediversity.split.by = NULL,
                     screpertoire.clonesizediversity.exportTable = FALSE,
                     screpertoire.clonesizediversity.n.boots = 100,
                     screpertoire.clonesizediversity.return.boots = FALSE,
                     screpertoire.clonesizediversity.skip.boots = FALSE,
                     
                     # CreateSeuratObject
                     seurat.createseuratobject.min.cells = 3,
                     
                     # subset - scRNA-seq
                     seurat.subset.min.nFeature_RNA = 200,
                     seurat.subset.max.nFeature_RNA = 6000,
                     seurat.subset.max.percent.mt = 10,
                     
                     # NormalizeData
                     seurat.normalizedata.assay = NULL,
                     seurat.normalizedata.normalization.method = "LogNormalize",
                     seurat.normalizedata.scale.factor = 10000,
                     seurat.normalizedata.margin = 1,
                     
                     # FindVariableFeatures
                     seurat.findvariablefeatures.selection.method = "vst",
                     seurat.findvariablefeatures.loess.span = 0.3,
                     seurat.findvariablefeatures.clip.max = "auto",
                     # seurat.findvariablefeatures.mean.function = FastExpMean,
                     # seurat.findvariablefeatures.dispersion.function = FastLogVMR,
                     seurat.findvariablefeatures.num.bin = 20,
                     seurat.findvariablefeatures.binning.method = "equal_width",
                     seurat.findvariablefeatures.nfeatures = 2000,
                     seurat.findvariablefeatures.mean.cutoff = c(0.1, 8),
                     seurat.findvariablefeatures.dispersion.cutoff = c(1, Inf),
                     
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
                     
                     # Seurat - SelectIntegrationFeatures
                     seurat.selectintegrationfeatures.nfeatures = 2000,
                     seurat.selectintegrationfeatures.assay = NULL,
                     seurat.selectintegrationfeatures.fvf.nfeatures = 2000,
                     
                     # Seurat - FindIntegrationAnchors
                     seurat.findintegrationanchors.assay = NULL,
                     seurat.findintegrationanchors.reference = NULL,
                     seurat.findintegrationanchors.scale = TRUE,
                     seurat.findintegrationanchors.normalization.method = "LogNormalize",
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
                     seurat.integratedata.new.assay.name = "integrated",
                     seurat.integratedata.normalization.method = "LogNormalize",
                     seurat.integratedata.features = NULL,
                     seurat.integratedata.features.to.integrate = NULL,
                     seurat.integratedata.dims = 30,
                     seurat.integratedata.k.weight = 100,
                     seurat.integratedata.weight.reduction = NULL,
                     seurat.integratedata.sd.weight = 1,
                     seurat.integratedata.sample.tree = NULL,
                     seurat.integratedata.preserve.order = FALSE,
                     seurat.integratedata.eps = 0,
                     
                     # harmony - RunHarmony
                     harmony.runharmony.dims = NULL,
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
                     
                     # Seurat - Integration FindVariableFeatures
                     seurat.integration.findvariablefeatures.nfeatures = 2000,
                     seurat.integration.findvariablefeatures.selection.method = "vst",
                     seurat.integration.findvariablefeatures.loess.span = 0.3,
                     seurat.integration.findvariablefeatures.clip.max = "auto",
                     # seurat.integration.findvariablefeatures.mean.function = FastExpMean,
                     # seurat.integration.findvariablefeatures.dispersion.function = FastLogVMR,
                     seurat.integration.findvariablefeatures.num.bin = 20,
                     seurat.integration.findvariablefeatures.binning.method = "equal_width",
                     
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
                     # seurat.integration.runumap.assay = DefaultAssay(object),
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
                     
                     # Seurat - FindClusters
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

                     # SingleR
                     singler.ref = "HumanPrimaryCellAtlas",
                     singler.labels = "HumanPrimaryCellAtlasLevel1",
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
                     
                     # combineExpression
                     screpertoire.combineexpression.cloneCall = "strict",
                     screpertoire.combineexpression.chain = "both",
                     screpertoire.combineexpression.group.by = "none",
                     screpertoire.combineexpression.proportion = TRUE,
                     screpertoire.combineexpression.filterNA = FALSE,
                     screpertoire.combineexpression.cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1),
                     screpertoire.combineexpression.addLabel = FALSE,
                     
                     # occupiedscRepertoire
                     screpertoire.occupiedscrepertoire.label = TRUE,
                     screpertoire.occupiedscrepertoire.facet.by = NULL,
                     screpertoire.occupiedscrepertoire.proportion = FALSE,
                     screpertoire.occupiedscrepertoire.na.include = FALSE,
                     screpertoire.occupiedscrepertoire.exportTable = FALSE,
                     
                     # clonalDiversity
                     screpertoire.clonaldiversity.chain = "both",
                     screpertoire.clonaldiversity.group.by = NULL,
                     screpertoire.clonaldiversity.x.axis = NULL,
                     screpertoire.clonaldiversity.split.by = NULL,
                     screpertoire.clonaldiversity.exportTable = FALSE,
                     screpertoire.clonaldiversity.n.boots = 100,
                     screpertoire.clonaldiversity.return.boots = FALSE,
                     screpertoire.clonaldiversity.skip.boots = FALSE,
                     
                     # clonalOverlay
                     screpertoire.clonaloverlay.freq.cutpoint = 0,
                     screpertoire.clonaloverlay.bins = 10,
                     screpertoire.clonaloverlay.facet = "GROUP",
                     
                     # clonalNetwork
                     screpertoire.clonalnetwork.identity = "seurat_clusters",
                     screpertoire.clonalnetwork.filter.clones = NULL,
                     screpertoire.clonalnetwork.filter.identity = NULL,
                     screpertoire.clonalnetwork.filter.proportion = NULL,
                     screpertoire.clonalnetwork.filter.graph = FALSE,
                     screpertoire.clonalnetwork.cloneCall = "aa",
                     screpertoire.clonalnetwork.chain = "both",
                     screpertoire.clonalnetwork.exportTable = FALSE,
                     screpertoire.clonalnetwork.exportClones = FALSE,
                     
                     # getKmers
                     immunarch.getkmers.col = c("aa", "nt"),
                     immunarch.getkmers.coding = TRUE,
                     
                     # kmer_profile
                     screpertoire.kmer_profile..method = c("freq"),
                     screpertoire.kmer_profile..remove.stop = TRUE){
  
  print("Initialising pipeline environment..")
  integration_method <- integration_method[1]
  pheno_data <- pheno_ini(pheno_file, pipeline = "scIMMUNE", isDir = T)
  pheno_data$SID <- paste(pheno_data$SAMPLE_ID, pheno_data$GROUP, pheno_data$CELL_TYPE, sep = "_")
  color_conditions <- color_ini()
  ctime <- time_ini()
  sample_colors <- gen_colors(color_conditions$manycolors, length(unique(pheno_data$SID)))
  names(sample_colors) <- unique(pheno_data$SID)
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir, recursive = T)
  
  hpca.se <- HumanPrimaryCellAtlasData()
  if(toupper(singler.ref) == toupper("HumanPrimaryCellAtlas")){
    singler.ref <- hpca.se
    if(toupper(singler.labels)[1] == toupper("HumanPrimaryCellAtlasLevel1")){
      singler.labels <- hpca.se$label.main
    }else if(toupper(singler.labels)[1] == toupper("HumanPrimaryCellAtlasLevel2")){
      singler.labels <- hpca.se$label.fine
    }
  }
  contig_dir <- paste(cdir, "/Contig/", sep = "")
  dir.create(contig_dir, recursive = T)
  
  sample_files <- list.files(input_dir, recursive = T, full.names = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep(gsub(".*\\/(.*)","\\1",pheno_file, ignore.case = T),sample_files, ignore.case = T, invert = T)]
  input_dir <- gsub("(.*\\/).*$","\\1",sample_files[1], ignore.case = T)
  vdjdb <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb")

  contig_pheno <- pheno_data[grep("contig", pheno_data$DATA_TYPE, ignore.case = T),]
  consensus_pheno <- pheno_data[grep("consensus", pheno_data$DATA_TYPE, ignore.case = T),]
  rna_pheno <- pheno_data[grep("rna", pheno_data$DATA_TYPE, ignore.case = T),]
  #######################################################################################################################################
  print("Reading data..")

  contig_list <- NULL
  for(i in 1:nrow(contig_pheno)){
    contig_list[[i]] <- sample_files[grep(contig_pheno$FILE[i],sample_files, ignore.case = T)]
    system(paste("cp ", contig_list[[i]], " ",contig_dir,"/", sep = ""))
    contig_list[[i]] <- read.csv(paste(contig_dir,"/",contig_pheno$FILE[i],sep = ""), header = T, sep = ",")
    names(contig_list)[i] <- contig_pheno[i,"SID"]
  }

  contig_underscore <- as.numeric(lapply(contig_list, function(x){ x<- length(grep("_", x$barcode))}))
  contig_hyphen <- as.numeric(lapply(contig_list, function(x){ x<- length(grep("-", x$barcode))}))

  if(length(which(contig_underscore == 0))>1 & length(which(contig_hyphen > 0)) == length(contig_list)){
    contig_list <- contig_list
  } else{
    for (i in seq_along(contig_list)) {
      contig_list[[i]] <- stripBarcode(contig_list[[i]], column = 1, connector = "_",
                                       num_connects = max(as.numeric(lapply(contig_list, function(x){ x<- length(grep("_", x$barcode))})))+1)
    }
  }

  consensus_list <- NULL
  # consensus_pheno <- pheno_data[grep("consensus", pheno_data$DATA_TYPE, ignore.case = T),]
  for(i in 1:nrow(consensus_pheno)){
    consensus_list[[i]] <- sample_files[grep(consensus_pheno$FILE[i],sample_files, ignore.case = T)]
    consensus_list[[i]] <- read.csv(paste(input_dir,"/",consensus_pheno$FILE[i],sep = ""), header = T, sep = ",")
    names(consensus_list)[i] <- consensus_pheno[i,"SID"]
  }

  consensus_underscore <- as.numeric(lapply(consensus_list, function(x){ x<- length(grep("_", x$barcode))}))
  consensus_hyphen <- as.numeric(lapply(consensus_list, function(x){ x<- length(grep("-", x$barcode))}))

  if(length(which(contig_underscore == 0))>1 & length(which(contig_hyphen > 0)) == length(consensus_list)){
    consensus_list <- consensus_list
  } else{

    for (i in seq_along(consensus_list)) {
      consensus_list[[i]] <- stripBarcode(consensus_list[[i]], column = 1, connector = "_",
                                          num_connects = max(as.numeric(lapply(consensus_list, function(x){ x<- length(grep("_", x$barcode))})))+1)
    }

  }

  if(nrow(rna_pheno) > 0){
    rna_list <- NULL
    for(i in 1:nrow(rna_pheno)){
      rna_list[[i]] <- sample_files[grep(rna_pheno$FILE[i],sample_files, ignore.case = T)]
    }
    rna_list <- rna_list[grep("_MACOSX",rna_list, ignore.case = T, invert = T)]
  }

  contig_monarch <- repLoad(contig_dir,
                            .mode = immunarch.repLoad.mode,
                            .coding = immunarch.repLoad.coding)
  for(i in 1:length(contig_monarch$data)){
    names(contig_monarch$data)[i] <- pheno_data[which(unlist(lapply(gsub("(.*)\\..*","\\1",pheno_data$FILE, ignore.case = T), function(x){grepl(x,names(contig_monarch$data)[i], ignore.case = T)})) == TRUE),"SID"]
  }
  contig_monarch$meta <- pheno_data[match(names(contig_monarch$data), pheno_data$SID),]
  contig_monarch$meta$Sample <- pheno_data[match(names(contig_monarch$data), pheno_data$SID),"SID"]

  meta_explore <- repExplore(contig_monarch$data, .method = "volume")
  p1 <- vis(meta_explore, .by = c("Sample"), .meta = contig_monarch$meta, .test = F) + xlab("Sample")
  p2 <- vis(meta_explore, .by = c("GROUP", "CELL_TYPE"), .meta = contig_monarch$meta, .test = F)

  p1plots <- p1 + p2

  exp_len_aa <- repExplore(contig_monarch$data, .method = "len", .col = "aa")
  exp_len_nt <- repExplore(contig_monarch$data, .method = "len", .col = "nt")
  exp_cnt <- repExplore(contig_monarch$data, .method = "count")
  exp_vol <- repExplore(contig_monarch$data, .method = "volume")

  p1 <- vis(exp_len_aa) + facet_wrap(~Sample) + ggtitle("DISTRIBUTION OF CDR3 LENGTH: AMINO ACID")
  p11 <- vis(exp_len_nt) + facet_wrap(~Sample) + ggtitle("DISTRIBUTION OF CDR3 LENGTH: NUCLEOTIDE")
  p2 <- vis(exp_cnt) + guides(color = guide_legend(ncol = ifelse(nrow(contig_monarch$meta) > 10, 2, 1)))
  p3 <- vis(exp_vol)

  p2plots <- p1 + p11
  p3plots <- p2 + p3

  imm_pr <- repClonality(contig_monarch$data, .method = "clonal.prop")
  imm_top <- repClonality(contig_monarch$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
  imm_rare <- repClonality(contig_monarch$data, .method = "rare")
  imm_hom <- repClonality(contig_monarch$data, .method = "homeo",.clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))

  p4plots <- vis(imm_top) + vis(imm_top, .by = "Sample", .meta = contig_monarch$meta)
  p5plots <- vis(imm_rare) + vis(imm_rare, .by = "Sample", .meta = contig_monarch$meta)
  p6plots <- vis(imm_hom) + vis(imm_hom, .by = c("Sample","GROUP"), .meta = contig_monarch$meta)

  imm_ov1 <- repOverlap(contig_monarch$data,
                        .method = immunarch.repoverlap.method1.method, # "morisita"
                        .col = immunarch.repoverlap.method1.col, # "nt+v+d+j"
                        .a = immunarch.repoverlap.method1.a,
                        .b = immunarch.repoverlap.method1.b,
                        .step = immunarch.repoverlap.method1.step,
                        .n.steps = immunarch.repoverlap.method1.n.steps,
                        .downsample = immunarch.repoverlap.method1.downsample,
                        .bootstrap = immunarch.repoverlap.method1.bootstrap,
                        .force.matrix = immunarch.repoverlap.method1.force.matrix)
  imm_ov2 <- repOverlap(contig_monarch$data,
                        .method = immunarch.repoverlap.method2.method, # "morisita"
                        .col = immunarch.repoverlap.method2.col, # "aa"
                        .a = immunarch.repoverlap.method2.a,
                        .b = immunarch.repoverlap.method2.b,
                        .step = immunarch.repoverlap.method2.step,
                        .n.steps = immunarch.repoverlap.method2.n.steps,
                        .downsample = immunarch.repoverlap.method2.downsample,
                        .bootstrap = immunarch.repoverlap.method2.bootstrap,
                        .force.matrix = immunarch.repoverlap.method2.force.matrix)
  

  p1 <- vis(imm_ov1, .text.size = 6) + ggtitle(paste("CLONOTYPES SHARED\n",stringr::str_to_title(immunarch.repoverlap.method1.method)," OVERLAP INDEX (",immunarch.repoverlap.method1.col,")", sep = ""))
  p2 <- vis(imm_ov2, .text.size = 6) + ggtitle(paste("CLONOTYPES SHARED\n",stringr::str_to_title(immunarch.repoverlap.method2.method)," (",immunarch.repoverlap.method2.col,")", sep = ""))
  p7plots  <- p1 + p2

  # Build public repertoire table using CDR3 nucleotide sequences
  gene_stats <- gene_stats()
  gene_stats <- gene_stats[grep("hs", gene_stats$alias,ignore.case = T),]
  all_genes <- colnames(gene_stats)[(!grepl("alias|species", colnames(gene_stats), ignore.case = T)) &
                                      (gene_stats != 0)]
  all_gene_usage <- NULL
  for(i in 1:length(all_genes)){
    all_gene_usage <- rbind(all_gene_usage, data.frame(Gene_Type = all_genes[i], geneUsage(contig_monarch$data, .gene = paste("hs.",all_genes[i], sep = ""))))
  }

  n <- 20
  groups <- unique(contig_monarch$meta$CELL_TYPE)
  names <- NULL
  p9plots <- NULL

  for(i in 1:length(groups)){
    for(j in 1:length(all_genes)){
      current_gu <- geneUsage(contig_monarch$data[which(toupper(names(contig_monarch$data)) %in% toupper(unlist(contig_monarch$meta[contig_monarch$meta$CELL_TYPE == groups[i],"Sample"])))], .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
      current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA|None",x, ignore.case = T))) > 0})),]
      current_gu[is.na(current_gu)] <- 0
      current_gu$Names <- gsub("NA;|None;","",current_gu$Names)
      current_gu$Names <- gsub(";NA|;None","",current_gu$Names)
      columnnames <- colnames(current_gu)
      current_gu <- split(current_gu, current_gu$Names)
      current_gu <- lapply(current_gu, function(x){
        if(nrow(x) > 1){
          x <- data.frame(Names = unique(x$Names), t(colSums(x[,grep("^Names$", colnames(x), invert = T)])))
        }else{
          x <- x
        }
      })
      current_gu <- do.call(rbind.data.frame, current_gu)
      current_gu <- as_tibble(current_gu)
      class(current_gu) <- c("immunr_gene_usage","immunr_gene_usage","tbl_df","tbl","data.frame")
      if(nrow(current_gu) > n){
        current_gu <- current_gu[order(rowMeans(current_gu[,grep("Name", colnames(current_gu), ignore.case = T, invert = T)]), decreasing = T),]
        current_gu <- current_gu[1:n,]
        current_gu <- current_gu[!is.na(current_gu$Names),]
      }else{
        n <- nrow(current_gu)
      }
      p1 <- NULL
      p1 <- vis(current_gu) + ggtitle(paste("TOP",n,": GENE USAGE: ", groups[i], ", ", toupper(all_genes[j]),
                                            "\n(All results will be plotted for those with total of <",n," hits)",sep = ""))
      if(i == 1 & j == 1){
        p9plots[[1]] <- p1
        names(p9plots)[1] <- paste(groups[i], "_", toupper(all_genes[j]), sep = "")
      }else{
        p9plots[[length(p9plots)+1]] <- p1
        names(p9plots)[length(p9plots)] <- paste(groups[i], "_", toupper(all_genes[j]), sep = "")
      }
      n <- 20
    }
  }

  p10plots <- NULL
  p11plots <- NULL
  p12plots <- NULL
  
  for(j in 1:length(all_genes)){
    current_gu <- geneUsage(contig_monarch$data,
                            .gene = paste("hs.", all_genes[j], sep = ""),
                            .quant = immunarch.geneusage.quant,
                            .ambig = immunarch.geneusage.ambig,
                            .type = immunarch.geneusage.type,
                            .norm = immunarch.geneusage.norm)
    current_gu <- current_gu[!is.na(current_gu$Names),]
    current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
    current_gu[is.na(current_gu)] <- 0
    imm_gu_m1 <- geneUsageAnalysis(current_gu,
                                   .method = immunarch.geneusageanalysis.method1.method, # "js"
                                   .base = immunarch.geneusageanalysis.method1.base,
                                   .norm.entropy = immunarch.geneusageanalysis.method1.norm.entropy,
                                   .cor = immunarch.geneusageanalysis.method1.cor,
                                   .do.norm = immunarch.geneusageanalysis.method1.do.norm,
                                   .laplace = immunarch.geneusageanalysis.method1.laplace,
                                   .k = immunarch.geneusageanalysis.method1.k,
                                   .eps = immunarch.geneusageanalysis.method1.eps,
                                   .perp = immunarch.geneusageanalysis.method1.perp,
                                   .theta = immunarch.geneusageanalysis.method1.theta)
    imm_gu_m2 <- geneUsageAnalysis(current_gu,
                                   .method = immunarch.geneusageanalysis.method2.method, # "cor"
                                   .base = immunarch.geneusageanalysis.method2.base,
                                   .norm.entropy = immunarch.geneusageanalysis.method2.norm.entropy,
                                   .cor = immunarch.geneusageanalysis.method2.cor,
                                   .do.norm = immunarch.geneusageanalysis.method2.do.norm,
                                   .laplace = immunarch.geneusageanalysis.method2.laplace,
                                   .k = immunarch.geneusageanalysis.method2.k,
                                   .eps = immunarch.geneusageanalysis.method2.eps,
                                   .perp = immunarch.geneusageanalysis.method2.perp,
                                   .theta = immunarch.geneusageanalysis.method2.theta)
    p1 <- vis(imm_gu_m1, .title = paste("Gene usage ",toupper(immunarch.geneusageanalysis.method1.method),": ", toupper(all_genes[j]), sep = ""), .leg.title = toupper(immunarch.geneusageanalysis.method1.method), .text.size = 4)
    p2 <- vis(imm_gu_m2, .title = paste("Gene usage ",toupper(immunarch.geneusageanalysis.method2.method),": ", toupper(all_genes[j]), sep = ""), .leg.title = toupper(immunarch.geneusageanalysis.method2.method), .text.size = 4)
    p <- p1 + p2
    p10plots[[j]] <- p
    names(p10plots)[j] <- all_genes[j]
    
    row.names(current_gu) <- current_gu$Names
    p <- NULL
    p <- vis(geneUsageAnalysis(current_gu,
                               .method = immunarch.geneusageanalysis.method3.method, # "cosine+hclust",
                               .base = immunarch.geneusageanalysis.method3.base,
                               .norm.entropy = immunarch.geneusageanalysis.method3.norm.entropy,
                               .cor = immunarch.geneusageanalysis.method3.cor,
                               .do.norm = immunarch.geneusageanalysis.method3.do.norm,
                               .laplace = immunarch.geneusageanalysis.method3.laplace,
                               .k = immunarch.geneusageanalysis.method3.k,
                               .eps = immunarch.geneusageanalysis.method3.eps,
                               .perp = immunarch.geneusageanalysis.method3.perp,
                               .theta = immunarch.geneusageanalysis.method3.theta), .rect = T)+ggtitle(paste("Gene usage ",toupper(immunarch.geneusageanalysis.method3.method),": ", toupper(all_genes[j]), sep = ""))
    p11plots[[j]] <- p
    names(p11plots)[j] <- all_genes[j]
    
    if(length(contig_monarch$data) > 3){
      imm_cl_pca <- geneUsageAnalysis(current_gu, "js+pca+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
      imm_cl_mds <- geneUsageAnalysis(current_gu, "js+mds+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .verbose = F)
      # imm_cl_tsne <- geneUsageAnalysis(current_gu, "js+tsne+kmeans", .k = length(unique(contig_monarch$meta$CELL_TYPE)), .perp = .01, .verbose = F)
      p1 <- vis(imm_cl_pca, .plot = "clust") + ggtitle(paste("PCA: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+PCA+K-MEANS)", sep = ""))
      p2 <- vis(imm_cl_mds, .plot = "clust") + ggtitle(paste("MDS: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+MDS+K-MEANS)", sep = ""))
      # p3 <- vis(imm_cl_tsne, .plot = "clust") + ggtitle(paste("tSNE: ", toupper(all_genes[j]), "\n(JS-DIVERGENCE+tSNE+K-MEANS)", sep = ""))
      p12plots[[j]] <- p1+p2
      names(p12plots)[j] <- all_genes[j]
    }
  }

  p13plots <- NULL
  for(i in 1:length(contig_monarch$data)){
    p1 <- vis(spectratype(contig_monarch$data[[i]],
                          .quant = immunarch.spectratype.method1.quant, # "id"
                          .col = immunarch.spectratype.method1.col))+ # "nt"
      ggtitle(paste("SAMPLE: ",names(contig_monarch$data)[i],sep = "")) +
      xlab(toupper(immunarch.spectratype.method1.col))
    p2 <- vis(spectratype(contig_monarch$data[[i]],
                          .quant = immunarch.spectratype.method2.quant, # "count"
                          .col = immunarch.spectratype.method2.col)) + # "aa+v"
      xlab(toupper(immunarch.spectratype.method2.col))
    p <- p1 + p2
    p13plots[[i]] <- p
    names(p13plots)[i] <- names(contig_monarch$data)[i]
  }

  # Diversity estimation
  diversity_estimates <- NULL
  diversity_methods <- c("chao1", "hill", "div", "gini.simp", "inv.simp", "d50","raref")
  p14plots <- NULL

  for(i in 1:length(diversity_methods)){
    p14plots[[i]] <- vis(repDiversity(contig_monarch$data,
                                      .method = diversity_methods[i],
                                      .col = immunarch.repdiversity.col,
                                      .max.q = immunarch.repdiversity.max.q,
                                      .min.q = immunarch.repdiversity.min.q,
                                      .q = immunarch.repdiversity.q,
                                      .step = immunarch.repdiversity.step,
                                      .quantile = immunarch.repdiversity.quantile,
                                      .extrapolation = immunarch.repdiversity.extrapolation,
                                      .perc = immunarch.repdiversity.perc,
                                      .norm = immunarch.repdiversity.norm,
                                      .do.norm = immunarch.repdiversity.do.norm,
                                      .laplace = immunarch.repdiversity.laplace))
    names(p14plots)[i] <- diversity_methods[i]
  }

  groups <- unique(contig_monarch$meta$CELL_TYPE)

  n <- 10
  p15plots <- NULL
  for(i in 1:length(contig_monarch$data)){
    tc1 <- trackClonotypes(contig_monarch$data,
                           list(i, n),
                           .col = immunarch.trackclonotypes.method1.col, # "nt"
                           .norm = immunarch.trackclonotypes.method1.norm)
    tc2 <- trackClonotypes(contig_monarch$data,
                           list(i, n),
                           .col = immunarch.trackclonotypes.method2.col, # "aa"
                           .norm = immunarch.trackclonotypes.method2.norm)
    p1 <- vis(tc1)+ggtitle(paste(names(contig_monarch$data)[i], ": TOP ",n," CLONOTYPES (",toupper(immunarch.trackclonotypes.method1.col),")", sep = ""))+theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=1))
    p2 <- vis(tc2)+ggtitle(paste(names(contig_monarch$data)[i], ": TOP ",n," CLONOTYPES (",toupper(immunarch.trackclonotypes.method2.col),")", sep = ""))+theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=1))
    p <- p1/p2
    print(p)
    p15plots[[i]] <- p
    names(p15plots)[i] <- names(contig_monarch$data)[i]
  }

  result_vdjdb <- dbAnnotate(contig_monarch$data, vdjdb, "CDR3.aa", "cdr3")

  # kmer

  p16data <- repOverlap(contig_monarch$data,
                        .col = "v+d+j",
                        .method = immunarch.repoverlap.method3.method,
                        .a = immunarch.repoverlap.method3.a,
                        .b = immunarch.repoverlap.method3.b,
                        .step = immunarch.repoverlap.method3.step,
                        .n.steps = immunarch.repoverlap.method3.n.steps,
                        .downsample = immunarch.repoverlap.method3.downsample,
                        .bootstrap = immunarch.repoverlap.method3.bootstrap,
                        .force.matrix = immunarch.repoverlap.method3.force.matrix)

  ov <- scale(p16data)
  ov <- t(scale(t(ov)))
  ov[is.na(ov)] <- 0
  p17plots <- complex_heatmap(ov, col = color_conditions$BlueYellowRed, legendtitle = "Z-Score", col_title = "Repetoire Overlaps (VDJ Regions)")

  print("Processing..")

  contig_types <- c("BCR","TCR")
  contig_table <- NULL

  bcr_list <- NULL

  if(length(grep("B_Cell|BCell|B.*Cell|BCR", pheno_data$CELL_TYPE, ignore.case = T)) > 0){
    bcr_list <- pheno_data[grepl("B_Cell|BCell|B.*Cell|BCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),"SID"]
  }

  # Very slow, take note of this
  bcr_pheno <- pheno_data[grepl("B_Cell|BCell|B.*Cell|BCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),]
  if(length(bcr_list) > 0){
    contig_bcr <- combineBCR(contig_list[which(names(contig_list) %in% bcr_list)],
                             samples = bcr_pheno[match(bcr_list,bcr_pheno$SID),"SID"],
                             ID = bcr_pheno[match(bcr_list,bcr_pheno$SID),"PAIR_ID"],
                             call.related.clones = screpertoire.combinebcr.call.related.clones,
                             threshold = screpertoire.combinebcr.threshold,
                             removeNA = screpertoire.combinebcr.removeNA,
                             removeMulti = screpertoire.combinebcr.removeMulti,
                             filterMulti = screpertoire.combinebcr.filterMulti)
    for(i in 1:length(contig_bcr)){
      contig_bcr[[i]]$ID <- bcr_pheno[match(bcr_list,bcr_pheno$SID),"INDIVIDUAL_ID"][i]
    }
    names(contig_bcr) <- bcr_pheno[match(bcr_list,bcr_pheno$SID),"SID"]
    # ID = bcr_pheno[match(bcr_list,bcr_pheno$SAMPLE_ID),"INDIVIDUAL_ID"])
    contig_table$BCR <- contig_bcr
  }else{
    contig_table$BCR <- NULL
  }

  print("Processing..")

  tcr_list <- NULL
  if(length(grep("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T)) > 0){
    tcr_list <- pheno_data[grepl("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),"SID"]
  }

  contig_tcr <- NULL
  tcr_pheno <- pheno_data[grepl("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),]
  if(length(tcr_list) > 0){
    contig_tcr <- contig_list[which(names(contig_list) %in% tcr_list)]
      contig_tcr <- combineTCR(contig_tcr,
                                   samples = tcr_pheno[match(names(contig_tcr),tcr_pheno$SID),"SID"],
                                   ID = tcr_pheno[match(names(contig_tcr),tcr_pheno$SID),"PAIR_ID"],
                                   removeNA = screpertoire.combinetcr.removeNA,
                                   removeMulti = screpertoire.combinetcr.removeMulti,
                                   filterMulti = screpertoire.combinetcr.filterMulti)
      for(i in 1:length(contig_tcr)){
        contig_tcr[[i]]$SID <- tcr_pheno[match(names(contig_tcr),tcr_pheno$SID),"SID"][i]
      }
      contig_table$TCR <- contig_tcr
  }

  # contig_bcr
  # contig_tcr
  p18plots <- NULL
  p19data <- NULL
  p20plots <- NULL
  p21plots <- NULL
  p22plots <- NULL
  p23plots <- NULL
  p24plots <- NULL
  p25data <- NULL
  p26plots <- NULL

  k <- 1
  name_ref <- data.frame(DATA_TYPE = c("BCR","TCR"),
                         DESC = c("BCR","TCR"))

  for(i in 1:length(contig_table)){

    current_type <- names(contig_table)[i]
    current_name <- name_ref[match(current_type, name_ref$DATA_TYPE),"DESC"]

    if(length(contig_table[[i]]) > 0){
      contig_table[[i]] <- lapply(contig_table[[i]], function(x){
        x <- data.frame(x, contig_pheno[match(unique(x$sample), contig_pheno$SID),])
      })

      number_contig <- quantContig(contig_table[[i]],
                                   exportTable = T,
                                   cloneCall = screpertoire.quantcontig.cloneCall,
                                   chain = screpertoire.quantcontig.chain,
                                   scale = screpertoire.quantcontig.scale,
                                   group.by = screpertoire.quantcontig.group.by,
                                   split.by = screpertoire.quantcontig.split.by,
                                   order = screpertoire.quantcontig.order)

      p1 <- abundanceContig(contig_table[[i]],
                            cloneCall = screpertoire.abundancecontig.cloneCall,
                            chain = screpertoire.abundancecontig.chain,
                            scale = FALSE,
                            group.by = screpertoire.abundancecontig.group.by,
                            split.by = screpertoire.abundancecontig.split.by,
                            order = screpertoire.abundancecontig.order,
                            exportTable = FALSE) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              plot.title = element_text(size = 20, face = "bold")) +
        guides(fill = guide_legend(ncol = 1))

      p2 <- abundanceContig(contig_table[[i]],
                            cloneCall = screpertoire.abundancecontig.cloneCall,
                            chain = screpertoire.abundancecontig.chain,
                            scale = TRUE,
                            group.by = screpertoire.abundancecontig.group.by,
                            split.by = screpertoire.abundancecontig.split.by,
                            order = screpertoire.abundancecontig.order,
                            exportTable = FALSE) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              plot.title = element_text(size = 20, face = "bold")) +
        guides(fill = guide_legend(ncol = 1))

      p18plots[[length(p18plots)+1]] <- (p1|p2)+plot_annotation(title = paste(current_name, " CLONOTYPE ABUNDANCE (",toupper(screpertoire.abundancecontig.cloneCall),",Chain = ",screpertoire.abundancecontig.chain,")", sep = ""),theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      names(p18plots)[length(p18plots)] <- current_name

      current <- abundanceContig(contig_table[[i]], cloneCall = screpertoire.abundancecontig.cloneCall,
                                 chain = screpertoire.abundancecontig.chain,
                                 scale = FALSE,
                                 group.by = screpertoire.abundancecontig.group.by,
                                 split.by = screpertoire.abundancecontig.split.by,
                                 order = screpertoire.abundancecontig.order,
                                 exportTable = TRUE)
      p19data[[length(p19data)+1]] <- current[order(current$Abundance, decreasing = T),]
      names(p19data)[length(p19data)] <- current_name

      p1 <- lengthContig(contig_table[[i]],
                         cloneCall = screpertoire.lengthcontig.choice1.cloneCall,
                         chain = screpertoire.lengthcontig.choice1.chain,
                         group.by = screpertoire.lengthcontig.choice1.group.by,
                         split.by = screpertoire.lengthcontig.choice1.split.by,
                         order = screpertoire.lengthcontig.choice1.order,
                         scale = screpertoire.lengthcontig.choice1.scale,
                         exportTable = FALSE)+
        ggtitle("CDR3 AA")+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              axis.text.x = element_text(angle = 0, hjust=1),
              plot.title = element_text(size = 20, face = "bold")) +
        guides(fill = guide_legend(ncol = 1))+
        scale_fill_tableau()
      p1 <- adjust_theme(p1)

      p2 <- lengthContig(contig_table[[i]],
                         cloneCall = screpertoire.lengthcontig.choice2.cloneCall,
                         chain = screpertoire.lengthcontig.choice2.chain,
                         group.by = screpertoire.lengthcontig.choice2.group.by,
                         split.by = screpertoire.lengthcontig.choice2.split.by,
                         order = screpertoire.lengthcontig.choice2.order,
                         scale = screpertoire.lengthcontig.choice2.scale,
                         exportTable = FALSE)+
        ggtitle("CDR3 NT")+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              axis.text.x = element_text(angle = 0, hjust=1),
              plot.title = element_text(size = 20, face = "bold"))+
        guides(fill = guide_legend(ncol = 1))+
        scale_fill_tableau()
      p2 <- adjust_theme(p2)

      p20plots[[length(p20plots)+1]] <- (p1/p2)+plot_annotation(title = paste(current_name," CONTIG LENGTH", sep = ""),
                                                                theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      names(p20plots)[length(p20plots)] <- current_name

      n <- 10
      temp <- compareClonotypes(contig_table[[i]],
                                cloneCall = screpertoire.compareclonotypes.choice1.cloneCall,
                                chain = screpertoire.compareclonotypes.choice1.chain,
                                samples = screpertoire.compareclonotypes.choice1.samples,
                                numbers = screpertoire.compareclonotypes.choice1.numbers,
                                split.by = screpertoire.compareclonotypes.choice1.split.by,
                                graph = screpertoire.compareclonotypes.choice1.graph,
                                exportTable = TRUE)
      # temp <- compareClonotypes_debugged(contig_table[[i]], cloneCall = "aa")
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)
      p21plots[[k]] <- compareClonotypes(contig_table[[i]],
                                         clonotypes = unique(as.character(temp$Clonotypes)),
                                         cloneCall = screpertoire.compareclonotypes.choice1.cloneCall,
                                         chain = screpertoire.compareclonotypes.choice1.chain,
                                         samples = screpertoire.compareclonotypes.choice1.samples,
                                         numbers = screpertoire.compareclonotypes.choice1.numbers,
                                         split.by = screpertoire.compareclonotypes.choice1.split.by,
                                         graph = screpertoire.compareclonotypes.choice1.graph,
                                         exportTable = FALSE)
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]]+
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = ceiling(length(unique(temp$Clonotypes))/12))) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (",toupper(screpertoire.compareclonotypes.choice1.cloneCall),")", sep = ""))
      names(p21plots)[k] <- paste(toupper(screpertoire.compareclonotypes.choice1.cloneCall),"_",current_name, sep = "")
      k <- k+1

      temp <- compareClonotypes(contig_table[[i]],
                                cloneCall = screpertoire.compareclonotypes.choice2.cloneCall,
                                chain = screpertoire.compareclonotypes.choice2.chain,
                                samples = screpertoire.compareclonotypes.choice2.samples,
                                numbers = screpertoire.compareclonotypes.choice2.numbers,
                                split.by = screpertoire.compareclonotypes.choice2.split.by,
                                graph = screpertoire.compareclonotypes.choice2.graph,
                                exportTable = TRUE)
      temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
      temp <- temp[order(temp$Proportion, decreasing = T),]
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)

      p21plots[[k]] <- compareClonotypes(contig_table[[i]],
                                         clonotypes = unique(as.character(temp$Clonotypes)),
                                         cloneCall = screpertoire.compareclonotypes.choice2.cloneCall,
                                         chain = screpertoire.compareclonotypes.choice2.chain,
                                         samples = screpertoire.compareclonotypes.choice2.samples,
                                         numbers = screpertoire.compareclonotypes.choice2.numbers,
                                         split.by = screpertoire.compareclonotypes.choice2.split.by,
                                         graph = screpertoire.compareclonotypes.choice2.graph,
                                         exportTable = FALSE)
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]] +
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = ceiling(length(unique(temp$Clonotypes))/12))) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (",toupper(screpertoire.compareclonotypes.choice2.cloneCall),")", sep = ""))
      names(p21plots)[k] <- paste(toupper(screpertoire.compareclonotypes.choice2.cloneCall),"_",current_name, sep = "")
      k <- k+1

      temp <- compareClonotypes(contig_table[[i]],
                                cloneCall = screpertoire.compareclonotypes.choice3.cloneCall,
                                chain = screpertoire.compareclonotypes.choice3.chain,
                                samples = screpertoire.compareclonotypes.choice3.samples,
                                numbers = screpertoire.compareclonotypes.choice3.numbers,
                                split.by = screpertoire.compareclonotypes.choice3.split.by,
                                graph = screpertoire.compareclonotypes.choice3.graph,
                                exportTable = TRUE)
      temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
      temp <- temp[order(temp$Proportion, decreasing = T),]
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)

      p21plots[[k]] <- compareClonotypes(contig_table[[i]],
                                         clonotypes = unique(as.character(temp$Clonotypes)),
                                         cloneCall = screpertoire.compareclonotypes.choice3.cloneCall,
                                         chain = screpertoire.compareclonotypes.choice3.chain,
                                         samples = screpertoire.compareclonotypes.choice3.samples,
                                         numbers = screpertoire.compareclonotypes.choice3.numbers,
                                         split.by = screpertoire.compareclonotypes.choice3.split.by,
                                         graph = screpertoire.compareclonotypes.choice3.graph,
                                         exportTable = FALSE)
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]] +
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (",toupper(screpertoire.compareclonotypes.choice3.cloneCall),")", sep = ""))
      names(p21plots)[k] <- paste(toupper(screpertoire.compareclonotypes.choice3.cloneCall),"_",current_name, sep = "")
      k <- k+1

      temp <- compareClonotypes(contig_table[[i]],
                                cloneCall = screpertoire.compareclonotypes.choice4.cloneCall,
                                chain = screpertoire.compareclonotypes.choice4.chain,
                                samples = screpertoire.compareclonotypes.choice4.samples,
                                numbers = screpertoire.compareclonotypes.choice4.numbers,
                                split.by = screpertoire.compareclonotypes.choice4.split.by,
                                graph = screpertoire.compareclonotypes.choice4.graph,
                                exportTable = TRUE)
      temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
      temp <- temp[order(temp$Proportion, decreasing = T),]
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)

      p21plots[[k]] <- compareClonotypes(contig_table[[i]],
                                         clonotypes = unique(as.character(temp$Clonotypes)),
                                         cloneCall = screpertoire.compareclonotypes.choice4.cloneCall,
                                         chain = screpertoire.compareclonotypes.choice4.chain,
                                         samples = screpertoire.compareclonotypes.choice4.samples,
                                         numbers = screpertoire.compareclonotypes.choice4.numbers,
                                         split.by = screpertoire.compareclonotypes.choice4.split.by,
                                         graph = screpertoire.compareclonotypes.choice4.graph,
                                         exportTable = FALSE)
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]] +
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (",toupper(screpertoire.compareclonotypes.choice4.cloneCall),")", sep = ""))
      names(p21plots)[k] <- paste(toupper(screpertoire.compareclonotypes.choice4.cloneCall),"_",current_name, sep = "")
      k <- k+1

      p <- NULL
      p <- vizGenes(contig_table[[i]],
                    gene = screpertoire.vizgenes.gene,
                    chain = screpertoire.vizgenes.chain,
                    plot = screpertoire.vizgenes.plot,
                    y.axis = screpertoire.vizgenes.y.axis,
                    order = screpertoire.vizgenes.order,
                    scale = screpertoire.vizgenes.scale,
                    group.by = screpertoire.vizgenes.group.by,
                    split.by = screpertoire.vizgenes.split.by,
                    exportTable = FALSE)+
        ggtitle(paste(current_name,": ",toupper(screpertoire.vizgenes.gene)," GENE USAGE", sep = ""))+
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              plot.title = element_text(size = 15))
      p22plots[[i]] <- p
      names(p22plots)[i] <- current_name

      p1 <- NULL
      p2 <- NULL
      p1 <- clonalHomeostasis(contig_table[[i]],
                              cloneTypes = screpertoire.clonalhomeostasis.cloneTypes,
                              cloneCall = screpertoire.clonalhomeostasis.cloneCall,
                              chain = screpertoire.clonalhomeostasis.chain,
                              group.by = screpertoire.clonalhomeostasis.group.by,
                              split.by = screpertoire.clonalhomeostasis.split.by,
                              exportTable = FALSE) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 20, face = "bold"))+
        scale_x_discrete(labels= names(contig_table[[i]]))
      p2 <- clonalProportion(contig_table[[i]],
                             split = screpertoire.clonalproportion.split,
                             cloneCall = screpertoire.clonalproportion.cloneCall,
                             chain = screpertoire.clonalproportion.chain,
                             group.by = screpertoire.clonalproportion.group.by,
                             split.by = screpertoire.clonalproportion.split.by,
                             exportTable = FALSE) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45, hjust=1))
      p23plots[[i]] <- p1 + p2 + plot_annotation(title = paste(current_name, sep = ""), theme =  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      names(p23plots)[i] <- current_name

      p24plots[[i]] <- clonalOverlap(contig_table[[i]],
                                     cloneCall = screpertoire.clonaloverlap.cloneCall,
                                     method = screpertoire.clonaloverlap.method,
                                     chain = screpertoire.clonaloverlap.chain,
                                     split.by = screpertoire.clonaloverlap.split.by,
                                     exportTable = FALSE) +
        scale_fill_continuous_tableau(palette = "Blue-Teal") +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
              plot.title = element_text(size = 15, face = "bold")) +
        ggtitle(paste(current_name, ": ",toupper(screpertoire.clonaloverlap.method)," FOR CLONAL OVERLAP (",toupper(screpertoire.clonaloverlap.cloneCall),")", sep = ""))
      p24plots[[i]] <- adjust_theme(p24plots[[i]], title_size = 15)
      names(p24plots)[i] <- current_name

      t4 <- clonesizeDistribution(contig_table[[i]],
                                  cloneCall = screpertoire.clonesizedistribution.cloneCall,
                                  chain = screpertoire.clonesizedistribution.chain,
                                  method = screpertoire.clonesizedistribution.method,
                                  threshold = screpertoire.clonesizedistribution.threshold,
                                  group.by = screpertoire.clonesizedistribution.group.by,
                                  split.by = screpertoire.clonesizedistribution.split.by,
                                  exportTable = TRUE)
      hclust <- hclust(as.dist(t4), method = screpertoire.clonesizedistribution.hclust.method)
      p25data[[i]] <- as.dendrogram(hclust)
      names(p25data)[i] <- current_name # JENSEN-SHANNON DISTANCE CLUSTERING (VDJC GENE + CDR3 NT)

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]],
                                                        cloneCall = "aa",
                                                        group.by = "sample",
                                                        chain = screpertoire.clonesizediversity.chain,
                                                        x.axis = screpertoire.clonesizediversity.x.axis,
                                                        split.by = screpertoire.clonesizediversity.split.by,
                                                        exportTable = screpertoire.clonesizediversity.exportTable,
                                                        n.boots = screpertoire.clonesizediversity.n.boots,
                                                        return.boots = screpertoire.clonesizediversity.return.boots,
                                                        skip.boots = screpertoire.clonesizediversity.skip.boots)+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (CDR3 AA)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("AA_", current_name, sep = "")

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]],
                                                        cloneCall = "gene",
                                                        group.by = "sample",
                                                        chain = screpertoire.clonesizediversity.chain,
                                                        x.axis = screpertoire.clonesizediversity.x.axis,
                                                        split.by = screpertoire.clonesizediversity.split.by,
                                                        exportTable = screpertoire.clonesizediversity.exportTable,
                                                        n.boots = screpertoire.clonesizediversity.n.boots,
                                                        return.boots = screpertoire.clonesizediversity.return.boots,
                                                        skip.boots = screpertoire.clonesizediversity.skip.boots)+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (VDJC GENE)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("GENE_", current_name, sep = "")

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]],
                                                        cloneCall = "nt",
                                                        group.by = "sample",
                                                        chain = screpertoire.clonesizediversity.chain,
                                                        x.axis = screpertoire.clonesizediversity.x.axis,
                                                        split.by = screpertoire.clonesizediversity.split.by,
                                                        exportTable = screpertoire.clonesizediversity.exportTable,
                                                        n.boots = screpertoire.clonesizediversity.n.boots,
                                                        return.boots = screpertoire.clonesizediversity.return.boots,
                                                        skip.boots = screpertoire.clonesizediversity.skip.boots)+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (CDR3 NT)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("NT_", current_name, sep = "")

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]],
                                                        cloneCall = "strict",
                                                        group.by = "sample",
                                                        chain = screpertoire.clonesizediversity.chain,
                                                        x.axis = screpertoire.clonesizediversity.x.axis,
                                                        split.by = screpertoire.clonesizediversity.split.by,
                                                        exportTable = screpertoire.clonesizediversity.exportTable,
                                                        n.boots = screpertoire.clonesizediversity.n.boots,
                                                        return.boots = screpertoire.clonesizediversity.return.boots,
                                                        skip.boots = screpertoire.clonesizediversity.skip.boots)+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (VDJC GENE + CDR3 NT)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("GENE_NT_", current_name, sep = "")

    }
  }

  print("Running scRNA-Seq..")

  if(length(rna_list) > 0){
    data_current <- NULL
    for(j in 1:nrow(rna_pheno)){
      current_name <- rna_pheno[j,"SID"]
      current_file <- NULL
      current_file <- rna_list[[grep(rna_pheno[j,"FILE"], rna_list, ignore.case = T)]]
      if(length(grep("\\.RDS$", current_file, ignore.case = T)) > 0){
        data_current[[j]] <- readRDS(current_file)
      }else{
        data_current[[j]] <- Read10X_h5(current_file)
      }

      if(length(data_current[[j]]) <= 5){
        temp <- data_current[[j]]$`Gene Expression`
      }else{
        temp <- data_current[[j]]
      }
      if(length(grep("Seurat", class(data_current[[j]]), ignore.case = T)) == 0){
        data_current[[j]] <- CreateSeuratObject(counts = temp, project = project_name,
                                                min.cells = seurat.createseuratobject.min.cells)
      }
      names(data_current)[j] <- current_name
      DefaultAssay(data_current[[j]]) <- "RNA"
      Idents(data_current[[j]]) <- "orig.ident"
      data_current[[j]][["percent.mt"]] <- PercentageFeatureSet(data_current[[j]], pattern = "^MT-")
      data_current[[j]] <- subset(data_current[[j]],
                                  subset = nFeature_RNA > seurat.subset.min.nFeature_RNA &
                                  nFeature_RNA < seurat.subset.max.nFeature_RNA &
                                  percent.mt < seurat.subset.max.percent.mt)
      data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, rna_pheno[match(current_name, rna_pheno$SID),])
      data_current[[j]] <- NormalizeData(data_current[[j]],
                                         assay = seurat.normalizedata.assay,
                                         normalization.method = seurat.normalizedata.normalization.method,
                                         scale.factor = seurat.normalizedata.scale.factor,
                                         margin = seurat.normalizedata.margin)
      data_current[[j]] <- FindVariableFeatures(data_current[[j]],
                                                selection.method = seurat.findvariablefeatures.selection.method,
                                                loess.span = seurat.findvariablefeatures.loess.span,
                                                clip.max = seurat.findvariablefeatures.clip.max,
                                                # mean.function = seurat.findvariablefeatures.mean.function,
                                                # dispersion.function = seurat.findvariablefeatures.dispersion.function,
                                                num.bin = seurat.findvariablefeatures.num.bin,
                                                binning.method = seurat.findvariablefeatures.binning.method,
                                                nfeatures = seurat.findvariablefeatures.nfeatures,
                                                mean.cutoff = seurat.findvariablefeatures.mean.cutoff,
                                                dispersion.cutoff = seurat.findvariablefeatures.dispersion.cutoff)
      print(paste("Complete: ", rna_pheno[j,"SID"], "..", sep = ""))
    }

    if(length(data_current) > 1){
      data_current <- lapply(data_current, function(x) {
        x <- ScaleData(x, verbose=FALSE,
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
                       min.cells.to.block = seurat.scaledata.min.cells.to.block)
        x <- RunPCA(x,
                    verbose = FALSE,
                    assay = seurat.runpca.assay,
                    features = seurat.runpca.features,
                    npcs = seurat.runpca.npcs,
                    rev.pca = seurat.runpca.rev.pca,
                    weight.by.var = seurat.runpca.weight.by.var,
                    ndims.print = seurat.runpca.ndims.print,
                    nfeatures.print = seurat.runpca.nfeatures.print,
                    reduction.name = seurat.runpca.reduction.name,
                    reduction.key = seurat.runpca.reduction.key,
                    seed.use = seurat.runpca.seed.use)
      })

      if(toupper(integration_method) == "SEURAT" | is.null(integration_method)){
        red_method <- "pca"
        integrative_features <- SelectIntegrationFeatures(object.list = data_current,
                                                          nfeatures = seurat.selectintegrationfeatures.nfeatures,
                                                          assay = seurat.selectintegrationfeatures.assay,
                                                          fvf.nfeatures = seurat.selectintegrationfeatures.fvf.nfeatures)
        data_anchors <- FindIntegrationAnchors(object.list = data_current,
                                               anchor.features = integrative_features,
                                               assay = seurat.findintegrationanchors.assay,
                                               reference = seurat.findintegrationanchors.reference,
                                               scale = seurat.findintegrationanchors.scale,
                                               normalization.method = seurat.findintegrationanchors.normalization.method,
                                               sct.clip.range = seurat.findintegrationanchors.sct.clip.range,
                                               reduction = seurat.findintegrationanchors.reduction,
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
                              new.assay.name = seurat.integratedata.new.assay.name,
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
        DefaultAssay(data) <- "integrated"
        rm(data_anchors)
      }else if(toupper(integration_method) == "HARMONY"){
        data <- merge(data_current[[1]], data_current[c(2:length(data_current))]) # , add.cell.ids = names(data_current)
        data <- NormalizeData(data,
                              verbose = FALSE,
                              assay = seurat.normalizedata.assay,
                              normalization.method = seurat.normalizedata.normalization.method,
                              scale.factor = seurat.normalizedata.scale.factor,
                              margin = seurat.normalizedata.margin) %>%
          FindVariableFeatures(selection.method = seurat.findvariablefeatures.selection.method,
                               loess.span = seurat.findvariablefeatures.loess.span,
                               clip.max = seurat.findvariablefeatures.clip.max,
                               # mean.function = seurat.findvariablefeatures.mean.function,
                               # dispersion.function = seurat.findvariablefeatures.dispersion.function,
                               num.bin = seurat.findvariablefeatures.num.bin,
                               binning.method = seurat.findvariablefeatures.binning.method,
                               nfeatures = seurat.findvariablefeatures.nfeatures,
                               mean.cutoff = seurat.findvariablefeatures.mean.cutoff,
                               dispersion.cutoff = seurat.findvariablefeatures.dispersion.cutoff) %>%
          ScaleData(verbose = FALSE,
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
                    min.cells.to.block = seurat.scaledata.min.cells.to.block) %>%
          RunPCA(verbose = FALSE,
                 assay = seurat.runpca.assay,
                 features = seurat.runpca.features,
                 npcs = seurat.runpca.npcs,
                 rev.pca = seurat.runpca.rev.pca,
                 weight.by.var = seurat.runpca.weight.by.var,
                 ndims.print = seurat.runpca.ndims.print,
                 nfeatures.print = seurat.runpca.nfeatures.print,
                 reduction.name = seurat.runpca.reduction.name,
                 reduction.key = seurat.runpca.reduction.key,
                 seed.use = seurat.runpca.seed.use)
        
        data <- RunHarmony(data,
                           reduction = "pca",
                           dims.use = harmony.runharmony.dims,
                           group.by.vars = harmony.runharmony.group.by.vars,
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
                           project.dim = harmony.runharmony.project.dim)
        red_method <- "harmony"
        DefaultAssay(data) <- "RNA"
      }
      rm(data_current)

    }else if(length(data_current) == 1){
      data <- data_current[[1]]
      rm(data_current)
      red_method <- "pca"
      DefaultAssay(data) <- "RNA"
    }

    print("data:")
    print(data)
    data <- FindVariableFeatures(data,
                                 nfeatures = seurat.integration.findvariablefeatures.nfeatures,
                                 selection.method = seurat.integration.findvariablefeatures.selection.method,
                                 loess.span = seurat.integration.findvariablefeatures.loess.span,
                                 clip.max = seurat.integration.findvariablefeatures.clip.max,
                                 # mean.function = seurat.integration.findvariablefeatures.mean.function,
                                 # dispersion.function = seurat.integration.findvariablefeatures.dispersion.function,
                                 num.bin = seurat.integration.findvariablefeatures.num.bin,
                                 binning.method = seurat.integration.findvariablefeatures.binning.method)
    data <- ScaleData(data,
                      verbose = FALSE,
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
                      min.cells.to.block = seurat.integration.scaledata.min.cells.to.block)
    data <- RunPCA(data,
                   verbose = FALSE,
                   features = seurat.integration.runpca.features,
                   npcs = seurat.integration.runpca.npcs,
                   rev.pca = seurat.integration.runpca.rev.pca,
                   weight.by.var = seurat.integration.runpca.weight.by.var,
                   ndims.print = seurat.integration.runpca.ndims.print,
                   nfeatures.print = seurat.integration.runpca.nfeatures.print,
                   reduction.key = seurat.integration.runpca.reduction.key,
                   seed.use = seurat.integration.runpca.seed.use,
                   approx = seurat.integration.runpca.approx)
    
    selected_pcs <- NULL
    selected_pcs <- find_selected_pcs(data,
                                      pcs = seurat.integration.runumap.pcs,
                                      pca = "harmony",
                                      cname = "integrated data")
    
    data <- RunUMAP(data,
                    reduction = red_method,
                    dims = 1:selected_pcs,
                    features = seurat.integration.runumap.features,
                    graph = seurat.integration.runumap.graph,
                    # assay = seurat.integration.runumap.assay,
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
    
    data <- FindNeighbors(data,
                          reduction = red_method,
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
                         edge.file.name = seurat.integration.findclusters.edge.file.name)
    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(data$seurat_clusters)))
    names(cluster_colors) <- c(sort(as.numeric(as.character(unique(data$seurat_clusters)))))

    print("Running SingleR..")
    Idents(data) <- "seurat_clusters"
    DefaultAssay(data) <- "RNA"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
                       clusters =  data$seurat_clusters,
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
    data$CELL_TYPE <- clu_ann$labels[match(data$seurat_clusters,row.names(clu_ann))]
    data@meta.data[which(is.na(data$CELL_TYPE)),"CELL_TYPE"] <- "Unidentifiable"
    Idents(data) <- data$CELL_TYPE

    celltype_colors <- gen_colors(color_conditions$monet, length(unique(data$CELL_TYPE)))
    names(celltype_colors) <- sort(as.character(unique(data$CELL_TYPE)))

    plotx <- gen10x_plotx(data, selected = c("PCA","UMAP"), include_meta = T)
    plotx$CLUSTER <- plotx$seurat_clusters

    group_colors <- gen_colors(color_conditions$tableau20, length(unique(data$GROUP)))
    names(group_colors) <- sort(as.character(unique(data$GROUP)))

    p27plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "SID", plot_title = project_name,
                             point_size = 1, numeric = F,label_size = 10,
                             col = sample_colors, annot = F, legend_position = "right")
    p27plots <- adjust_theme(p27plots)

    p29plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER", plot_title = project_name,
                             point_size = 1, numeric = T,label_size = 10,
                             col = cluster_colors, annot = TRUE, legend_position = "right")
    p29plots <- adjust_theme(p29plots)

    p30plots <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", group_by = "GROUP",
                                  color_by = "CLUSTER", isfacet = T, xlabel = "UMAP_1", ylabel = "UMAP_2",point_size = 1,
                                  title = project_name,col = cluster_colors)
    p30plots <- adjust_theme(p30plots)

    p31plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = project_name,
                             point_size = 1,label_size = 10,
                             col = celltype_colors, annot = TRUE, legend_position = "right")
    p31plots <- adjust_theme(p31plots)

    p32plots <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", group_by = "GROUP",
                                  color_by = "CELL_TYPE", isfacet = T, xlabel = "UMAP_1", ylabel = "UMAP_2",
                                  point_size = 1,
                                  title = project_name,col = celltype_colors)
    p32plots <- adjust_theme(p32plots)

    orig_names <- row.names(data@meta.data)
    row.names(data@meta.data) <- paste(gsub("_CONTIG$|_SCRNA$","",data$SID, ignore.case = T),
                                       data$PAIR_ID,
                                       row.names(data@meta.data), sep = "_")
    row.names(data@meta.data) <- gsub("_[0-9]+$","",row.names(data@meta.data), ignore.case = T)
    cmultiplexs <- rna_pheno[grep("ALL.*CELL", rna_pheno$CELL_TYPE, ignore.case = T),"SID"]
    cmultiplex_tbcr <- pheno_data[which(pheno_data$PAIR_ID %in% unique(pheno_data[which(pheno_data$SID %in% cmultiplexs),"PAIR_ID"]) &
                                          grepl("CONTIG",pheno_data$DATA_TYPE, ignore.case = T)),"SID"]
    if(length(cmultiplexs) > 0){
      for(i in 1:length(cmultiplexs)){
        row.names(data@meta.data) <- gsub(cmultiplexs[i],gsub("_ALL.*CELL.*","",cmultiplexs[i]),row.names(data@meta.data), ignore.case = T)
      }
    }

    current_contigs <- NULL
    for(i in 1:length(contig_table)){
      current <- contig_table[[i]]
      for(j in 1:length(current)){
        print(paste(i,",",j,sep = ""))
        if(length(which(toupper(names(current)[j]) %in% toupper(cmultiplex_tbcr))) > 0){
          current[[j]]$barcode <- gsub("_BCR|_TCR","",current[[j]]$barcode, ignore.case = T)
        }
        current[[j]]$barcode <- gsub("_BCR|_TCR","",current[[j]]$barcode, ignore.case = T)
        current[[j]]$barcode <- gsub("_CONTIG|_CONSENSUS","",current[[j]]$barcode, ignore.case = T)
        current[[j]]$barcode <- gsub("_[0-9]+$","",current[[j]]$barcode, ignore.case = T)
        data@meta.data[which(row.names(data@meta.data) %in% current[[j]]$barcode),"DATA_TYPE"] <- unique(current[[j]][, grep("CELL_TYPE", colnames(current[[j]]), ignore.case = T)])[1]
      }
      current_contigs <- c(current_contigs, current)
    }
    data <- combineExpression(current_contigs,
                              data,
                              cloneCall = screpertoire.combineexpression.cloneCall,
                              chain = screpertoire.combineexpression.chain,
                              group.by = screpertoire.combineexpression.group.by,
                              proportion = screpertoire.combineexpression.proportion,
                              filterNA = screpertoire.combineexpression.filterNA,
                              cloneTypes = screpertoire.combineexpression.cloneTypes,
                              addLabel = screpertoire.combineexpression.addLabel)
    data@meta.data[which(!toupper(data$DATA_TYPE) %in% toupper(c("BCR","TCR","B","T","T-AB","TGD"))),"DATA_TYPE"] <- "OTHERS"

    if(screpertoire.combineexpression.proportion == TRUE){
      current_groups <- c("Hyperexpanded (0.1 < X <= 1)",
                          "Large (0.01 < X <= 0.1)",
                          "Medium (0.001 < X <= 0.01)",
                          "Small (1e-04 < X <= 0.001)",
                          "Rare (0 < X <= 1e-04)", NA)
    }else{
      current_groups <- c("Hyperexpanded (X > 100)",
                          "Large (20 < X <= 100)",
                          "Medium (5 < X <= 20)",
                          "Small (1 < X <= 5)",
                          "Single (0 < X <= 1)", NA)
    }

    Idents(data) <- data$CELL_TYPE
    plotx <- gen10x_plotx(data, selected = c("PCA","UMAP"), include_meta = T)

    clone_cols <- c(rev(c("#1D47A0","#8BC3FA","#D1FBEC","#F3B750","#EB5A35")),"grey")
    names(clone_cols) <- current_groups

    p33plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "cloneType", plot_title = project_name,
                             point_size = 1, numeric = F,label_size = 10,
                             col = clone_cols, annot = "cloneType", legend_position = "right")
    p33plots <- adjust_theme(p33plots)

    datatype_colors <- gen_colors(color_conditions$bright, length(unique(data$DATA_TYPE)))
    names(datatype_colors) <- sort(as.character(unique(data$DATA_TYPE)))

    p28plots <- own_facet_scatter(plotx, feature1 = "UMAP_1", feature2 = "UMAP_2", group_by = "DATA_TYPE",
                                  color_by = "DATA_TYPE", isfacet = T, xlabel = "UMAP_1", ylabel = "UMAP_2",point_size = 1,
                                  title = project_name,col = datatype_colors)
    p28plots <- adjust_theme(p28plots)

    # current <- occupiedscRepertoire(data, x.axis = "CELL_TYPE", exportTable = T, proportion = F)
    current <- data.frame(table(data@meta.data[,c("cloneType","CELL_TYPE", "GROUP", "DATA_TYPE")]))
    current <- current[current$Freq > 0,]
    current <- split(current, current$CELL_TYPE)
    for(k in 1:length(current)){
      temp <- current[[k]]
      temp <- split(temp, temp$GROUP)
      for(m in 1:length(temp)){
        temp2 <- temp[[m]]
        temp2 <- split(temp2, temp2$DATA_TYPE)
        temp2 <- lapply(temp2, function(x){
          x <- data.frame(x, Prop = x$Freq/sum(x$Freq))
        })
        temp2 <- do.call(rbind.data.frame, temp2)
        temp[[m]] <- temp2
      }
      temp <- do.call(rbind.data.frame, temp)
      current[[k]] <- temp
    }
    current <- do.call(rbind.data.frame, current)

    print("Running occupiedscRepertoire..")
    Idents(data) <- "CELL_TYPE"
    p1 <- occupiedscRepertoire(data,
                               x.axis = "CELL_TYPE",
                               proportion = FALSE,
                               label = screpertoire.occupiedscrepertoire.label,
                               facet.by = screpertoire.occupiedscrepertoire.facet.by,
                               na.include = screpertoire.occupiedscrepertoire.na.include,
                               exportTable = screpertoire.occupiedscrepertoire.exportTable)+
      scale_fill_manual(values = clone_cols)+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 15),
            legend.position = "none",
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1.2, "cm"),
            axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 20),
            # axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))

    p2 <- occupiedscRepertoire(data,
                               x.axis = "CELL_TYPE",
                               proportion = TRUE,
                               label = screpertoire.occupiedscrepertoire.label,
                               facet.by = screpertoire.occupiedscrepertoire.facet.by,
                               na.include = screpertoire.occupiedscrepertoire.na.include,
                               exportTable = screpertoire.occupiedscrepertoire.exportTable)+
      scale_fill_manual(values = clone_cols)+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 15),
            # legend.position = "right",
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1.2, "cm"),
            axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 20),
            # axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))

    p34plots <- p1+p2

    p1 <- ggplot(current, aes_string(x = "CELL_TYPE", y = "Freq", fill = "cloneType")) +
      geom_bar(stat = "identity") + ylab("Single Cells") + theme_classic() +
      theme(axis.title.x = element_blank()) + geom_text(aes(label = Freq), position = position_stack(vjust = 0.5))+
      facet_wrap(~GROUP+DATA_TYPE)+
      scale_fill_manual(values = clone_cols)+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 15),
            legend.position = "none",
            strip.background = element_blank(),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1.2, "cm"),
            axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))

    p2 <- ggplot(current, aes_string(x = "CELL_TYPE", y = "Prop", fill = "cloneType")) +
      geom_bar(stat = "identity") + ylab("Single Cells") + theme_classic() +
      theme(axis.title.x = element_blank()) + geom_text(aes(label = Freq), position = position_stack(vjust = 0.5))+
      facet_wrap(~GROUP+DATA_TYPE)+
      scale_fill_manual(values = clone_cols)+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 15),
            # legend.position = "right",
            strip.background = element_blank(),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.key.size = unit(1.2, "cm"),
            axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))

    p35plots <- p1+p2

    combined_meta <- expression2List(data, split.by = "CELL_TYPE")
    combined_meta <- combined_meta[unlist(lapply(combined_meta, nrow)) > 0]
    p36plots <- NULL
    p36plots$GENE_NT <- clonalDiversity(combined_meta,
                                        cloneCall = "strict",
                                        chain = screpertoire.clonaldiversity.chain,
                                        group.by = screpertoire.clonaldiversity.group.by,
                                        x.axis = screpertoire.clonaldiversity.x.axis,
                                        split.by = screpertoire.clonaldiversity.split.by,
                                        exportTable = screpertoire.clonaldiversity.exportTable,
                                        n.boots = screpertoire.clonaldiversity.n.boots,
                                        return.boots = screpertoire.clonaldiversity.return.boots,
                                        skip.boots = screpertoire.clonaldiversity.skip.boots) +
      ggtitle("CLONAL DIVERISTY: VDJC GENE + CDR3 NUCLEOTIDE") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    p36plots$GENE <- clonalDiversity(combined_meta,
                                     cloneCall = "gene",
                                     chain = screpertoire.clonaldiversity.chain,
                                     group.by = screpertoire.clonaldiversity.group.by,
                                     x.axis = screpertoire.clonaldiversity.x.axis,
                                     split.by = screpertoire.clonaldiversity.split.by,
                                     exportTable = screpertoire.clonaldiversity.exportTable,
                                     n.boots = screpertoire.clonaldiversity.n.boots,
                                     return.boots = screpertoire.clonaldiversity.return.boots,
                                     skip.boots = screpertoire.clonaldiversity.skip.boots) +
      ggtitle("CLONAL DIVERISTY: VDJC GENE") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    p36plots$NT <- clonalDiversity(combined_meta,
                                   cloneCall = "nt",
                                   chain = screpertoire.clonaldiversity.chain,
                                   group.by = screpertoire.clonaldiversity.group.by,
                                   x.axis = screpertoire.clonaldiversity.x.axis,
                                   split.by = screpertoire.clonaldiversity.split.by,
                                   exportTable = screpertoire.clonaldiversity.exportTable,
                                   n.boots = screpertoire.clonaldiversity.n.boots,
                                   return.boots = screpertoire.clonaldiversity.return.boots,
                                   skip.boots = screpertoire.clonaldiversity.skip.boots) +
      ggtitle("CLONAL DIVERISTY: CDR3 NUCLEOTIDE") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    p36plots$AA <- clonalDiversity(combined_meta,
                                   cloneCall = "aa",
                                   chain = screpertoire.clonaldiversity.chain,
                                   group.by = screpertoire.clonaldiversity.group.by,
                                   x.axis = screpertoire.clonaldiversity.x.axis,
                                   split.by = screpertoire.clonaldiversity.split.by,
                                   exportTable = screpertoire.clonaldiversity.exportTable,
                                   n.boots = screpertoire.clonaldiversity.n.boots,
                                   return.boots = screpertoire.clonaldiversity.return.boots,
                                   skip.boots = screpertoire.clonaldiversity.skip.boots) +
      ggtitle("CLONAL DIVERISTY: CDR3 AMINO ACID") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    Idents(data) <- "seurat_clusters"
    p37plots <- clonalOverlay(data,
                              reduction = "umap",
                              freq.cutpoint = screpertoire.clonaloverlay.freq.cutpoint,
                              bins = screpertoire.clonaloverlay.bins,
                              facet = screpertoire.clonaloverlay.facet) +
      guides(color = "none")+ scale_color_manual(values = cluster_colors)
    p37plots <- adjust_theme(p37plots, strip_size = 25)

    Idents(data) <- "CELL_TYPE"
    p38plots <- clonalOverlay(data,
                              reduction = "umap",
                              freq.cutpoint = screpertoire.clonaloverlay.freq.cutpoint,
                              bins = screpertoire.clonaloverlay.bins,
                              facet = screpertoire.clonaloverlay.facet) +
      guides(color = "none")+ scale_color_manual(values = celltype_colors)
    p38plots <- adjust_theme(p38plots, strip_size = 25)

    p39plots <- clonalNetwork(data,
                              reduction = "umap",
                              identity = screpertoire.clonalnetwork.identity,
                              filter.clones = screpertoire.clonalnetwork.filter.clones,
                              filter.identity = screpertoire.clonalnetwork.filter.identity,
                              filter.proportion = screpertoire.clonalnetwork.filter.proportion,
                              filter.graph = screpertoire.clonalnetwork.filter.graph,
                              cloneCall = screpertoire.clonalnetwork.cloneCall,
                              chain = screpertoire.clonalnetwork.chain,
                              exportTable = screpertoire.clonalnetwork.exportTable,
                              exportClones = screpertoire.clonalnetwork.exportClones) + scale_color_manual(values = cluster_colors)
    p39plots <- adjust_theme(p39plots)
  }

  #######################################################################################################################################
  somePDFPath = paste(cdir,"1URSA_PLOT_scIMMUNE_NUMBER_UNIQUE_CLONOTYPES_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=5,pointsize=12)
  print(p1plots)
  dev.off()

  somePDFPath = paste(cdir,"2URSA_PLOT_scIMMUNE_DISTR_CDR3_LENGTHS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=5,pointsize=12)
  print(p2plots)
  dev.off()

  write.csv(imm_pr, paste(cdir,"3URSA_TABLE_scIMMUNE_CLONE_PROPORTIONS_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)
  write.csv(imm_top, paste(cdir,"3URSA_TABLE_scIMMUNE_TOP_N_MOST_ABUNDANT_CLONOTYPES_PORPORTION_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)
  write.csv(imm_rare, paste(cdir,"3URSA_TABLE_scIMMUNE_RELATIVE_ABUNDANCE_RARE_CLONOTYPES_PORPORTION_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)
  write.csv(imm_hom, paste(cdir,"3URSA_TABLE_scIMMUNE_CLONALSPACE_HOMEOSTASIS_ALL_SAMPLES_ANALYSIS_",project_name,".csv", sep = ""), quote = F, row.names = T)

  somePDFPath = paste(cdir,"4URSA_PLOT_scIMMUNE_TOP_CLONAL_PROPORTIONS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  print(p4plots)
  dev.off()

  somePDFPath = paste(cdir,"5URSA_PLOT_scIMMUNE_RARE_CLONAL_PROPORTIONS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  print(p5plots)
  dev.off()

  somePDFPath = paste(cdir,"6URSA_PLOT_scIMMUNE_CLONALSPACE_HOMEOSTASIS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  print(p6plots)
  dev.off()

  somePDFPath = paste(cdir,"7URSA_PLOT_scIMMUNE_REPERTOIRE_OVERLAPS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  print(p7plots)
  dev.off()

  write.csv(all_gene_usage, paste(cdir,"8URSA_TABLE_scIMMUNE_GENE_USAGE_",project_name,".csv", sep = ""), quote = F, row.names = T)

  somePDFPath = paste(cdir,"9URSA_PLOT_scIMMUNE_GENE_USAGE_BY_CELL_TYPES_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7,pointsize=12)
  for(i in 1:length(p9plots)){
    plotx <- p9plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"10URSA_PLOT_scIMMUNE_GENE_USAGE_DIVERGENCE_CORRELATION_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  for(i in 1:length(p10plots)){
    plotx <- p10plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"11URSA_PLOT_scIMMUNE_CLUSTERING_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  for(i in 1:length(p11plots)){
    plotx <- p11plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"12URSA_PLOT_scIMMUNE_GENE_WISE_SAMPLE_DISTANCE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  for(i in 1:length(p12plots)){
    plotx <- p12plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"13URSA_PLOT_scIMMUNE_SPECTRATYPE_CLONALTYPE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  for(i in 1:length(p13plots)){
    plotx <- p13plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"14URSA_PLOT_scIMMUNE_DIVERSITY_ESTIMATES_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7,pointsize=12)
  for(i in 1:length(p14plots)){
    plotx <- p14plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"15URSA_PLOT_scIMMUNE_ALLUVIAL_TOP_10_MOST_ABUNDANT_CLONOTYPES_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=12.5,pointsize=12)
  for(i in 1:length(p15plots)){
    plotx <- p15plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"15URSA_PLOT_scIMMUNE_TOP_10_AA_KMER_SIZE_10_TO_20_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=9,pointsize=12)
  n <- 10
  for(i in 10:20){
    plotx <- contig_monarch[['data']]
    kmers <- getKmers(plotx, i, .col = immunarch.getkmers.col, .coding = immunarch.getkmers.coding)
    kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
    kmers <- kmers[order(unlist(kmers[,2]), decreasing = T),]
    p1 <- NULL
    p2 <- NULL
    p <- NULL
    p1 <- vis(kmers, .position = "stack", .head = n)
    p2 <- vis(kmers, .head = n, .position = "dodge")
    p <- p1/p2
    print(p)
  }
    dev.off()

    somePDFPath = paste(cdir,"15URSA_PLOT_scIMMUNE_AA_SEQUENCE_MOTIF_KMER_SIZE_",project_name,".pdf", sep = "")
    pdf(file=somePDFPath, width=12, height=6, pointsize=12)
    for(j in 1:length(contig_monarch)){
      plotx <- contig_monarch$data[[j]]
      kmers <- getKmers(plotx, i, .col = immunarch.getkmers.col, .coding = immunarch.getkmers.coding)
      kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
      kmers <- kmers[order(unlist(kmers[,2]), decreasing = T),]
      kmers_aa_stats <- kmer_profile(kmers,.method = screpertoire.kmer_profile..method, .remove.stop = screpertoire.kmer_profile..remove.stop)
      colnames(kmers_aa_stats) <- paste("KMER_POS_", 1:ncol(kmers_aa_stats), sep = "")

      p1 <- NULL
      p2 <- NULL
      p <- NULL
      p1 <- vis(kmers_aa_stats) + ggtitle(paste("POSITION FREQUENCY MATRIX: ", names(contig_monarch$data)[j], sep = ""))+
        scale_color_manual(values = gen_colors(color_conditions$tableau20,nrow(kmers_aa_stats)))+
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1,vjust = 1))
      p2 <- vis(kmers_aa_stats, .plot = "seq")+scale_fill_manual(values = gen_colors(color_conditions$tenx,nrow(kmers_aa_stats)))
      p <- p1+p2
      print(p)
    }
    dev.off()

  somePDFPath = paste(cdir,"16URSA_PLOT_scIMMUNE_CIRCOS_DIAGRAM_REPERTOIRE_OVERLAPS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=6, height=6, pointsize=12)
  vis(p16data, .plot = "circos", annotationTrack = c("grid", "axis"), preAllocateTracks = 1, grid.col = sample_colors, transparency = 0.2)
  title(paste(project_name,": Repertoire Overlaps", sep = ""), cex = 15)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, 0.6, CELL_META$sector.index,cex = 1,
                facing = "bending.outside", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA)
  dev.off()

  somePDFPath = paste(cdir,"17URSA_PLOT_scIMMUNE_HEATMAP_REPERTOIRE_OVERLAPS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5, pointsize=12)
  print(p17plots)
  dev.off()

  somePDFPath = paste(cdir,"18URSA_PLOT_scIMMUNE_CLONOTYPE_ABUNDANCE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3, pointsize=12)
  for(i in 1:length(p18plots)){
    plotx <- p18plots[[i]]
    print(plotx)
  }
  dev.off()

  for(i in 1:length(p19data)){
    x <- p19data[[i]]
    write.csv(x, paste(cdir,"19URSA_TABLE_scIMMUNE_CONTIG_ABUNDANCE_",names(p19data)[i],"_",project_name,".csv", sep = ""), quote = F, row.names = T)
  }

  somePDFPath = paste(cdir,"20URSA_PLOT_scIMMUNE_CDR3_CONTIG_LENGTH_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=10, pointsize=12)
  for(i in 1:length(p20plots)){
    plotx <- p20plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"21URSA_PLOT_scIMMUNE_ALLUVIAL_CLONOTYPES_PROPORTIONS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=16, height=16, pointsize=12)
  for(i in 1:length(p21plots)){
    plotx <- p21plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"22URSA_PLOT_scIMMUNE_CONTIG_V_GENE_USAGE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6, pointsize=12)
  for(i in 1:length(p22plots)){
    plotx <- p22plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"23URSA_PLOT_scIMMUNE_CLONAL_HOMEOSTASIS_PROPORTION_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3, pointsize=12)
  for(i in 1:length(p23plots)){
    plotx <- p23plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"24URSA_PLOT_scIMMUNE_MORISITA_INDEX_CLONAL_OVERLAP_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=8, pointsize=12)
  for(i in 1:length(p24plots)){
    plotx <- p24plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"25SCA_scImmun_HIERARCHICAL_CLUSTERING_JENSEN_SHANNON_DISTANCE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=10, pointsize=12)
  par(mar=c(3,4,1,6))
  for(i in 1:length(p25data)){
    plotx <- p25data[[i]]
    print(plot(plotx, horiz = TRUE))
  }
  dev.off()

  somePDFPath = paste(cdir,"26URSA_PLOT_scIMMUNE_CLONAL_SAMPLE_DIVERSITY_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7, pointsize=12)
  for(i in 1:length(p26plots)){
    plotx <- p26plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"27URSA_PLOT_scIMMUNE_UMAP_BY_SAMPLE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8, pointsize=12)
  print(p27plots)
  dev.off()

  somePDFPath = paste(cdir,"28URSA_PLOT_scIMMUNE_UMAP_BY_DATA_TYPE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8, pointsize=12)
  print(p28plots)
  dev.off()

  somePDFPath = paste(cdir,"29URSA_PLOT_scIMMUNE_UMAP_BY_CLUSTER_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8, pointsize=12)
  print(p29plots)
  dev.off()

  somePDFPath = paste(cdir,"30URSA_PLOT_scIMMUNE_UMAP_BY_GROUP_CLUSTER_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8, pointsize=12)
  print(p30plots)
  dev.off()

  somePDFPath = paste(cdir,"31URSA_PLOT_scIMMUNE_UMAP_BY_CELL_TYPE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8, pointsize=12)
  print(p31plots)
  dev.off()

  somePDFPath = paste(cdir,"32URSA_PLOT_scIMMUNE_UMAP_BY_GROUP_CELL_TYPE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=9, pointsize=12)
  print(p32plots)
  dev.off()

  somePDFPath = paste(cdir,"33URSA_PLOT_scIMMUNE_UMAP_BY_CLONOTYPE_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=6.7, pointsize=12)
  print(p33plots)
  dev.off()

  somePDFPath = paste(cdir,"34URSA_PLOT_scIMMUNE_CLONOTYPES_BY_CELL_TYPE_CLONECALL_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=6, pointsize=12)
  print(p34plots)
  dev.off()

  somePDFPath = paste(cdir,"35URSA_PLOT_scIMMUNE_CLONOTYPES_BY_CELL_TYPE_GROUP_DATA_TYPE_CLONECALL_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=18, height=10, pointsize=12)
  print(p35plots)
  dev.off()

  somePDFPath = paste(cdir,"36URSA_PLOT_scIMMUNE_CLONO_DIVERISITY_BY_CELL_TYPE_CLONECALL_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7, pointsize=12)
  for(i in 1:length(p36plots)){
    plotx <- p36plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"37URSA_PLOT_scIMMUNE_CLONAL_OVERLAY_CLUSTERS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=18, height=6.7, pointsize=12)
  print(p37plots)
  dev.off()

  somePDFPath = paste(cdir,"38URSA_PLOT_scIMMUNE_CLONAL_OVERLAY_CELL_TYPES_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=18, height=6.7, pointsize=12)
  print(p38plots)
  dev.off()

  somePDFPath = paste(cdir,"39URSA_PLOT_scIMMUNE_CLONAL_NETWORK_CLUSTERS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=14, height=8, pointsize=12)
  print(p39plots)
  dev.off()

  saveRDS(data, paste(cdir,"40URSA_DATA_scIMMUNE_INTEGRATED_DATA_",project_name,".RDS", sep = ""))
  print("Completed!")
}



