############################################################################################################
# Ursa R Package
# Phecda: CyTOF
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
#' @import akmedoids
#' @import Biobase
#' @import ConsensusClusterPlus
#' @import flowCore
#' @import FlowSOM
#' @import flowWorkspace
#' @import ggcyto
#' @import limma
#' @import openCyto
#' @import Rtsne
#' @import umap
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format CyTOF pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification CyTOF
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_CyTOF' by default.
#' @param input_dir Directory to all input files. Current working directory by
#' default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format
#' files.
#' @export
#'
CyTOFPip <- function(project_name = "Ursa_CyTOF",
                    input_dir = "./",
                    output_dir = "./",
                    pheno_file){
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "CYTOF", isDir = T)
  color_conditions <- color_ini()
  ctime <- time_ini()
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir)
  sample_files <- list.files(input_dir, recursive = T, full.names = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep(gsub(".*\\/(.*)","\\1",pheno_file, ignore.case = T),sample_files, ignore.case = T, invert = T)]
  input_dir <- gsub("(.*\\/).*$","\\1",sample_files[1], ignore.case = T)

  files_submit <- list.files(pattern = "*.fcs", path = input_dir, full.names = T, ignore.case = T, recursive = T)

  print("Reading FCS files..")

  data_submit <- NULL
  pheno_data$SID <- paste(pheno_data$SAMPLE_ID,pheno_data$GROUP, sep = "_")

  for(i in 1:nrow(pheno_data)){
    if (length(grep(pheno_data[i,grep("FILE", colnames(pheno_data),ignore.case = T)], files_submit, ignore.case = T)) > 0){
      current_name <- pheno_data[i,"SID"]
      current <- read.FCS(files_submit[grep(pheno_data[i,grep("FILE", colnames(pheno_data),ignore.case = T)], files_submit, ignore.case = T)], transformation = FALSE, truncate_max_range = FALSE)
      col_current <- as.character(parameters(current)$desc)
      if(length(which(is.na(col_current))) > 0){
        col_current[which(is.na(col_current))] <- as.character(parameters(current)$name)[which(is.na(col_current))]
      }
      current <- data.frame(FILE = current_name, flowCore::exprs(current))
      colnames(current) <- c("FILE", col_current)
      data_submit <- rbind(data_submit, current)
    }
  }

  sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SID)))
  names(sample_colors) <- unique(pheno_data$SID)

  colnames(data_submit) <- gsub("^X\\.","",colnames(data_submit), perl = TRUE)
  colnames(data_submit) <- gsub("^CR5","CXCR5",colnames(data_submit), perl = TRUE)
  data_submit$FILE <- gsub(".*?/","",data_submit$FILE, perl = TRUE)
  data_submit$FILE <- gsub("^\\s+|\\s+$|\\s+|\\s+$|\\.fcs","",data_submit$FILE, perl = TRUE)
  data_submit$FILE <- toupper(data_submit$FILE)
  colnames(data_submit) <- gsub("\\.|_|-|\\s+","_",colnames(data_submit), perl = TRUE)
  colnames(data_submit)[grep("FILE", colnames(data_submit), ignore.case = T, invert = T)] <- toupper(colnames(data_submit))[grep("FILE", colnames(data_submit), ignore.case = T, invert = T)]

  EXPR_PARA <- data_submit[,grep("FILE|time|cell_length|dna|Event_length|Center|Offset|Width|Residual", colnames(data_submit), ignore.case = T)]

  cofactor <- 5
  data <- asinh(data_submit[,grep("FILE|time|cell_length|dna|Event_length|Center|Offset|Width|Residual", colnames(data_submit), ignore.case = T, invert = T)] / cofactor)
  if(length(unique(pheno_data$BATCH)) > 1){
  data <- apply(data, 2, function(x){x <- lm(x ~ pheno_data[match(data_submit$FILE,pheno_data$SID),"BATCH"])$residual})
  }
  data <- data.frame(FILE = data_submit$FILE, data)

  for (i in grep("FILE",colnames(data), ignore.case = T, invert = T)){
    data[,i] <- as.numeric(as.character(data[,i]))
  }

  counts <- data.frame(plyr::count(data$FILE))
  colnames(counts) <- c("SID","CELL_COUNT")
  pheno_data <- plyr::join(pheno_data, counts, by = "SID")

  sample_list <- unique(data$FILE)
  sample_list
  colnames(data)
  marker_list <- as.character(colnames(data)[grep("FILE|population|tsne|SNEx|SNEy|PeakX|PeakY|time|length|DNA", colnames(data), ignore.case = T, invert = T)])

  ggdf <- pheno_data
  ggdf <- ggdf[order(ggdf[,grep("Individual.*ID", colnames(ggdf), ignore.case = T)], decreasing = T),]
  ggdf[,grep("Sample.*ID", colnames(ggdf), ignore.case = T)] <- factor(ggdf[,grep("Sample.*ID", colnames(ggdf), ignore.case = T)], levels = unique(c(ggdf[,grep("Sample.*ID", colnames(ggdf), ignore.case = T)])))

  individual_colors <- gen_colors(color_conditions$general, length(unique(pheno_data$INDIVIDUAL_ID)))
  names(individual_colors) <- unique(pheno_data$INDIVIDUAL_ID)

  print("Generating summary statistics..")

  p1plots <- ggplot(ggdf, aes(x = SAMPLE_ID, y = CELL_COUNT, fill = INDIVIDUAL_ID, label = CELL_COUNT)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_bw()+
    ylab("Cell Count")+ ggtitle(project_name)+
    scale_fill_manual(values = individual_colors) +
    guides(fill = guide_legend(reverse = TRUE))
  p1plots <- adjust_theme(p1plots)

  somePNGPath <- paste(cdir,"1URSA_PLOT_CyTOF_SAMPLE_CELL_SUMMARY_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = length(unique(ggdf$SAMPLE_ID))*200, units = "px", res = 300)
  print(p1plots)
  dev.off()

  ggdf <- data.frame(SAMPLE_ID = pheno_data[match(data$FILE, pheno_data$SID),"SAMPLE_ID"],
                     data[,grep("FILE", colnames(data), ignore.case = T, invert = T)])

  ggdf <- reshape2::melt(ggdf, id.var = "SAMPLE_ID", value.name = "Expression", variable.name = "Marker")
  ggdf$GROUP <- pheno_data[match(ggdf$SAMPLE_ID, pheno_data$SAMPLE_ID),"GROUP"]

  group_colors <- gen_colors(color_conditions$bright[c(1,3:(length(color_conditions$bright)),2)], length(unique(pheno_data$GROUP)))
  names(group_colors) <- unique(pheno_data$GROUP)

  p2plots <- ggplot(ggdf, aes(x = Expression, y = SAMPLE_ID, color = GROUP, group = SAMPLE_ID)) +
    geom_density_ridges(alpha = 0.7) +
    facet_wrap(~ Marker, ncol = 8, scales = "free") +
    theme_bw() + ylab("Density") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 11),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 14),
          strip.background = element_blank(),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.key = element_rect(size = 6),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(1, "cm")) +
    guides(color = guide_legend(ncol = 1)) +
    scale_color_manual(values = group_colors)

  somePNGPath <- paste(cdir,"2URSA_PLOT_CyTOF_MARKER_DENSITIES_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = (ceiling(length(unique(ggdf$Marker))/8)+length(unique(ggdf$SAMPLE_ID)))*200, units = "px", res = 300)
  print(p2plots)
  dev.off()

  files <- unique(data$FILE)
  files
  plot_median <- NULL
  expr <- NULL

  for(i in 1:length(files)){
    current <- data[which(data$FILE == files[i]),]
    current <- current[,grep("FILE|peak|tsne|snex|sney|population",colnames(current), invert = T, ignore.case = T)]
    colnames(current)
    for (j in 1:ncol(current)){
      current[,j] <- as.numeric(as.character(current[,j]))
      expr <- c(expr,median(current[,j])) # current[,j] > 0
    }
    plot_median <- rbind(plot_median, expr)
    expr <- NULL
  }

  plot_median <- data.frame(t(plot_median))
  row.names(plot_median) <- marker_list
  colnames(plot_median) <- toupper(gsub("\\.fcs","",files,ignore.case = T))

  mds <- plotMDS(plot_median, plot = FALSE)
  pca_out <- prcomp(t(plot_median), center = TRUE, scale. = FALSE)
  ggdf <- data.frame(SID = colnames(plot_median), MDS1 = mds$x, MDS2 = mds$y, PC1 = pca_out$x[,1], PC2 = pca_out$x[,2])
  ggdf <- plyr::join(ggdf, pheno_data, by = "SID")

  p3plots <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = GROUP, label = SAMPLE_ID)) +
    geom_point(aes(size = CELL_COUNT), alpha = 0.8) +
    geom_label_repel(show.legend = F) +
    scale_size_continuous(range = c(6, 12))+
    # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
    theme_bw() +
    scale_color_manual(values = group_colors)
  p3plots <- adjust_theme(p3plots)

  somePNGPath <- paste(cdir,"3URSA_PLOT_CyTOF_MDS_SAMPLE_DISTANCE_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 2500, units = "px", res = 300)
  print(p3plots)
  dev.off()

  p4plots <- ggplot(ggdf, aes(x = PC1, y = PC2, color = GROUP, label = SAMPLE_ID)) +
    geom_point(aes(size = CELL_COUNT), alpha = 0.8) +
    geom_label_repel(show.legend = F) +
    scale_size_continuous(range = c(6, 12))+
    # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
    theme_bw() +
    scale_color_manual(values = group_colors)
  p4plots <- adjust_theme(p4plots)

  somePNGPath <- paste(cdir,"4URSA_PLOT_CyTOF_PCA_SAMPLE_DISTANCE_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 2500, units = "px", res = 300)
  print(p4plots)
  dev.off()

  current <- plot_median
  current <- (scale((current)))
  current <- t(scale(t(current)))
  colnames(current) <- pheno_data[match(colnames(current), pheno_data$SID),"SAMPLE_ID"]

  p5plots <- complex_heatmap(as.matrix(current), col = jet2.col(n=100), legendtitle = "ZScore", col_title = project_name)

  somePNGPath <- paste(cdir,"5URSA_PLOT_CyTOF_HEATMAP_SAMPLE_SCALED_ASINH_MEDIAN_EXPRESSION_",project_name,".png", sep = "")
  png(somePNGPath, width=3000, height=2000, units = "px", res = 400)
  print(p5plots)
  dev.off()

  plotx <- melt(data.frame(Marker = row.names(plot_median), plot_median))
  colnames(plotx) <- c("MARKER","SID","ASINH_EXPRESSION")
  plotx <- cbind(plotx, pheno_data[match(gsub("\\.|\\s+","-",plotx$SID),gsub("\\.|\\s+","-",pheno_data$SID)),c("SAMPLE_ID","GROUP","BATCH")])

  plotx$MARKER <- factor(plotx$MARKER, levels = sort(unique(plotx$MARKER)))
  p6plots <- ggboxplot(plotx, x = "GROUP", y = "ASINH_EXPRESSION", fill = "GROUP", palette = "jco", add = "jitter", repel = TRUE)
  p6plots <- p6plots + facet_wrap(~ MARKER, scales = "free", ncol = 6) + stat_compare_means() + theme_classic() +
    scale_fill_manual(values = group_colors)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))+
    scale_color_manual(values = sample_colors)
  p6plots <- adjust_theme(p6plots)

  somePNGPath <- paste(cdir,"6URSA_PLOT_CyTOF_BOXPLOT_MARKERS_BY_GROUP_",project_name,".png", sep = "")
  png(somePNGPath, width=5000, height=ceiling(length(unique(plotx$MARKER))/6)*700, units = "px", res = 300)
  print(p6plots)
  dev.off()

  p7data <- as.dendrogram(hclust(as.dist(1-cor((current)))))

  somePNGPath <- paste(cdir,"7URSA_PLOT_CyTOF_DENDROGRAM_SAMPLES_",project_name,".png", sep = "")
  png(somePNGPath, width=3000, height=2000, units = "px", res = 300)
  par(mar=c(3,4,1,6))
  print(plot(p7data, horiz = TRUE))
  dev.off()

  NRS <- function(x, ncomp = 3){
    pr <- prcomp(x, center = TRUE, scale. = FALSE)
    score <- rowSums(outer(rep(1, ncol(x)),pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
    return(score)
  }

  nrs_sample <- NULL
  current <- data[,grep("FILE|SNE|PEak|population",colnames(data), ignore.case = T, invert = T)]
  for(i in 1:length(unique(pheno_data$SID))){
    nrs_sample <- rbind(nrs_sample, NRS(current[which(data$FILE == unique(pheno_data$SID)[i]),]))
  }

  rownames(nrs_sample) <- unique(pheno_data$SID)
  nrs_sample <- data.frame(nrs_sample)
  nrs <- colMeans(nrs_sample, na.rm = TRUE)
  markers_ord <- names(sort(nrs, decreasing = TRUE))
  nrs_sample$SID <- rownames(nrs_sample)
  ggdf <- melt(nrs_sample, id.var = "SID",
               value.name = "nrs", variable.name = "Markers")
  colnames(ggdf) <- c("SID","Markers","NRScore")
  ggdf$Markers <- factor(ggdf$Markers, levels = markers_ord)
  ggdf <- plyr::join(ggdf, pheno_data, by = "SID")

  p8plots <- ggplot(ggdf, aes(x = Markers, y = NRScore)) +
    # geom_point(aes(color = SAMPLE_ID), alpha = 0.8,
    # position = position_jitter(width = 0.3, height = 0)) +
    geom_boxplot(aes(fill = GROUP), alpha = 0.8, outlier.color = NA) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = sample_colors)+
    scale_fill_manual(values = group_colors)
  p8plots <- adjust_theme(p8plots, xangle = 45, xsize = 12, hejust = 1, vejust = 1)

  somePNGPath <- paste(cdir,"8URSA_PLOT_CyTOF_NRS_RANKING_ALL_MARKERS_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =2000, units = "px", res = 300)
  print(p8plots)
  dev.off()

  files <- unique(data$FILE)
  data0 <- NULL

  for(i in 1:length(files)){
    data0[[i]] <- data[which(data$FILE == files[i]),grep("FILE|SNE|PEak|population",colnames(data), ignore.case = T, invert = T)]
  }

  for (i in 1:length(files)){
    meta <- data.frame(name=colnames(data0[[i]]),desc=colnames(data0[[i]]))
    meta$range <- apply(apply(data0[[i]],2,range),2,diff)
    meta$minRange <- apply(data0[[i]],2,min)
    meta$maxRange <- apply(data0[[i]],2,max)
    data0[[i]] <- new("flowFrame",exprs=as.matrix(data0[[i]]),parameters=AnnotatedDataFrame(meta))
  }

  print("Dimension reduction and clustering..")

  data_fs = as(data0,"flowSet")
  flowCore::pData(data_fs)$name <- files
  som_input <- ReadInput(data_fs)
  set.seed(59)
  som <- BuildSOM(som_input)
  codes <- som$map$codes
  nmc <- 90
  somePDFPath <- paste(cdir,"9URSA_PLOT_CyTOF_CONSENSUS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=12,pointsize=12)
  mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 50,
                             pItem = 0.9, pFeature = 1, plot = NULL,
                             clusterAlg = "hc", innerLinkage = "complete", finalLinkage = "complete", distance = "euclidean", seed = 1234)
  dev.off()
  Kvec = 2:nmc
  x1 = 0.1; x2 = 0.9
  PAC = rep(NA,length(Kvec))
  names(PAC) = paste("K=",Kvec,sep="")
  for(i in Kvec){
    M = mc[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }
  optK = Kvec[which.min(PAC)]
  PAC

  PAC <- data.frame(K = as.numeric(as.character(gsub("K=","",names(PAC)))), PAC = as.numeric(as.character(PAC)))
  ck1 <- elbow_point(PAC$K, PAC$PAC)
  ck1 <- ceiling(ck1$x)
  PAC$Difference <- NULL
  PAC[1,"Difference"] <- NA
  for (i in 2:nrow(PAC)){
    PAC[i,"Difference"] <- PAC[i-1,"PAC"] - PAC[i,"PAC"]
  }

  # cPAC <- PAC[(max(which(PAC$PAC > quantile(PAC$PAC, 0.95)))+1):nrow(PAC),]
  DeltaY <- diff(PAC$PAC)
  PAC_turn <- which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
  ck2 <- min(PAC_turn[PAC_turn>5]) + 1
  if(ck1 != ck2){
    chosen_k <- min(ck1, ck2)
  }else{
    chosen_k <- ck1
  }

  code_clustering1 <- mc[[chosen_k]]$consensusClass
  cell_clustering1 <- code_clustering1[som$map$mapping[,1]]
  data$population <- cell_clustering1
  row.names(data) <- paste("Cell",1:nrow(data), sep = "")
  PAC$Label <- ifelse(PAC$K == chosen_k, paste("Chosen K = ", chosen_k, sep = ""), "")

  p9plots <- ggplot(PAC, aes(x= K, y= PAC, label = Label)) + geom_line(colour = "grey") +
    ggtitle("Chosen K Number of Clusters Based on PAC Method")+
    theme_classic(base_size = 20)+ geom_text_repel(
      max.overlaps = Inf,force=1,
      point.padding = 0, # additional padding around each point
      min.segment.length = 0, # draw all line segments
      max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
      box.padding = 0.3, size = 12, colour = "red")+
    geom_point(colour = ifelse(PAC$K == chosen_k, "red", "grey"), size = ifelse(PAC$K == chosen_k, 5, 2))
  p9plots <- adjust_theme(p9plots)

  somePNGPath <- paste(cdir,"9URSA_PLOT_CyTOF_PAC_ELBOW_PLOT_CHOSEN_CLUSTER_NUMBER_",chosen_k,"_CLUSTERS_",project_name,".png", sep = "")
  png(somePNGPath, width=4000, height=3000, res = 300)
  print(p9plots)
  dev.off()

  plot_median <- NULL
  expr <- NULL
  cell_number <- NULL
  pop_list <- 1:chosen_k
  data1 <- data.frame(FILE = data$FILE,
                      data[,grep("FILE|SNE|PEak|population",colnames(data), ignore.case = T, invert = T)],
                      population = cell_clustering1)

  for(i in 1:length(pop_list)){
    current <- data1[which(data1$population == pop_list[i]),grep("FILE|SNE|PEak|population",colnames(data1), ignore.case = T, invert = T)]
    cell_number <- c(cell_number,nrow(current))
    for (j in 1:ncol(current)){
      expr <- c(expr,median(current[,j])) # current[,j] > 0
    }
    plot_median <- rbind(plot_median, expr)
    expr <- NULL
  }

  row.names(plot_median) <- paste("Cluster_", pop_list, ":",cell_number,sep = "")
  colnames(plot_median) <- gsub("^.*?_","",marker_list)
  plot_median <- data.frame(plot_median)
  dista <- hclust(as.dist(1-cor(t(plot_median))), method = "complete")
  plot_median <- scale(plot_median)
  plot_median <- t(scale(t(plot_median)))
  plot_median <- as.matrix(plot_median)
  plot_median_csv <- plot_median[dista$order,]
  plot_median_csv <- plot_median_csv[,hclust(as.dist(1-cor((plot_median_csv))), method = "complete")$order]
  melt_plot_median <- melt(plot_median_csv)
  colnames(melt_plot_median) <- c("Cluster","Marker","Expression")
  melt_plot_median$Cluster_Order <- rep(1:nrow(plot_median_csv),ncol(plot_median_csv))
  melt_plot_median$Marker_Order <- rep(1:ncol(plot_median_csv),each = nrow(plot_median_csv))
  melt_plot_median$Size <- gsub(".*:(.*)","\\1",melt_plot_median$Cluster)
  melt_plot_median$Cluster <- gsub("Cluster_(.*):.*","\\1",melt_plot_median$Cluster)

  p10plots <- complex_heatmap(plot_median, col = jet2.col(n = 100), legendtitle = "ZScore", col_title = project_name)
  somePNGPath <- paste(cdir,"10URSA_PLOT_CyTOF_HEATMAP_CLUSTERS_ASINH_MEDIAN_",project_name,".png", sep = "")
  png(somePNGPath, width=ncol(plot_median)*120, height=nrow(plot_median)*300, units = "px", res = 400)
  print(p10plots)
  dev.off()

  plot_median <- NULL
  expr <- NULL
  cell_number <- NULL
  pop_list <- 1:chosen_k
  data1 <- data.frame(FILE = data$FILE,
                      data[,grep("FILE|SNE|PEak|population",colnames(data), ignore.case = T, invert = T)],
                      population = cell_clustering1)
  files <- unique(data$FILE)
  for(i in 1:length(pop_list)){
    for(j in 1:length(files)){
      current <- data1[which(data1$population == pop_list[i] & data1$FILE == files[j]),grep("FILE|SNE|PEak|population",colnames(data1), ignore.case = T, invert = T)]
      if(nrow(current) > 0){
        cell_number <- c(cell_number,nrow(current))
        for (k in 1:ncol(current)){
          expr <- c(expr,median(current[,k])) # current[,k] > 0
        }
        plot_median <- rbind(plot_median, data.frame(FILE = files[j], CLUSTER = pop_list[i], MARKER = colnames(current), expr))
        expr <- NULL
      }
    }
  }

  colnames(plot_median) <- c("SID","CLUSTER","MARKER","ASINH_EXPRESSION")
  plot_median <- cbind(plot_median, pheno_data[match(gsub("\\.|\\s+","-",plot_median$SID),gsub("\\.|\\s+","-",pheno_data$SID)),c("SAMPLE_ID","GROUP","BATCH")])

  p11plots <- ggplot(plot_median) +
    geom_boxplot(aes(x = MARKER, y = ASINH_EXPRESSION, fill = GROUP), position = position_dodge(),
                 alpha = 1,
                 outlier.color = NA) +
    geom_point(aes(x = MARKER, y = ASINH_EXPRESSION), alpha = 0.8, size = 0.2) +
    facet_wrap(~ CLUSTER, scales = "free", ncol = 2) +
    theme_classic() +
    scale_fill_manual(values = group_colors)+
    scale_color_manual(values = sample_colors)
  p11plots <- adjust_theme(p11plots, xangle = 45, xsize = 8, title_size = 15, hejust = 1, vejust = 1)

  somePNGPath <- paste(cdir,"11URSA_PLOT_CyTOF_BOXPLOT_CLUSTERS_BY_GROUP_",project_name,".png", sep = "")
  png(somePNGPath, width = ceiling(length(unique(plot_median$MARKER)))*120, height = ceiling(length(unique(plot_median$CLUSTER))/2)*700, units = "px", res = 300)
  print(p11plots)
  dev.off()

  print("Dimension reduction may take sometime depending on the size of the data and number of samples..")

  set.seed(59)
  pca_out <- prcomp(data[,grep("FILE|population|BARCODE|^CD45$", colnames(data), ignore.case = T, invert = T)], center = TRUE, scale. = FALSE)

  set.seed(59)
  umap_out <- umap::umap(data[,grep("FILE|population|BARCODE|^CD45$", colnames(data), ignore.case = T, invert = T)], pca = FALSE, perplexity = 30)

  plot_out <- data.frame(PC1 = pca_out$x[,1], PC2 = pca_out$x[,2],
                         UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2])

  plot_out$Cluster <- factor(data$population, levels = c(sort(as.numeric(as.character(unique(data$population))))))
  plot_out$FILE <- data$FILE
  plot_out$GROUP <- pheno_data[match(plot_out$FILE, pheno_data$SID),"GROUP"]

  cluster_colors <- gen_colors(color_conditions$tenx, length(unique(data$population)))
  names(cluster_colors) <- c(sort(as.numeric(as.character(unique(data$population)))))

  p12plots <- ggplot(plot_out,  aes(x = UMAP1, y = UMAP2, color = GROUP)) +
    geom_point(size = 0.95, alpha = 1) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    scale_color_manual(values = group_colors) +
    ggtitle(paste(project_name,": UMAP BY GROUP", sep = ""))
  p12plots <- adjust_theme(p12plots)

  somePNGPath <- paste(cdir,"12URSA_PLOT_CyTOF_UMAP_BY_GROUP_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =2500, units = "px", res = 300)
  print(p12plots)
  dev.off()

  p13plots <- own_facet_scatter(plot_out, feature1 = "UMAP1", feature2 = "UMAP2", isfacet = T, legend_pos = "none",
                                title = paste(project_name,": UMAP BY SAMPLE", sep = ""), col = sample_colors,ncol = 4,
                                color_by = "FILE", group_by = "FILE", xlabel = "UMAP1", ylabel = "UMAP2", strip_size = 12)
  p13plots <- adjust_theme(p13plots, legend = "none", strip_size = 12, xsize = 15)

  somePNGPath <- paste(cdir,"13URSA_PLOT_CyTOF_UMAP_BY_SAMPLE_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = ceiling(length(unique(plot_out$FILE))/4)*800, units = "px", res = 300)
  print(p13plots)
  dev.off()

  p14plots <- plot_bygroup(plot_out, x = "UMAP1", y = "UMAP2", group = "Cluster",
                           plot_title = paste(project_name,": UMAP BY CLUSTER", sep = ""),
                           col = cluster_colors, numeric = T, annot = T, point_size = 0.8, legendsize = 20, label_size = 8)
  p14plots <- adjust_theme(p14plots)

  somePNGPath <- paste(cdir,"14URSA_PLOT_CyTOF_UMAP_BY_CLUSTER_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 2500, units = "px", res = 300)
  print(p14plots)
  dev.off()

  p15plots <- plot_bygroup(plot_out, x = "PC1", y = "PC2", group = "Cluster",
                           plot_title = paste(project_name,": PCA BY CLUSTER", sep = ""),
                           col = cluster_colors, numeric = T, annot = F, point_size = 0.8, legendsize = 20, label_size = 8)
  p15plots <- adjust_theme(p15plots)
  somePNGPath <- paste(cdir,"15URSA_PLOT_CyTOF_PCA_BY_CLUSTER_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 2500, units = "px", res = 300)
  print(p15plots)
  dev.off()

  marker_list <- unique(colnames(data[,grep("FILE|population", colnames(data), ignore.case = T, invert = T)]))
  melt_data <- melt(data.frame(data, UMAP1 = plot_out$UMAP1, UMAP2 = plot_out$UMAP2),
                    measure.vars = colnames(data)[grep("FILE|population", colnames(data), ignore.case = T, invert = T)],
                    id.vars = c("UMAP1","UMAP2", "FILE"))
  colnames(melt_data) <- c("UMAP1","UMAP2","FILE","Marker","Asinh_Expression")
  melt_data$GROUP <- pheno_data[match(melt_data$FILE, pheno_data$SID),"GROUP"]
  p16data <- melt_data

  for(i in 1:length(group_colors)){
    plotx <- p16data[which(p16data$GROUP == names(group_colors)[i]),]
    p <- NULL
    p <- ggplot(plotx,  aes(x = UMAP1, y = UMAP2, color = Asinh_Expression)) +
      facet_wrap(~Marker, ncol = 6) +
      geom_point(size = 0.8) + theme_classic() + ggtitle(names(group_colors)[i]) +
      scale_color_gradientn(colours = colorRampPalette(color_conditions$BlueYellowRed)(100))
    p <- adjust_theme(p, xsize = 10)
    somePNGPath <- paste(cdir,"16URSA_PLOT_CyTOF_UMAP_ASINH_EXPRESSION_MARKER_",names(group_colors)[i],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
    print(p)
    dev.off()
  }

  data$population <- cell_clustering1
  sample_list <- sort(as.character(unique(data$FILE)))
  sample_list
  summary <- data.frame(population = sort(unique(data$population)))

  for (i in 1:length(sample_list)){
    current <- data[which(data$FILE == sample_list[i]),]
    summary <- plyr::join(summary,plyr::count(current, c("FILE","population"))[c("population","freq")],by = "population", type = "left")
    colnames(summary)[which(colnames(summary) == "freq")] <- unique(current$FILE)
  }
  summary[is.na(summary)] <- 0

  colnames(summary)[which(colnames(summary) == "population")] <- "Cluster"
  summary_percent = data.frame(summary)
  summary_percent$Total_Cells_Number <- rowSums(summary_percent[,which(!colnames(summary) == "Cluster")])

  disease_types <- as.character(unique(pheno_data$GROUP))

  for (i in 1:length(disease_types)){
    mm <- c("Cluster",as.character(pheno_data[match(colnames(summary)[grep("Cluster",colnames(summary), invert = T, ignore.case = T)],as.character(pheno_data$SID)),"GROUP"]))
    if(length(grep(as.character(disease_types[i]),mm, ignore.case = T)) > 1){
      current_rs <- rowSums(summary[,grep(as.character(disease_types[i]),mm, ignore.case = T)])
    } else{
      current_rs <- as.numeric(as.character(summary[,grep(as.character(disease_types[i]),mm, ignore.case = T)]))
    }
    summary_percent[,ncol(summary_percent)+1] <- current_rs
    colnames(summary_percent)[ncol(summary_percent)] <- paste("GROUP_",disease_types[i],"_Total_Cells", sep = "")
    summary_percent[,ncol(summary_percent)+1] <- (current_rs/summary_percent[,"Total_Cells_Number"])*100
    colnames(summary_percent)[ncol(summary_percent)] <- paste("GROUP_",disease_types[i],"_Percent", sep = "")
  }

  for(l in grep("Cluster|Total_Cells|Strata_|Percent",colnames(summary_percent), ignore.case = TRUE, invert = TRUE)){
    summary_percent[,l] <- summary_percent[,l]/summary_percent[,"Total_Cells_Number"]
  }

  melt_summary <- melt(summary_percent[,grep("Percent|Total_Cells|Strata_",colnames(summary_percent), ignore.case = TRUE, invert = T)], id.vars = c("Cluster"))
  colnames(melt_summary) <- c("Cluster","SID","Proportion")
  melt_summary <- cbind(melt_summary, pheno_data[match(gsub("\\.|\\s+","\\-",melt_summary$SID), gsub("\\.|\\s+","\\-",pheno_data$SID)),grep("SID", colnames(pheno_data), ignore.case = T, invert = T)])
  melt_summary$Cluster <- factor(melt_summary$Cluster, levels = sort(unique(as.numeric(as.character(melt_summary$Cluster)))))

  p17plots <- ggplot(melt_summary, aes(Cluster, Proportion, fill = SAMPLE_ID))+
    geom_bar(stat="identity", alpha=0.8)+
    coord_polar()+
    scale_fill_viridis(discrete=TRUE)+
    ggtitle(paste(project_name,"\nFrequency of Samples in Each Cluster", sep = ""))+
    theme_bw(base_size = 28)+
    theme(plot.margin = unit(c(3,3,3,3), "cm"),
          plot.title = element_text(size=25, face = "bold", hjust = 0.5),
          legend.title=element_text(size=20, face = "bold"),
          legend.text=element_text(size=18),
          legend.key.size = unit(2, 'lines'),
          axis.title.x = element_text(colour="black", size = 20, face = "bold", vjust = -10),
          axis.title.y = element_text(colour="black", size = 20, face = "bold", vjust = 10),
          strip.text = element_text(size = 20, face = "bold"),
          axis.text.x=element_text(colour="black", size = 20),
          axis.text.y=element_text(colour="black", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())

  somePNGPath <- paste(cdir,"17URSA_PLOT_CyTOF_PROPORTION_OF_SAMPLES_IN_EACH_CLUSTER_",project_name,".png", sep = "")
  png(somePNGPath, width = 5000, height = 4800, units = "px", res = 300)
  print(p17plots)
  dev.off()

  p18plots <- ggplot(melt_summary,aes(SAMPLE_ID,Proportion,fill=Cluster))+
    geom_bar(stat = "identity", position = "fill", color = "black", size = 1)+ #, color = "black", size = 1
    theme_classic()+
    scale_fill_manual(values = cluster_colors)+
    guides(color=guide_legend(title="CLUSTERS", ncol = 2),fill = guide_legend(title="CLUSTERS", ncol = 2))
  p18plots <- adjust_theme(p18plots, xangle = 45, hejust = 1, vejust = 1, strip_size = 1)

  somePNGPath <- paste(cdir,"18URSA_PLOT_CyTOF_PROPORTION_OF_CLUSTERS_IN_EACH_SAMPLE_",project_name,".png", sep = "")
  png(somePNGPath, width = length(unique(melt_summary$SAMPLE_ID))*400, height = 2500, units = "px", res = 300)
  print(p18plots)
  dev.off()

  code_sizes <- table(factor(som$map$mapping[,1], levels = 1:nrow(codes)))
  code_sizes <- as.numeric(code_sizes)
  set.seed(59)
  umap_out <- umap(codes)
  tsne_out <- Rtsne(codes)
  pca_out <- prcomp(codes, center = TRUE, scale. = FALSE)
  codes_dr <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2],
                         tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                         PC1 = pca_out$x[, 1], PC2 = pca_out$x[, 2])
  codes_dr$Cluster <- factor(code_clustering1)
  codes_dr$Size <- code_sizes
  codes_dr$SOM_ID <- as.numeric(as.character(row.names(codes_dr)))

  som_data <- data
  som_data$population <- som$map$mapping[,1]
  som_result <- data.frame(SOM_ID = sort(unique(som_data$population)))
  som_result$SOM_Phenotype <- "UNDEFINED"

  # choose top n to display
  top_n_markers <- 5
  # markers_to_note <- c("CD4","CD8","CD56","CD45RA","CD45RO")

  current0 = apply(som_data[,grep(paste(marker_list, collapse = "|"), colnames(som_data), ignore.case = TRUE)], 2, list)
  current0 = lapply(current0, function(x){x = (as.numeric(unlist(x))[as.numeric(unlist(x))>0])})
  median_c0 <- as.numeric(lapply(current0, median))

  som_result$Cluster <- codes_dr[match(som_result$SOM_ID, codes_dr$SOM_ID),"code_clustering1"]

  node_result <- data.frame(Cluster_ID = sort(unique(data$population)))
  node_result$Cluster_Phenotype <- "UNDEFINED"

  for(i in 1:nrow(node_result)){
    node_name3 <- NULL
    current_node <- as.numeric(as.character(node_result[i,"Cluster_ID"]))
    current3 = apply(data[data$population==current_node,grep(paste(marker_list, collapse = "|"), colnames(data), ignore.case = TRUE)], 2, list)
    current3 = lapply(current3, function(x){x = (as.numeric(unlist(x))[as.numeric(unlist(x))>1])})
    length_c3 = as.numeric(unlist(lapply(current3, length)))
    median_c3 <- as.numeric(lapply(current3, median))
    median_c3c0 <- (median_c3+1)/(median_c0+1)
    node_list3 = marker_list[ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE)][order(median_c3c0[which(ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE))], decreasing = T)]
    if (length(node_list3) == 0){
      node_name3 = "UNDEFINED"
    }else{
      if(length(node_list3) > top_n_markers){node_list3 <- node_list3[1:top_n_markers]}
      node_name3 <- paste(node_name3,paste(node_list3, collapse = "+"), "+", sep = "")

    }
    node_result[i,"Cluster_Phenotype"] <- node_name3
    node_list3 <- NULL

  }

  codes_dr$Cluster_Phenotype <- node_result[match(codes_dr$Cluster, node_result$Cluster_ID),"Cluster_Phenotype"]
  colnames(codes_dr)[which(colnames(codes_dr) == "Cluster")] <- "Cluster"

  for(i in 1:nrow(som_result)){
    node_name3 <- NULL
    current_node <- as.numeric(as.character(som_result[i,"SOM_ID"]))
    current3 = apply(som_data[som_data$population==current_node,grep(paste(marker_list, collapse = "|"), colnames(som_data), ignore.case = TRUE)], 2, list)
    current3 = lapply(current3, function(x){x = (as.numeric(unlist(x))[as.numeric(unlist(x))>0])})
    length_c3 = as.numeric(unlist(lapply(current3, length)))
    median_c3 <- as.numeric(lapply(current3, median))
    median_c3c0 <- (median_c3+1)/(median_c0+1)
    node_list3 = marker_list[ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE)][order(median_c3c0[which(ifelse(median_c3c0 > 1.25 & length_c3 > 3 & is.na(median_c3c0) == FALSE, TRUE, FALSE))], decreasing = T)]
    if (length(node_list3) == 0){
      node_name3 = "UNDEFINED"
    }else{
      if(length(node_list3) > top_n_markers){node_list3 <- node_list3[1:top_n_markers]}
      node_name3 <- paste(node_name3,paste(node_list3, collapse = "+"), "+", sep = "")

    }
    som_result[i,"SOM_Phenotype"] <- node_name3
    node_list3 <- NULL

  }

  codes_dr$SOM_Phenotype <- som_result[match(codes_dr$SOM_ID, som_result$SOM_ID),"SOM_Phenotype"]
  # For user to download their SOM Phenotypes
  # codes_dr_simp <- codes_dr[,c("SOM_ID","SOM_Phenotype","Cluster","Cluster_Phenotype","Size")]
  # colnames(codes_dr_simp) <- c("SOM_ID","SOM_Phenotype","NODE_ID","NODE_Phenotype","Cluster_Size")
  codes_dr$SOM_ID <- factor(codes_dr$SOM_ID, levels = sort(unique(as.numeric(as.character(codes_dr$SOM_ID)))))
  colnames(codes_dr)[grep("code_clustering", colnames(codes_dr), ignore.case = T)] <- "Cluster"
  codes_dr$Cluster <- factor(codes_dr$Cluster, levels = sort(unique(as.numeric(as.character(codes_dr$Cluster)))))
  p20data <- codes_dr[,grep("UMAP|tSNE|PCA|Label", colnames(codes_dr), ignore.case = T, invert = T)]
  write.csv(p20data, paste(cdir, "20URSA_PLOT_TABLE_DATA_CLUSTER_AUTO_PHENOTYPES_PREDICTION_",project_name,".csv",sep = ""),row.names = F)

  codes_dr$Label <- paste("SOM:",codes_dr$SOM_ID, sep = ":")
  p19plots <- ggplot(codes_dr, aes(x = tSNE1, y = tSNE2,color = Cluster, size = Size, label = Label)) +
    geom_point(alpha = 0.7) + theme_classic() +
    scale_point_size_continuous(range = c(9,18))+
    geom_label_repel(box.padding = 0.2, max.overlaps = Inf, size = 2.2, show.legend = F)+
    scale_color_manual(values = cluster_colors) +
    theme(legend.position = "right",
          plot.title = element_text(size=18, face = "bold", hjust = 0.5))+
    ggtitle(paste(project_name,"\nSOM Nodes tSNE", sep = ""))
  p19plots <- adjust_theme(p19plots)

  somePNGPath <- paste(cdir,"19URSA_PLOT_CyTOF_tSNE_SOM_NODES_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 400)
  print(p19plots)
  dev.off()

  p21plots <- ggplot(codes_dr, aes(x = UMAP1, y = UMAP2,color = Cluster, size = Size, label = Label)) +
    geom_point(alpha = 0.7) + theme_classic() +
    scale_point_size_continuous(range = c(9,18))+
    geom_label_repel(box.padding = 0.2, max.overlaps = Inf, size = 2.2, show.legend = F)+
    scale_color_manual(values = cluster_colors) +
    theme(legend.position = "right",
          plot.title = element_text(size=18, face = "bold", hjust = 0.5))+
    ggtitle(paste(project_name,"\nSOM Nodes UMAP", sep = ""))
  p21plots <- adjust_theme(p21plots)

  somePNGPath <- paste(cdir,"21URSA_PLOT_CyTOF_UMAP_SOM_NODES_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 400)
  print(p21plots)
  dev.off()

  p22plots <- ggplot(codes_dr, aes(x = PC1, y = PC2, color = Cluster, size = Size, label = Label)) +
    geom_point(alpha = 0.7) + theme_classic() +
    scale_point_size_continuous(range = c(9,18))+
    geom_label_repel(box.padding = 0.2, max.overlaps = Inf, size = 2.2, show.legend = F)+
    scale_color_manual(values = cluster_colors) +
    theme(legend.position = "right",
          plot.title = element_text(size=18, face = "bold", hjust = 0.5))+
    ggtitle(paste(project_name,"\nSOM Nodes PCA", sep = ""))
  p22plots <- adjust_theme(p22plots)

  somePNGPath <- paste(cdir,"22URSA_PLOT_CyTOF_PCA_SOM_NODES_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 400)
  print(p22plots)
  dev.off()

  saveRDS(data, paste(cdir,"23URSA_DATA_CyTOF_",project_name,".RDS", sep = ""))
  print("Completed!")

}
