############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Merak: Flow
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
#' @import CytoExploreR
#' @import CytoExploreRData
#' @import flowCore
#' @import FlowSOM
#' @import flowWorkspace
#' @import ggcyto
#' @import limma
#' @import openCyto
#' @import robustbase
#' @import Seurat
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format Flow Cytometry pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification Flow
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_Flow' by default.
#' @param input_dir Directory to all input files. Current working directory by
#' default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format
#' files.
#' @export
#'
FlowPip <- function(project_name = "Ursa_Flow",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file){
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "FLOW", isDir = T)
  color_conditions <- color_ini()
  sample_colors <- gen_colors(color_conditions$monet, length(unique(pheno_data$SAMPLE_ID)))
  names(sample_colors) <- unique(pheno_data$SAMPLE_ID)
  ctime <- time_ini()
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir)
  sample_files <- list.files(input_dir, recursive = T, full.names = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep(gsub(".*\\/(.*)","\\1",pheno_file, ignore.case = T),sample_files, ignore.case = T, invert = T)]
  input_dir <- gsub("(.*\\/).*$","\\1",sample_files[1], ignore.case = T)

  fs_data <- read.flowSet(path = input_dir, pattern = ".fcs", alter.names = TRUE, transformation = FALSE) # linearize + scale to [0,1]
  sampleNames(fs_data) <- pheno_data[match(sampleNames(fs_data), pheno_data$FILE),grep("SAMPLE.*ID", colnames(pheno_data), ignore.case = T)]
  phenoData(fs_data)$name <- pheno_data[match(phenoData(fs_data)$name, pheno_data$FILE),grep("SAMPLE.*ID", colnames(pheno_data), ignore.case = T)]
  phenoData(fs_data)$GROUP <- pheno_data[match(sampleNames(fs_data), pheno_data$SAMPLE_ID),"GROUP"]
  phenoData(fs_data)$BATCH <- pheno_data[match(sampleNames(fs_data), pheno_data$SAMPLE_ID),"BATCH"]

  data <- NULL
  for(i in 1:length(fs_data)){
    data <- rbind(data,data.frame(SAMPLE_ID = sampleNames(fs_data)[i], fs_data[[i]]@exprs))
  }

  current_names <- fs_data[[1]]@parameters@data[match(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],
                                                      fs_data[[1]]@parameters@data$name),"desc"]
  colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)] <- paste(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],ifelse(is.na(current_names), "", current_names), sep = "_")
  colnames(data) <- gsub("(.*)_$","\\1",colnames(data))
  data <- melt(data)
  colnames(data) <- c("SAMPLE_ID","CHANNEL","ASINH_COUNT")
  cofactor <- 150
  data$ASINH_COUNT <- asinh(data$ASINH_COUNT/cofactor)
  pheno_data <- pheno_data[order(pheno_data$GROUP, decreasing = F),]
  data$GROUP <- pheno_data[match(data$SAMPLE_ID,pheno_data$SAMPLE_ID),"GROUP"]
  data$SAMPLE_ID <- factor(data$SAMPLE_ID, unique(pheno_data$SAMPLE_ID))

  p1plots <- NULL
  channels <- as.character(unique(data$CHANNEL))
  somePDFPath <- paste(cdir,"1URSA_PLOT_ASINH_DENSITY_PLOT_CHANNELS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  for(j in 1:length(channels)){
    current <- data[data$CHANNEL == channels[j],]
    p1plots[[j]] <- ggplot(data[data$CHANNEL == channels[j],],
                           aes(x = ASINH_COUNT, y = SAMPLE_ID, fill = SAMPLE_ID, color = SAMPLE_ID)) +
      theme_classic()+
      geom_density_ridges(alpha = 0.8) +
      ggtitle(paste("Channel: ", channels[j], sep = "")) +
      scale_fill_manual(values = sample_colors) +
      scale_color_manual(values = sample_colors)
    p1plots[[j]] <- adjust_theme(p1plots[[j]])
    names(p1plots)[j] <- channels[j]

    # somePNGPath <- paste(cdir,"1URSA_PLOT_FLOW_ASINH_DENSITY_PLOT_CHANNEL_",channels[[j]],"_",project_name,".png", sep = "")
    # png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
    print(p1plots[[j]])
    # dev.off()
  }
  dev.off()

  print("Running automated gating..")
  gs <- NULL
  lastGate <- NULL
  p2plots <- NULL
  data_current <- NULL
  data_summary <- NULL
  k <- 1
  somePDFPath <- paste(cdir,"2URSA_PLOT_FLOW_GATING_STRATEGY_PART1_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  for(i in 1:length(fs_data)){
    cname <- sampleNames(fs_data)[i]
    print(paste("Running gating for ", cname, "..", sep = ""))
    gs[[i]] <- GatingSet(fs_data[i])
    gs[[i]] <- cyto_transform(gs[[i]], type = "arcsinh")
    if(length(colnames(fs_data[[i]])[grep("FSC.*A|SSC.*A", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnl <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
      chnlx <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("SSC.*A", colnames(fs_data), ignore.case = T)]
      local_min <- density(fs_data[[i]]@exprs[,chnl])$x[which(diff(sign(diff(density(fs_data[[i]]@exprs[,chnl])$y)))==+2)+1]
      current_ff <- gh_pop_get_data(gs[[i]])
      g <- openCyto:::.mindensity(current_ff, channels = chnl, filterId = "Exclude Debris", gate_range=c(local_min[1]-500,local_min[2]-1000))
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g)
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "root")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 0.7) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 1: Exclude Debris",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(flowCore::exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(flowCore::exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      print("after creating ggplot2 object")
      lastGate <- "Exclude Debris"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      # somePNGPath <- paste(cdir,"2URSA_PLOT_FLOW_GATING_STRATEGY_PART1_",chnlx,"_VS_",chnly,"_",cname,".png", sep = "")
      # png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
      print("before plot")
      print(p2plots[[k]])
      print("after plot")

      # dev.off()
      k <- k + 1
    }

    if(length(colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnl <- colnames(fs_data[[i]])[grep("FSC.*A|FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]
      chnlx <- colnames(fs_data[[i]])[grep("FSC.*A", colnames(fs_data[[i]]), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]

      g <- openCyto:::.singletGate(fs_data[[i]], channels = chnl, filterId = "Exclude Doublets or Multiplets")
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g, parent = "Exclude Debris")
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "Exclude Debris")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 1.5) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 2: Exclude Doublets or Multiplets",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(flowCore::exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(flowCore::exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Exclude Doublets or Multiplets"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      # somePNGPath <- paste(cdir,"2URSA_PLOT_FLOW_GATING_STRATEGY_PART2_",chnlx,"_VS_",chnly,"_",cname,".png", sep = "")
      # png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
      print(p2plots[[k]])
      # dev.off()
      k <- k + 1
    }

    if(length(colnames(fs_data[[i]])[grep("SSC.*W|SSC.*H", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnlx <- colnames(fs_data[[i]])[grep("SSC.*W", colnames(fs_data[[i]]), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("SSC.*H", colnames(fs_data[[i]]), ignore.case = T)]

      g <- gate_flowclust_2d(fs_data[[i]], xChannel = chnlx, yChannel = chnly, K = 3, filterId = "Single Cells Gate")
      g@filterId <- "Single Cells Gate"
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g, parent = "Exclude Doublets or Multiplets")
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "Exclude Doublets or Multiplets")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 1.5) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 2: Exclude Doublets or Multiplets",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(flowCore::exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(flowCore::exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Single Cells Gate"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      # somePNGPath <- paste(cdir,"2URSA_PLOT_FLOW_GATING_STRATEGY_PART3_",chnlx,"_VS_",chnly,"_",cname,".png", sep = "")
      # png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
      print(p2plots[[k]])
      # dev.off()
      k <- k + 1
    }

    if(length(colnames(fs_data[[i]])[grep("FSC.*H|FSC.*W", colnames(fs_data[[i]]), ignore.case = T)]) == 2){
      chnlx <- colnames(fs_data[[i]])[grep("FSC.*W", colnames(fs_data[[i]]), ignore.case = T)]
      chnly <- colnames(fs_data[[i]])[grep("FSC.*H", colnames(fs_data[[i]]), ignore.case = T)]

      g <- gate_flowclust_2d(fs_data[[i]], xChannel = chnlx, yChannel = chnly, K = 3, filterId = "Single Cells Gate 2")
      g@filterId <- "Single Cells Gate 2"
      mylimits <- ggcyto_par_set(limits = "instrument")
      gs_pop_add(gs[[i]][[1]], g, parent = "Single Cells Gate")
      recompute(gs[[i]][[1]])
      plotx <- gh_pop_get_data(gs[[i]], "Single Cells Gate")
      p2plots[[k]] <- autoplot(plotx, x = chnlx,y = chnly, bin = 300, strip.text = "gate") +
        geom_gate(g, colour = "red", size = 1.5) +
        geom_stats(size = 8,adjust = 0.5, label.padding = unit(0.05, "lines"), digits = 4) +
        theme_classic(base_size = 20) + ggtitle(paste("STEP 4: Single Cells Gate 2",sep = ""))+
        scale_fill_gradientn(colours = color_conditions$flow)+
        ggcyto_par_set(limits = list(x = c(-10,max(flowCore::exprs(current_ff)[,chnlx])),
                                     y = c(-10, max(flowCore::exprs(current_ff)[,chnly]))))+
        theme(strip.text = element_text(size = 20, face = "bold"),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      lastGate <- "Single Cells Gate 2"
      names(p2plots)[k] <- paste(cname,lastGate, sep = ":")
      # somePNGPath <- paste(cdir,"2URSA_PLOT_FLOW_GATING_STRATEGY_PART4_",chnlx,"_VS_",chnly,"_",cname,".png", sep = "")
      # png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
      print(p2plots[[k]])
      # dev.off()
      k <- k + 1
    }

    data_current[[i]] <- gh_pop_get_data(gs[[i]][[1]], y = lastGate, inverse.transform = T)
    ps <- data.frame(gs_pop_get_count_fast(gs[[i]][[1]]))
    ps$Percent_of_parent <- ps$Count/ps$ParentCount
    data_summary <- rbind(data_summary, ps)
  }
  dev.off()

  data_summary$name <- factor(data_summary$name, levels = sort(unique(data_summary$name)))
  write.csv(data_summary, paste(cdir, "URSA_TABLE_FLOW_FILTERING_SUMMARY_",project_name,".csv",sep = ""),row.names = F)

  gating_colors <- gen_colors(color_conditions$manycolors, length(unique(data_summary$Population)))
  names(gating_colors) <- unique(data_summary$Population)

  group_colors <- gen_colors(color_conditions$general, length(unique(pheno_data$GROUP)))
  names(group_colors) <- unique(pheno_data$GROUP)

  data_current <- as(data_current,"flowSet")
  sampleNames(data_current) <- sampleNames(fs_data)
  flowCore::pData(data_current)$name <- sampleNames(fs_data)

  data <- NULL
  for(i in 1:length(data_current)){
    data <- rbind(data,data.frame(SAMPLE_ID = flowCore::pData(data_current)$name[i], asinh(flowCore::exprs(data_current[[i]])/cofactor)))
  }

  current_names <- fs_data[[1]]@parameters@data[match(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],fs_data[[1]]@parameters@data$name),"desc"]
  colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)] <- paste(colnames(data)[grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)],ifelse(is.na(current_names), "", current_names), sep = "_")
  colnames(data) <- gsub("(.*)_$","\\1",colnames(data))

  melt_data <- melt(data)
  colnames(melt_data) <- c("SAMPLE_ID","CHANNEL","ASINH_COUNT")
  # cofactor <- 5
  # melt_data$ASINH_COUNT <- asinh(melt_data$ASINH_COUNT/cofactor)
  melt_data$GROUP <- pheno_data[match(melt_data$SAMPLE_ID,pheno_data$SAMPLE_ID),"GROUP"]
  melt_data$SAMPLE_ID <- factor(melt_data$SAMPLE_ID, unique(pheno_data$SAMPLE_ID))

  channels <- as.character(unique(melt_data$CHANNEL))
  p3plots <- NULL
  somePDFPath <- paste(cdir,"3URSA_PLOT_FLOW_POST_FILTERING_ASINH_DENSITY_PLOT_CHANNELS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  for(j in 1:length(channels)){
    p3plots[[j]] <- ggplot(melt_data[melt_data$CHANNEL == channels[j],], aes(x = ASINH_COUNT, y = SAMPLE_ID, fill = SAMPLE_ID, color = SAMPLE_ID)) +
      geom_density_ridges(alpha = 0.8) +
      theme_classic()+
      ggtitle(paste("Channel: ", channels[j], sep = "")) +
      scale_fill_manual(values = sample_colors) +
      scale_color_manual(values = sample_colors)
    p3plots[[j]] <- adjust_theme(p3plots[[j]])
    names(p3plots)[j] <- channels[j]
    # somePNGPath <- paste(cdir,"3URSA_PLOT_FLOW_POST_FILTERING_ASINH_DENSITY_PLOT_CHANNEL_",channels[[j]],"_",project_name,".png", sep = "")
    # png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
    print(p1plots[[j]])
  }
  dev.off()

  p4plots <- ggplot(data_summary, aes(x = name, y = Count, fill = Population, group = Population)) +
    geom_bar(stat="identity",alpha=0.9, size=0.5, colour = "black") +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab("Cell Count")+xlab("Samples") +
    scale_fill_manual(values = gating_colors)
  p4plots <- adjust_theme(p4plots, legend = "bottom") + guides(fill=guide_legend(ncol=1))

  somePDFPath <- paste(cdir,"4URSA_PLOT_FLOW_CELL_COUNT_SUMMARY_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p4plots)
  dev.off()

  current <- split(data_summary, data_summary$name)
  current <- lapply(current, function(x){
    current <- data.frame(name = unique(x$name),Total_Cell = max(x$Count),
                          Filtered_Cell_Count = min(x$Count))
  })
  current <- do.call(rbind.data.frame, current)
  pheno_data$Total_Cell <- current[match(pheno_data$SAMPLE_ID, current$name),"Total_Cell"]
  pheno_data$Filtered_Cell_Count <- current[match(pheno_data$SAMPLE_ID, current$name),"Filtered_Cell_Count"]

  # if(is.null(selected_markers)){
  data <- data[,grep("FSC.*W|FSC.*A|FSC.*H|SSC.*A|SSC.*H|SSC.*W|Time|Hoechst",
                     colnames(data), ignore.case = T, invert = T)]
  data <- data[,grep("_",colnames(data), ignore.case = T)]

  print("Running samplewise median expression and distance..")

  files <- unique(data$SAMPLE_ID)
  files
  plot_median <- NULL
  expr <- NULL

  for(i in 1:length(files)){
    current <- data[which(data$SAMPLE_ID == files[i]),]
    current <- current[,grep("SampleID|Sample.ID|SAMPLE_ID|Time|file",colnames(current), invert = T, ignore.case = T)]
    colnames(current)
    for (j in 1:ncol(current)){
      current[,j] <- as.numeric(as.character(current[,j]))
      expr <- c(expr,median(current[current[,j] > 0,j]))
    }
    plot_median <- rbind(plot_median, expr)
    expr <- NULL
  }

  plot_median <- data.frame(t(plot_median))
  row.names(plot_median) <- colnames(data[,grep("SAMPLE.*ID", colnames(data), ignore.case = T, invert = T)])
  colnames(plot_median) <- files

  mds <- plotMDS(plot_median, plot = FALSE)
  pca_out <- prcomp(t(plot_median), center = TRUE, scale. = FALSE)
  ggdf <- data.frame(SAMPLE_ID = colnames(plot_median), MDS1 = mds$x, MDS2 = mds$y, PC1 = pca_out$x[,1], PC2 = pca_out$x[,2])
  ggdf$GROUP <- pheno_data[match(ggdf$SAMPLE_ID,pheno_data$SAMPLE_ID),"GROUP"]
  ggdf$CELL_COUNT <- pheno_data[match(ggdf$SAMPLE_ID, pheno_data$SAMPLE_ID),"Filtered_Cell_Count"]

  p5plots <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = GROUP, size = CELL_COUNT, label = SAMPLE_ID)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(range = c(6, 12))+
    geom_label_repel(show.legend = F) +
    # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
    theme_bw() +
    ggtitle("Multidimensional Scaling Plot") +
    scale_color_manual(values = group_colors)
  p5plots <- adjust_theme(p5plots)

  somePDFPath <- paste(cdir,"5URSA_PLOT_FLOW_SAMPLE_MDS_ASINH_MEDIAN_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p5plots)
  dev.off()

  p6plots <- ggplot(ggdf, aes(x = PC1, y = PC2, color = GROUP, size = CELL_COUNT, label = SAMPLE_ID)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(range = c(6, 12))+
    geom_label_repel(show.legend = F) +
    # geom_text(aes(label=SAMPLE_ID), size = 2.5, position = position_jitter(width = 0.05, height = 0))+
    theme_bw() +
    ggtitle("Principle Component Analysis Plot") +
    scale_color_manual(values = group_colors)
  p6plots <- adjust_theme(p6plots)

  somePDFPath <- paste(cdir,"6URSA_PLOT_FLOW_SAMPLE_PCA_ASINH_MEDIAN_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p6plots)
  dev.off()

  p7data <- plot_median
  p7data <- (scale((p7data)))
  p7data <- t(scale(t(p7data)))

  somePDFPath <- paste(cdir,"7URSA_PLOT_FLOW_HEATMAP_SAMPLE_SCALED_ASINH_MEDIAN_EXPRESSION_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7,pointsize=12)
  print(complex_heatmap(p7data, col = color_conditions$BlueYellowRed, legendtitle = "ZScore", col_title = project_name))
  dev.off()

  p8data <- plot_median
  colnames(p8data) <- paste(colnames(p8data), pheno_data[match(colnames(p8data), pheno_data$SAMPLE_ID),"GROUP"])
  p8data <- as.dendrogram(hclust(as.dist(1-cor((p8data)))))

  somePDFPath <- paste(cdir,"8URSA_PLOT_FLOW_DENDROGRAM_SAMPLES_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7,pointsize=12)
  par(mar=c(3,4,1,6))
  print(plot(p8data, horiz = TRUE))
  dev.off()

  NRS <- function(x, ncomp = 3){
    pr <- prcomp(x, center = TRUE, scale. = FALSE)
    score <- rowSums(outer(rep(1, ncol(x)),pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
    return(score)
  }

  nrs_sample <- NULL
  current <- data[,grep("SAMPLE.*ID",colnames(data), ignore.case = T, invert = T)]
  for(i in 1:length(unique(pheno_data$SAMPLE_ID))){
    nrs_sample <- rbind(nrs_sample, NRS(current[which(data$SAMPLE_ID == unique(pheno_data$SAMPLE_ID)[i]),]))
  }

  rownames(nrs_sample) <- unique(pheno_data$SAMPLE_ID)
  nrs_sample <- data.frame(nrs_sample)
  nrs <- colMeans(nrs_sample, na.rm = TRUE)
  markers_ord <- names(sort(nrs, decreasing = TRUE))
  nrs_sample$SAMPLE_ID <- rownames(nrs_sample)
  ggdf <- melt(nrs_sample, id.var = "SAMPLE_ID",
               value.name = "nrs", variable.name = "Markers")
  colnames(ggdf) <- c("SAMPLE_ID","Markers","NRScore")
  ggdf$Markers <- factor(ggdf$Markers, levels = markers_ord)
  ggdf <- merge(ggdf, pheno_data, by.x = "SAMPLE_ID", by.y = "SAMPLE_ID")

  p9plots <- ggplot(ggdf, aes(x = Markers, y = NRScore)) +
    # geom_point(aes(color = SAMPLE_ID), alpha = 0.8,
    # position = position_jitter(width = 0.3, height = 0)) +
    geom_boxplot(aes(fill = GROUP), alpha = 0.8, outlier.color = NA) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = sample_colors)+
    scale_fill_manual(values = group_colors)
  p9plots <- adjust_theme(p9plots, xangle = 45, hejust = 1, vejust = 1, xsize = 15)

  somePDFPath <- paste(cdir,"9URSA_PLOT_FLOW_BOXPLOT_NRS_SCORE_ASINH_MARKER_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=6.7,pointsize=12)
  print(p9plots)
  dev.off()

  samples <- unique(data$SAMPLE_ID)
  data0 <- NULL

  for(i in 1:length(samples)){
    data0[[i]] <- data[which(data$SAMPLE_ID == samples[i]),grep("SAMPLE_ID|SAMPLEID|SAMPLE.ID",colnames(data), ignore.case = T, invert = T)]
  }

  for (i in 1:length(samples)){
    meta <- data.frame(name=colnames(data0[[i]]),desc=colnames(data0[[i]]))
    meta$range <- apply(apply(data0[[i]],2,range),2,diff)
    meta$minRange <- apply(data0[[i]],2,min)
    meta$maxRange <- apply(data0[[i]],2,max)
    data0[[i]] <- new("flowFrame",exprs=as.matrix(data0[[i]]),parameters=AnnotatedDataFrame(meta))
  }


  print("Running dimension reduction and clustering..")
  fs_data = as(data0,"flowSet")
  flowCore::pData(fs_data)$name <- samples
  som_input <- ReadInput(fs_data)
  set.seed(59)
  som <- BuildSOM(som_input)
  codes <- som$map$codes
  nmc <- 90
  print("Running ConsensusClusterPlus..")
  somePDFPath <- paste(cdir,"10URSA_PLOT_FLOW_CONSENSUS_",project_name,".pdf", sep = "")
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
  ck2 <- min(PAC_turn[PAC_turn>1]) + 1
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

  p10plots <- ggplot(PAC, aes(x= K, y= PAC, label = Label)) + geom_line(colour = "grey") +
    ggtitle("Chosen K Number of Clusters Based on PAC Method")+
    theme_classic()+ geom_text_repel(
      max.overlaps = Inf,force=1,
      point.padding = 0, # additional padding around each point
      min.segment.length = 0, # draw all line segments
      max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
      box.padding = 0.3, size = 4, colour = "red")+
    geom_point(colour = ifelse(PAC$K == chosen_k, "red", "grey"), size = ifelse(PAC$K == chosen_k, 5, 2))

  somePDFPath <- paste(cdir,"10URSA_PLOT_FLOW_PAC_CHOSEN_",chosen_k,"_BASED_ON_CONSENSUS_CLUSTERS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p10plots)
  dev.off()

  plot_median <- NULL
  expr <- NULL
  cell_number <- NULL
  pop_list <- 1:chosen_k

  for(i in 1:length(pop_list)){
    current <- data[which(data$population == pop_list[i]),grep("SAMPLE.*ID|population",colnames(data), ignore.case = T, invert = T)]
    cell_number <- c(cell_number,nrow(current))
    for (j in 1:ncol(current)){
      expr <- c(expr,median(current[current[,j] > 0,j]))
    }
    plot_median <- rbind(plot_median, expr)
    expr <- NULL
  }

  row.names(plot_median) <- paste("Cluster_", pop_list, ":",cell_number,sep = "")
  colnames(plot_median) <- toupper(gsub("^.*?_","",colnames(data)[grep("SAMPLE.*ID|population",colnames(data), ignore.case = T, invert = T)]))
  plot_median <- data.frame(plot_median)
  dista <- hclust(as.dist(1-cor(t(plot_median))), method = "complete")
  plot_median <- scale(plot_median)
  plot_median <- t(scale(t(plot_median)))
  plot_median <- as.matrix(plot_median)

  data_meta <- data[,grep("SAMPLE.*ID|population", colnames(data), ignore.case = T)]
  x <- data.frame(t(data[,grep("SAMPLE.*ID|population", colnames(data), ignore.case = T, invert = T)]))
  x <- CreateSeuratObject(counts = x)
  x$PROJECT <- project_name
  x@meta.data <- cbind(x@meta.data,data_meta)
  x$orig.ident <- x$SAMPLE_ID
  x@assays$RNA@data <- x@assays$RNA@counts
  x <- ScaleData(x)
  x <- FindVariableFeatures(x)
  x <- RunPCA(x, features = VariableFeatures(x))
  x <- RunUMAP(x, reduction = "pca", dims = 1:ifelse(length(x@reductions$pca) < 30, length(x@reductions$pca), 30))

  plotx <- data.frame(UMAP_1 = x@reductions$umap@cell.embeddings[,"UMAP_1"],
                      UMAP_2 = x@reductions$umap@cell.embeddings[,"UMAP_2"],
                      PC_1 = x@reductions$pca@cell.embeddings[,"PC_1"],
                      PC_2 = x@reductions$pca@cell.embeddings[,"PC_2"],
                      CELLL_ID = row.names(x@meta.data))

  plotx$CLUSTER <- factor(data[match(plotx$CELLL_ID, row.names(data)), "population"], levels = c(unique(sort(as.numeric(as.character(data$population))))))
  plotx$SAMPLE_ID <- data$SAMPLE_ID

  cluster_colors <- gen_colors(color_conditions$tenx,length(unique(plotx$CLUSTER)))
  names(cluster_colors) <- levels(plotx$CLUSTER)

  p11plots <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTER",
                           plot_title = paste("UMAP - CLUSTER: ",project_name, sep = ""),
                           col = cluster_colors, annot = F, legend_position = "right",
                           numeric = T, point_size = 0.2, label_size = 6, legendsize = 15)

  somePDFPath <- paste(cdir,"11URSA_PLOT_FLOW_UMAP_CLUSTERS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p11plots)
  dev.off()

  p12plots <- ggplot(plotx,  aes(x = UMAP_1, y = UMAP_2, color = CLUSTER)) +
    geom_point(size = 0.2, alpha = 0.7) +
    theme_classic() + facet_wrap(~SAMPLE_ID, ncol = 6) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    scale_color_manual(values = cluster_colors) +
    ggtitle(paste("UMAP - CLUSTER (BY SAMPLES): ",project_name, sep = ""))
  p12plots <- adjust_theme(p12plots)

  # somePDFPath <- paste(cdir,"12URSA_PLOT_FLOW_UMAP_CLUSTERS_BY_SAMPLES_",project_name,".pdf", sep = "")
  # pdf(file=somePDFPath, width=12, height=ceiling(length(unique(plotx$SAMPLE_ID))/4)*2.4,pointsize=12)
  somePNGPath = paste(cdir,"12URSA_PLOT_FLOW_UMAP_CLUSTERS_BY_SAMPLES_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 800*ceiling(length(unique(plotx$SAMPLE_ID))/4),units = "px", res = 300)
  print(p12plots)
  dev.off()

  p13plots <- plot_bygroup(plotx, x = "PC_1", y = "PC_2", group = "CLUSTER",
                           plot_title = paste("PCA - CLUSTER: ",project_name, sep = ""),
                           col = cluster_colors, annot = F, legend_position = "right",
                           numeric = T, point_size = 0.2, label_size = 6, legendsize = 15)
  p13plots <- adjust_theme(p13plots)

  # somePDFPath <- paste(cdir,"13URSA_PLOT_FLOW_PCA_CLUSTERS_",project_name,".pdf", sep = "")
  # pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  somePNGPath = paste(cdir,"13URSA_PLOT_FLOW_PCA_CLUSTERS_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 3000,units = "px", res = 300)
  print(p13plots)
  dev.off()

  p14data <- plot_median
  somePDFPath <- paste(cdir,"14URSA_PLOT_FLOW_HEATMAP_CLUSTER_SCALED_ASINH_MEDIAN_EXPRESSION_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=8,pointsize=12)
  print(complex_heatmap(p14data, col = color_conditions$BlueYellowRed, legendtitle = "ZScore", col_title = project_name))
  dev.off()

  marker_list <- unique(colnames(data[,grep("SAMPLE.*ID|population", colnames(data), ignore.case = T, invert = T)]))
  melt_data <- melt(data.frame(data,UMAP_1 = plotx$UMAP_1, UMAP_2 = plotx$UMAP_2),
                    id.vars = c("UMAP_1","UMAP_2", "SAMPLE_ID","population"))
  colnames(melt_data) <- c("UMAP_1","UMAP_2","SAMPLE_ID","CLUSTER","Marker","Asinh_Expression")

  p15plots <- ggplot(melt_data,  aes(x = UMAP_1, y = UMAP_2, color = Asinh_Expression)) +
    facet_wrap(~Marker) +
    geom_point(size = 0.2) + theme_bw()+
    scale_color_gradientn(colors = gen_colors(c("blue","cyan","green","yellow","orange","red","red4"), 100))
  p15plots <- adjust_theme(p15plots)

  # somePDFPath <- paste(cdir,"15URSA_PLOT_FLOW_UMAP_MARKER_EXPRESSIONS_",project_name,".pdf", sep = "")
  # pdf(file=somePDFPath, width=10, height=3*ceiling(length(unique(melt_data$Marker))/3),pointsize=12)
  somePNGPath = paste(cdir,"15URSA_PLOT_FLOW_UMAP_MARKER_EXPRESSIONS_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height = 1200*ceiling(length(unique(melt_data$Marker))/3),units = "px", res = 300)
  print(p15plots)
  dev.off()

  total_count <- data.frame(table(data$population))
  current <- data.frame(table(data[,c("SAMPLE_ID","population")]))
  current <- current[current$Freq != 0,]
  colnames(current) <- c("SAMPLE_ID","CLUSTER","COUNT")
  current$CLUSTER_TOTAL_COUNT <- total_count[match(current$CLUSTER, total_count$Var1),"Freq"]
  current$PROPORTION <- current$COUNT/current$CLUSTER_TOTAL_COUNT

  node_proportion <- current

  p16plots <- ggplot(node_proportion, aes(CLUSTER, PROPORTION, fill = SAMPLE_ID))+
    geom_bar(stat="identity", alpha=0.8)+
    coord_polar()+
    scale_fill_viridis(discrete = T)+
    ggtitle(paste("Frequency of Samples in Each Cluster: ", project_name, sep = ""))+
    theme_classic()
  p16plots <- adjust_theme(p16plots)

  somePDFPath <- paste(cdir,"16URSA_PLOT_FLOW_SAMPLE_PROPORTIONS_IN_CLUSTERS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p16plots)
  dev.off()

  x$CLUSTER <- plotx[match(colnames(x), plotx$CELLL_ID),"CLUSTER"]
  marker_colors <- gen_colors(color_conditions$bright, nrow(x))
  names(marker_colors) <- row.names(x)
  p17plots <- VlnPlot(x, group.by = "CLUSTER", features = row.names(x), stack = T, flip = T) & xlab("CLUSTERS")

  somePDFPath <- paste(cdir,"17URSA_PLOT_FLOW_VIOLIN_MARKERS_IN_CLUSTERS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=7.5,pointsize=12)
  print(p17plots)
  dev.off()

  saveRDS(data, paste(cdir, "18URSA_DATA_FLOW_FILTERED_DATA_COMBINED_",project_name,".RDS", sep = ""))
  print("Completed!")

}
