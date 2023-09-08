############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Megrez: scCNV
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
#' @import adegenet
#' @import ape
#' @import fastcluster
#' @import GenomicRanges
#' @import ggtree
#' @import intervalaverage
#' @import umap
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format single-cell CNV (scCNV) pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification scCNV
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Ursa_scCNV' by default.
#' @param input_dir Directory to all input files. Current working directory by default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format files.
#' @param diploid.min.event_confidence Filter diploid cells with a minimum copy number
#' event confidence threshold. Default to 100. Refer to:
#' https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/algorithms/cnv_calling
#' @param nondiploid.min.event_confidence Filter non-diploid cells with a minimum copy number
#' event confidence threshold. Default to 100. Refer to:
#' https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/algorithms/cnv_calling
#' @param confidence_threshold_quantile CNV event confidence threshold by quantile percentage.
#' Default is above 25% quantile of data event confidence.
#' @param size_threshold_quantile CNV event length threshold by quantile percentage.
#' Default to 25% quantile of data CNV event length.
#' @param cnv_cell_threshold Filtering out CNV events that are present in less than N% of all cells.
#' Values range from 0 to 100. Default to 5.
#' @param min.cnv_event_size Filtering out CNV events with an interval of less than the indicated threshold.
#' Default to the CNV event length threshold (i.e., size_threshold_quantile) percent quantile of
#' the CNV event lengths.
#' @param min.cnv_event_confidence Filtering out CNV events with event confidence less than the indicated threshold.
#' Default to the CNV event confidence threshold (i.e., confidence_threshold_quantile) percent quantile of
#' the CNV event confidences.
#' @param min.mappable.size Filtering out CNV events with an mappable region interval of less than the indicated threshold.
#' Default to the CNV event length threshold (i.e., size_threshold_quantile) percent quantile of
#' the mappable region intervals.
#' @export
#'
scCNVPip <- function(project_name = "Ursa_scCNV",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file,
                     
                     # QC
                     diploid.min.event_confidence = 100,
                     nondiploid.min.event_confidence = 100,
                     
                     confidence_threshold_quantile = 25,
                     size_threshold_quantile = 25,
                     cnv_cell_threshold = 5,

                     # hclust
                     region_type.hclust.method = "complete",
                     
                     # Filtering CNV
                     min.cnv_event_size = NULL,
                     min.cnv_event_confidence = NULL,
                     min.mappable.size = NULL,
                     
                     # find.clusters
                     adegenet.find.clusters.clust = NULL,
                     adegenet.find.clusters.n.pca = 200,
                     adegenet.find.clusters.n.clust = NULL,
                     adegenet.find.clusters.method = "kmeans",
                     adegenet.find.clusters.stat = "BIC",
                     adegenet.find.clusters.choose.n.clust = FALSE,
                     adegenet.find.clusters.criterion = "goodfit",
                     adegenet.find.clusters.max.n.clust = 10,
                     adegenet.find.clusters.n.iter = 1e5,
                     adegenet.find.clusters.n.start = 10,
                     adegenet.find.clusters.center = TRUE,
                     adegenet.find.clusters.scale = TRUE,
                     adegenet.find.clusters.pca.select = "nbEig",
                     adegenet.find.clusters.perc.pca = NULL,
                     adegenet.find.clusters.dudi = NULL,
                     
                     # dapc
                     adegenet.dapc.n.pca = 200,
                     adegenet.dapc.n.da = 10,
                     adegenet.dapc.center = TRUE,
                     adegenet.dapc.scale = FALSE,
                     adegenet.dapc.var.contrib = TRUE,
                     adegenet.dapc.var.loadings = FALSE,
                     adegenet.dapc.pca.info = TRUE,
                     adegenet.dapc.pca.select = "nbEig",
                     adegenet.dapc.perc.pca = NULL,
                     adegenet.dapc.dudi = NULL,
                     
                     # umap
                     # umap.config = umap.defaults,
                     umap.method = "naive",
                     umap.preserve.seed = TRUE,

                     # hclust
                     cell_binary_cnv.hclust.method = "complete",
                     
                     # hclust
                     ploidy.hclust.method = 'complete',
                     
                     # fastme.bal
                     ape.fastme.bal.nni = TRUE,
                     ape.fastme.bal.spr = TRUE,
                     ape.fastme.bal.tbr = FALSE){
  
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "scCNV", isDir = T)
  ctime <- time_ini()
  ploidy_levels <- c("0","1","2","3","4","5","6","7","8","9",">=10")
  cell_ploidy <- c("diploid", "non-diploid", "noisy")
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir, recursive = T)
  
  sample_files <- list.files(input_dir, recursive = T, full.names = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep(gsub(".*\\/(.*)","\\1",pheno_file, ignore.case = T),sample_files, ignore.case = T, invert = T)]

  data <- NULL
  data_current <- NULL
  data_summary <- NULL
  data_cell_stats <- NULL
  data_mappable <- NULL
  data_group <- NULL
  hc_summary <- NULL

  for(i in 1:nrow(pheno_data)){
    print(pheno_data[i,"SAMPLE_ID"])
    print(paste("Processing sample: ", pheno_data[i,"SAMPLE_ID"], "..",sep = ""))
    current_dir <- list.files(path = input_dir, pattern = pheno_data[i,grep("FILE", colnames(pheno_data), ignore.case = T)], ignore.case = T, full.names = T, recursive = T)
    if(length(current_dir) > 0){
      current_bed <- current_dir[grep(current_dir, pattern = "\\cnv.*calls.*bed$", ignore.case = T)][1]

      print("Reading in cnv calls..")
      current <- read.table(current_bed)
      colnames(current) <- c("chrom","start","end","id","copy_number","event_confidence")
      if(length(grep("chr", unique(current$chrom), ignore.case = T)) > 0){
        print("Standardizing chromosome naming format..")
        current$chrom <- gsub("chr","",current$chrom, ignore.case = T)
      }
      current_summary <- current_dir[grep(current_dir, pattern = "summary_metrics.*\\.csv$", ignore.case = T)][1]
      current_summary <- read.csv(current_summary)
      current_mappable <- current_dir[grep(current_dir, pattern = "\\mappable_regions.*bed$", ignore.case = T)][1]
      print("Reading in mappable regions..")
      current_mappable <- read.table(current_mappable)
      colnames(current_mappable) <- c("chrom","start","end")
      if(length(grep("chr", unique(current_mappable$chrom), ignore.case = T)) > 0){
        print("Standardizing chromosome naming format..")
        current_mappable$chrom <- gsub("chr","",current_mappable$chrom, ignore.case = T)
      }
      data_mappable <- rbind(data_mappable, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current_mappable))
      current_group <- current[which(!((current$id %in% current_summary$cell_id) & (current$id <= max(current_summary$cell_id)))),]
      data_group <- rbind(data_group, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current_group))
      current$Size <- as.numeric(as.character(current$end)) - as.numeric(as.character(current$start))
      current <- split(current, current$chrom)
      current <- lapply(current, function(x){x <- x[order(x$start, decreasing = F),]})
      current <- do.call(rbind.data.frame,current)
      current <- current[which(current$id %in% current_summary$cell_id),]
      current$Ploidy <- current$copy_number
      current[which(current$Ploidy >=10),"Ploidy"] <- ">=10"
      current$Ploidy <- factor(current$Ploidy, levels = ploidy_levels)
      current$Region_type <- ifelse(current$copy_number == 2, "diploid","non-diploid")
      current[current_summary[match(current$id,current_summary$cell_id),"is_noisy"] == 1,"Region_type"] <- "noisy"
      chr_num <- sort(as.numeric(as.character(unique(current$chrom)[grep("X|Y",unique(current$chrom), ignore.case = T, invert = T)])))
      chr_num <- c(as.character(chr_num), unique(current$chrom)[grep("X|Y",unique(current$chrom), ignore.case = T)])

      temp <- split(current, current$chrom)
      temp <- lapply(temp, function(x){
        x <- data.table(x)
        colnames(x)[grep("start", colnames(x), ignore.case = T)] <- "start_old"
        colnames(x)[grep("end", colnames(x), ignore.case = T)] <- "end_old"
        x <- isolateoverlaps(x,interval_vars = c("start_old","end_old"),group_vars = "id")
      })

      paste("Finnished isolating overlaps ..")
      current <- do.call(rbind.data.frame,temp)
      current$range <- paste(current$chrom,":",current$start,"-",current$end, sep = "")
      duplicated_data <- NULL
      duplicated_data <- current[!(!(duplicated(current[,c("id","range")]) | duplicated(current[,c("id","range")], fromLast = TRUE))),]
      if(nrow(duplicated_data) > 0){
        duplicated_data <- split(duplicated_data, duplicated_data$id)
        duplicated_data <- lapply(duplicated_data, function(x){x <- split(x, x$range)})
        duplicated_data <- lapply(duplicated_data, function(x){lapply(x, function(y){
          y <- y[which.max(y$event_confidence),]
        })})
        duplicated_data <- lapply(duplicated_data, function(x){do.call(rbind.data.frame,x)})
        duplicated_data <- do.call(rbind.data.frame,duplicated_data)
        current <- current[(!(duplicated(current[,c("id","range")]) | duplicated(current[,c("id","range")], fromLast = TRUE))),]
        current <- rbind(current,duplicated_data)
        rm(duplicated_data)
      }

      rm(temp)

      current$Region_type_numeric <- current$Region_type
      current$Region_type_numeric <- gsub("non-diploid","2",current$Region_type_numeric)
      current$Region_type_numeric <- gsub("^diploid$","1",current$Region_type_numeric)
      current$Region_type_numeric <- gsub("noisy","-1",current$Region_type_numeric)

      cell_stats <- as.data.frame.matrix(table(data.frame(current[,c("id","Region_type")])))
      cell_stats$Cell_Ploidy <- colnames(cell_stats)[apply(cell_stats,1,which.max)]
      cell_stats$id <- row.names(cell_stats)
      current$Cell_Ploidy <- cell_stats[match(current$id, cell_stats$id),"Cell_Ploidy"]
      current$CID <- paste(current$id, current$range, sep = "-")
      temp <- current[which(tolower(current$Cell_Ploidy) == "diploid" & current$event_confidence > diploid.min.event_confidence),]
      temp <- rbind(temp, current[which(tolower(current$Cell_Ploidy) != "diploid" & current$event_confidence > nondiploid.min.event_confidence),])
      temp <- rbind(temp, current[which(current$range %in% temp$range & !current$CID %in% temp$CID),])
      temp <- reshape2::dcast(temp, formula = range~id, value.var = "Region_type_numeric")
      temp[is.na(temp)] <- 0
      row.names(temp) <- temp$range
      temp <- temp[,grep("range", colnames(temp), ignore.case = T, invert = T)]
      start_time <- Sys.time()
      temp <- as.matrix(temp)
      print("Running clustering..")
      hc <- hclust(dist(t(temp)), method = region_type.hclust.method)
      end_time <- Sys.time()
      print(paste("Clustering took ", round(end_time - start_time, digits = 3), " minutes..",sep = ""))
      hc_summary[[i]] <- hc
      h1names <- colnames(temp)[hc$order]
      current$id <- factor(current$id, levels = colnames(temp)[hc$order])
      cell_stats$Cell_Ploidy <- factor(cell_stats$Cell_Ploidy, levels = c("diploid","non-diploid","noisy"))
      cell_stats <- cell_stats[match(as.character(h1names),cell_stats$id),]
      cell_stats$id <- factor(cell_stats$id, levels = as.character(h1names))
      cell_stats$Mean_ploidy <- round(current_summary[match(cell_stats$id,current_summary$cell_id),"mean_ploidy"])
      current$chrom <- factor(current$chrom, levels = chr_num)
      current$Cell_type <- cell_stats[match(current$id, cell_stats$id),"Cell_Ploidy"]
      current$Cell_type <- factor(current$Cell_type, levels = c("diploid","non-diploid","noisy"))
      current$range <- factor(current$range, levels = unique(current$range))
      data_current[[i]] <- current
      names(data_current)[i] <- pheno_data[i,"SAMPLE_ID"]
      data <- rbind(data, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current))
      data_summary <- rbind(data_summary, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], current_summary))
      data_cell_stats <- rbind(data_cell_stats, data.frame(Sample = pheno_data[i,"SAMPLE_ID"], cell_stats))
      current <- NULL
      cell_stats <- NULL
      current_summary <- NULL
      current_mappable <- NULL
      current_group <- NULL
      print(paste("Done processing ", pheno_data[i,"SAMPLE_ID"]," ..", sep = ""))
    }
  }

  print(paste("Filtering CNV..", sep = ""))
  data_group$Size <- data_group$end - data_group$start
  data_mappable$Size <- data_mappable$end - data_mappable$start
  selected_data <- data[grep("non-diploid", data$Cell_type, ignore.case = T),]
  selected_data$cell_id <- paste(selected_data$Sample, selected_data$id, sep = "_")
  
  if(is.null(min.cnv_event_size)){
    min.cnv_event_size <- quantile(data_group$Size, size_threshold_quantile/100)
  }
  
  if(is.null(min.cnv_event_confidence)){
    min.cnv_event_confidence <- quantile(data_group$event_confidence, confidence_threshold_quantile/100)
  }
  
  if(is.null(min.mappable.size)){
    min.mappable.size <- quantile(data_mappable$Size, size_threshold_quantile/100)
  }
  
  selected_data_group <- data_group[which(data_group$Size >= min.cnv_event_size & data_group$event_confidence > min.cnv_event_confidence),]
  selected_data_group <- selected_data_group[which(!selected_data_group$copy_number == 2),]
  selected_data_mappable <- data_mappable[which(data_mappable$Size >= min.mappable.size),]
  final_data_group <- NULL
  for(i in 1:nrow(selected_data_mappable)){
    temp <- selected_data_group[which((selected_data_group$Sample == selected_data_mappable[i,"Sample"]) &
                                        (selected_data_group$chrom == selected_data_mappable[i,"chrom"]) &
                                        (selected_data_group$start >= selected_data_mappable[i,"start"]) &
                                        (selected_data_group$end <= selected_data_mappable[i,"end"])),]
    if(nrow(temp) > 0){
      # print(i)
      final_data_group <- rbind(final_data_group, temp)
    }
    temp <- NULL
  }

  total_cells <- length(unique(selected_data$cell_id))
  cnv_events <- unique(final_data_group[,grep("chrom|start|end", colnames(final_data_group), ignore.case = T)])
  final_cnv_events <- NULL
  for(i in 1:nrow(cnv_events)){
    temp <- selected_data[which((selected_data$chrom == cnv_events[i,"chrom"]) &
                                  (selected_data$start_old >= cnv_events[i,"start"]) &
                                  (selected_data$end_old <= cnv_events[i,"end"])),]
    if(nrow(temp) > 0){
      temp <- temp[!duplicated(temp[,grep("cell_id|chrom|start_old|end_old", colnames(temp), ignore.case = T)]),]
      # print(length(unique(temp$cell_id))/total_cells)
      if(length(unique(temp$cell_id))/total_cells > cnv_cell_threshold/100){
        final_cnv_events <- rbind(final_cnv_events, cnv_events[i,])
      }
    }
    temp <- NULL
  }

  final_cnvs <- NULL
  for(i in 1:nrow(final_cnv_events)){
    temp <- selected_data_mappable[which((selected_data_mappable$chrom == final_cnv_events[i,"chrom"]) &
                                           (selected_data_mappable$start <= final_cnv_events[i,"start"]) &
                                           (selected_data_mappable$end >= final_cnv_events[i,"end"])),]
    if(nrow(temp) > 0){
      final_cnvs <- rbind(final_cnvs, temp)
    }
    temp <- NULL
  }

  final_cnvs <- final_cnvs[!duplicated(final_cnvs[,grep("chrom|start|end", colnames(final_cnvs), ignore.case = T)]),grep("chrom|start|end|Size|range", colnames(final_cnvs), ignore.case = T)]
  final_cnvs$range <- paste(final_cnvs$chrom,":",final_cnvs$start,"-",final_cnvs$end, sep = "")

  melt_cell_binary_cnv <- NULL
  selected_cells <- unique(selected_data$cell_id)
  for(i in 1:nrow(final_cnvs)){
    temp <- selected_data[which((selected_data$chrom == final_cnvs[i,"chrom"]) &
                                  (selected_data$start_old >= final_cnvs[i,"start"]) &
                                  (selected_data$end_old <= final_cnvs[i,"end"])),]

    if(nrow(temp) > 0){
      temp <- unique(temp$cell_id)
      melt_cell_binary_cnv <- rbind(melt_cell_binary_cnv, data.frame(CNV_Event = final_cnvs[i,"range"], cell_id = temp , Count = 1))
    }
    temp <- NULL
  }

  cell_binary_cnv <- reshape2::dcast(melt_cell_binary_cnv, cell_id ~ CNV_Event, value.var = "Count")
  cell_binary_cnv[is.na(cell_binary_cnv)] <- 0
  row.names(cell_binary_cnv) <- cell_binary_cnv$cell_id
  cell_binary_cnv <- cell_binary_cnv[,grep("cell_id", colnames(cell_binary_cnv), ignore.case = T, invert = T)]

  print(paste("Running dimension reduction and clustering.. ", sep = ""))
  cell_groups <- find.clusters(cell_binary_cnv,
                               clust = adegenet.find.clusters.clust,
                               n.pca = adegenet.find.clusters.n.pca,
                               n.clust = adegenet.find.clusters.n.clust,
                               method = adegenet.find.clusters.method,
                               stat = adegenet.find.clusters.stat,
                               choose.n.clust = adegenet.find.clusters.choose.n.clust,
                               criterion = adegenet.find.clusters.criterion,
                               max.n.clust = adegenet.find.clusters.max.n.clust,
                               n.iter = adegenet.find.clusters.n.iter,
                               n.start = adegenet.find.clusters.n.start,
                               center = adegenet.find.clusters.center,
                               scale = adegenet.find.clusters.scale,
                               pca.select = adegenet.find.clusters.pca.select,
                               perc.pca = adegenet.find.clusters.perc.pca,
                               dudi = adegenet.find.clusters.dudi)
  dapc_out <- dapc(cell_binary_cnv,
                   grp = cell_groups$grp,
                   n.pca = adegenet.dapc.n.pca,
                   n.da = adegenet.dapc.n.da,
                   center = adegenet.dapc.center,
                   scale = adegenet.dapc.scale,
                   var.contrib = adegenet.dapc.var.contrib,
                   var.loadings = adegenet.dapc.var.loadings,
                   pca.info = adegenet.dapc.pca.info,
                   pca.select = adegenet.dapc.pca.select,
                   perc.pca = adegenet.dapc.perc.pca,
                   dudi = adegenet.dapc.dudi)
  # print(head(dapc_out))
  dc_cnvs <- dapc_out$var.contr
  umap_out <- umap::umap(cell_binary_cnv, 
                         # config = umap.config,
                         method = umap.method,
                         preserve.seed = umap.preserve.seed)
  umap_coords <- data.frame(UMAP_1 = umap_out$layout[,1], UMAP_2 = umap_out$layout[,2],
                            cell_id = row.names(umap_out$layout))
  umap_coords$Cluster <- dapc_out$assign[match(umap_coords$cell_id,row.names(dapc_out$posterior))]
  umap_coords <- umap_coords[order(umap_coords$Cluster, decreasing = F),]
  umap_coords$cell_id <- factor(umap_coords$cell_id, levels = unique(umap_coords$cell_id))

  hc_cnv <- hclust(dist(t(cell_binary_cnv)), method=cell_binary_cnv.hclust.method)
  binary_cnv_events <- cell_binary_cnv
  binary_cnv_events$cell_id <- row.names(cell_binary_cnv)
  binary_cnv_events <- reshape2::melt(binary_cnv_events)
  colnames(binary_cnv_events) <- c("cell_id","CNV_Event","Count")
  binary_cnv_events$Count <- factor(binary_cnv_events$Count)
  binary_cnv_events$CNV_Event <- factor(binary_cnv_events$CNV_Event, levels = colnames(cell_binary_cnv)[hc_cnv$order])
  binary_cnv_events$Cluster <- dapc_out$assign[match(binary_cnv_events$cell_id,row.names(dapc_out$posterior))]
  binary_cnv_events <- binary_cnv_events[order(binary_cnv_events$Cluster, decreasing = F),]
  binary_cnv_events$cell_id <- factor(binary_cnv_events$cell_id, levels = unique(binary_cnv_events$cell_id))

  top_events <- data.frame(dapc_out$var.contr)
  top_events <- top_events[order(top_events$LD1, decreasing = T),]
  ld1 <- row.names(top_events)[1:50]
  top_events <- top_events[order(top_events$LD2, decreasing = T),]
  ld2 <- row.names(top_events)[1:50]
  top_events <- c(unique(ld1), unique(ld2))[1:50]

  top_event_chroms <- data.frame(table(gsub("(.*?):.*","\\1",top_events)))
  colnames(top_event_chroms) <- c("chrom","Count")
  top_event_chrom_threshold <- 3
  selected_chroms <- top_event_chroms[which(top_event_chroms$Count >= top_event_chrom_threshold),"chrom"]
  select_bychrom_data <- selected_data[which(selected_data$chrom %in% selected_chroms),]
  select_bychrom_data$Cluster <- dapc_out$assign[match(select_bychrom_data$cell_id,row.names(dapc_out$posterior))]
  select_bychrom_data <- select_bychrom_data[order(select_bychrom_data$Cluster, decreasing = F),]
  select_bychrom_data$CID <- paste(select_bychrom_data$cell_id, select_bychrom_data$range, sep = "-")
  temp <- split(select_bychrom_data, select_bychrom_data$Cluster)
  for(i in 1:length(temp)){
    x <- temp[[i]]
    if(nrow(x) > 1){
      x$Ploidy <- gsub(">=10","10",x$Ploidy)
      y <- x[which(x$event_confidence > quantile(x$event_confidence, confidence_threshold_quantile/100)),]
      y <- rbind(y, x[which(x$range %in% y$range & !x$CID %in% y$CID),])
      if(nrow(y) > 100){
        x <- y
      }
      x <- dcast(x, range ~ cell_id, value.var = "Ploidy")
      row.names(x) <- x$range
      x <- x[,which(colnames(x) != "range")]
      x[is.na(x)] <- 0
      hc <- fastcluster::hclust(dist(t(x)), method=ploidy.hclust.method)
      x <- data.frame(Cluster = unique(temp[[i]]$Cluster), cell_id = colnames(x)[hc$order])
      temp[[i]] <- x
    }
  }

  temp <- do.call(rbind.data.frame, temp)
  datahcnames <- unique(temp$cell_id)
  select_bychrom_data <- select_bychrom_data[which(select_bychrom_data$cell_id %in% datahcnames),]
  select_bychrom_data$cell_id <- factor(select_bychrom_data$cell_id, levels = datahcnames)
  select_bychrom_data <- select_bychrom_data[,grep("cell_id|^range$|Ploidy|chrom|Cluster", colnames(select_bychrom_data), ignore.case = T)]
  select_bychrom_data <- select_bychrom_data[!duplicated(select_bychrom_data),]
  select_bychrom_data <- split(select_bychrom_data, select_bychrom_data$chrom)
  select_bychrom_data <- lapply(select_bychrom_data, function(x){
    x <- x[order(as.numeric(as.character(gsub(".*?:(.*)-(.*)","\\2",x$range)))),]
  })
  select_bychrom_data <- do.call(rbind.data.frame, select_bychrom_data)
  select_bychrom_data$range <- factor(select_bychrom_data$range, levels = unique(select_bychrom_data$range))
  select_bychrom_data_bar <- data.frame(cell_id = select_bychrom_data$cell_id, Cluster = select_bychrom_data$Cluster)
  select_bychrom_data_bar <- select_bychrom_data_bar[match(datahcnames, select_bychrom_data_bar$cell_id),]
  select_bychrom_data_bar <- unique(select_bychrom_data_bar)

  data_prop <- data.frame(cell_id = unique(selected_data$cell_id), Type = "non-diploid")
  data_prop$Cluster <- dapc_out$assign[match(data_prop$cell_id,row.names(dapc_out$posterior))]
  data_prop$Sample <- gsub("(.*)_.*","\\1",data_prop$cell_id)
  data_prop <- as.data.frame.matrix(table(data_prop[,c("Sample","Cluster")]))
  current <- data[grep("noisy", data$Cell_type, ignore.case = T, invert = T),c("Sample","id")]
  current <- current[!duplicated(current),]
  total_cells <- data.frame(table(current[,c("Sample")]))
  total_cells <- total_cells[match(row.names(data_prop), total_cells$Var1),]
  for(i in 1:nrow(data_prop)){
    data_prop[i,] <- data_prop[i,]/total_cells[i,"Freq"]
  }
  data_prop$Sample <- row.names(data_prop)
  data_prop <- melt(data_prop)
  colnames(data_prop) <- c("Sample","Cluster","Proportion")

  cn_clusters <- data.frame(cell_id = unique(selected_data$cell_id))
  cn_clusters$Cluster <- dapc_out$assign[match(cn_clusters$cell_id,row.names(dapc_out$posterior))]
  cn_clusters$copy_number <- round(data_summary$mean_ploidy)[match(cn_clusters$cell_id, paste(data_summary$Sample,data_summary$cell_id,sep = "_"))]
  cn_clusters <- cn_clusters[which(!is.na(cn_clusters$Cluster) & !is.na(cn_clusters$copy_number)),]
  dcast_cn_clusters <- reshape2::dcast(cn_clusters, cell_id ~ Cluster, value.var = "copy_number")
  dcast_cn_clusters[is.na(dcast_cn_clusters)] <- 0
  row.names(dcast_cn_clusters) <- dcast_cn_clusters$cell_id
  dcast_cn_clusters <- dcast_cn_clusters[,grep("cell_id", colnames(dcast_cn_clusters), ignore.case = T, invert = T)]
  tree_est <- fastme.bal(dist(t(dcast_cn_clusters)),
                         nni = ape.fastme.bal.nni,
                         spr = ape.fastme.bal.spr,
                         tbr = ape.fastme.bal.tbr)

  plotx <- split(cn_clusters, cn_clusters$Cluster)
  plotx <- lapply(plotx, function(x){
    x <- data.frame(Cluster = unique(x$Cluster),
                    Median_CN = median(x$copy_number),
                    Size = nrow(x))
  })
  plotx <- do.call(rbind.data.frame,plotx)
  plotx$Cluster <- as.character(plotx$Cluster)
  plotx <- plotx[match(tree_est$tip.label, plotx$Cluster),]
  clusters <- split(plotx[,"Cluster"], plotx$Cluster)
  tree_est <- groupOTU(tree_est, clusters)
  tree_est$plotx <- plotx

 ccolor_schemes <- ini_colorschemes(assay = "scCNV", data = list(data_cell_stats = data_cell_stats, umap_coords = umap_coords))

  results <- NULL
  results$project <- project_name
  results$assay <- "scCNV"
  results$data <- data_current
  results$cell_stats <- data_cell_stats
  results$dim <- umap_coords
  results$binary_cnv <- binary_cnv_events
  results$selected_chrom <- select_bychrom_data
  results$chrom_bar <- select_bychrom_data_bar
  results$proportion <- data_prop
  results$tree <- tree_est
  results$color_schemes <- ccolor_schemes

  plot_ploidy(results, cdir)

  p <- NULL
  p <- ggplot(results$cell_stats, aes(x = Mean_ploidy, y = Sample)) +
    geom_density_ridges(aes(fill = Sample)) + ggtitle(results$project) +
    xlab("Mean single-cell ploidy") + ylab("Samples") + theme_classic() +
    scale_fill_manual(values = results$color_schemes$sample_colors)
  p <- adjust_theme(p, xsize = 20, title_size = 25)

  somePNGPath <- paste(cdir,"2URSA_PLOT_scCNV_PLOT_PLOIDY_INFO_MEAN_CELL_PLOIDY_",results$project,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
  print(p)
  dev.off()

  somePNGPath <- paste(cdir,"3URSA_PLOT_scCNV_DAPC_COMPONENT12_",results$project,".png", sep = "")
  png(somePNGPath, width = 3000, height =2000, units = "px", res = 300)
  print(scatter(dapc_out, bg="white", col = results$color_schemes$cluster_colors, legend=T,
                scree.da=FALSE, inset.solid = 0.6,cex.lab = 1, label.inds = c("DAPC_1","DAPC_2")))
  dev.off()

  p <- NULL
  p <- plot_bygroup(results$dim, x = "UMAP_1", y = "UMAP_2", group = "Cluster", plot_title = results$project,
                    col = results$color_schemes$cluster_colors, annot = F, legend_position = "right", numeric = T,
                    point_size = 0.5, label_size = 8,legendsize = 15)

  somePNGPath <- paste(cdir,"4URSA_PLOT_scCNV_UMAP_COMPONENT12_",results$project,".png", sep = "")
  png(somePNGPath, width = 3000, height =2000, units = "px", res = 300)
  print(p)
  dev.off()

  plot_ploidy_binary(results, output_dir = cdir)
  plot_selected_chroms(results, output_dir = cdir)

  p <- NULL
  p <- ggplot(results$proportion, aes(x=Sample, y=Proportion, fill=Cluster, group = Cluster)) +
    geom_area(alpha=0.9, size=0.5, colour = "black")+
    theme_classic()+ggtitle(results$project)+
    ylab("Fraction of Cells")+
    scale_fill_manual(values = results$color_schemes$cluster_colors)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p <- adjust_theme(p, legend = "right", xsize = 20, xangle = 45, hejust = 1, vejust = 1)

  somePNGPath <- paste(cdir,"7URSA_PLOT_scCNV_CLUSTER_PROPORTION_BY_SAMPLES_",results$project,".png", sep = "")
  png(somePNGPath, width = 4000, height = 3000, units = "px", res = 300)
  print(p)
  dev.off()

  p <- NULL
  p <- ggtree(tree_est, layout="daylight",
              branch.length = "none", # c(tree_est$plotx$Size/300,rep(min(tree_est$plotx$Size/300),(nrow(tree_est$edge)+1) - nrow(tree_est$plotx)))
              aes(color = group,
                  size = c(tree_est$plotx$Size/50,rep(min(tree_est$plotx$Size/100),(nrow(tree_est$edge)+1) - nrow(tree_est$plotx))))) +
    geom_tiplab(hjust = -2, offset=.1) +
    scale_color_manual(values = results$color_schemes$cluster_colors) +
    theme(legend.position="none") +
    scale_size_continuous()

  somePNGPath <- paste(cdir,"8URSA_PLOT_scCNV_PHYLOGENETIC_TREE_CLUSTERS_",results$project,".png", sep = "")
  png(somePNGPath, width = 3000, height =2000, units = "px", res = 300)
  print(p)
  dev.off()

  saveRDS(results, paste(cdir,"9URSA_DATA_scCNV_RESULTS_",results$project,".RDS", sep = ""))
  print("Completed!")

}
