############################################################################################################
# Ursa: an automated multi-omics package for single-cell analysis
# Dubhe: scImmune
# Version: V1.1.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Date: 2022-02-17
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
#' @export
scImmunePip <- function(project_name = "Ursa_scImmune",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file,
                     integration_method = c("harmony","seurat"),
                     auto_move = T){
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
  dir.create(cdir)
  hpca.se <- HumanPrimaryCellAtlasData()
  contig_dir <- paste(cdir, "/Contig/", sep = "")
  dir.create(contig_dir)
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

  contig_monarch <- repLoad(contig_dir,.mode = "paired")
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

  imm_ov1 <- repOverlap(contig_monarch$data, .method = "morisita", .verbose = F, .col = "nt+v+d+j")
  imm_ov2 <- repOverlap(contig_monarch$data, .method = "morisita", .verbose = F, .col = "aa")

  p1 <- vis(imm_ov1, .text.size = 6) + ggtitle("CLONOTYPES SHARED\nMORISITA OVERLAP INDEX (NT+VDJ)")
  p2 <- vis(imm_ov2, .text.size = 6) + ggtitle("CLONOTYPES SHARED\nMORISITA OVERLAP INDEX (AA)")
  p7plots  <- p1 + p2

  # Build public repertoire table using CDR3 nucleotide sequences
  pr.ntv <- pubRep(contig_monarch$data, .col = "nt+v", .quant = "count", .verbose = F)
  pr.aav <- pubRep(contig_monarch$data, .col = "aa", .quant = "count", .verbose = F)
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
  for(j in 1:length(all_genes)){
    current_gu <- geneUsage(contig_monarch$data, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
    current_gu <- current_gu[!is.na(current_gu$Names),]
    current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
    imm_gu_js <- geneUsageAnalysis(current_gu, .method = "js", .verbose = F)
    imm_gu_cor <- geneUsageAnalysis(current_gu, .method = "cor", .verbose = F)
    p1 <- vis(imm_gu_js, .title = paste("Gene usage JS-divergence: ", toupper(all_genes[j]), sep = ""), .leg.title = "JS", .text.size = 4)
    p2 <- vis(imm_gu_cor, .title = paste("Gene usage correlation: ", toupper(all_genes[j]), sep = ""), .leg.title = "Cor", .text.size = 4)
    p <- p1 + p2
    p10plots[[j]] <- p
    names(p10plots)[j] <- all_genes[j]
  }

  p11plots <- NULL
  for(j in 1:length(all_genes)){
    current_gu <- geneUsage(contig_monarch$data, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
    current_gu <- current_gu[!is.na(current_gu$Names),]
    current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
    row.names(current_gu) <- current_gu$Names
    p <- vis(geneUsageAnalysis(current_gu, "cosine+hclust", .verbose = F), .rect = T)+ggtitle(paste("Optimal number of clusters\n(Clustered based on: ", toupper(all_genes[j]), ")", sep = ""))
    p11plots[[j]] <- p
    names(p11plots)[j] <- all_genes[j]
  }

  p12plots <- NULL

  if(length(contig_monarch$data) > 3){

    for(j in 1:length(all_genes)){
      current_gu <- geneUsage(contig_monarch$data, .gene = paste("hs.", all_genes[j], sep = ""), .norm = T)
      current_gu <- current_gu[!is.na(current_gu$Names),]
      current_gu <- current_gu[unlist(lapply(strsplit(current_gu$Names,split = ";"), function(x){length(grep("TRUE",!grepl("NA",x, ignore.case = T))) > 0})),]
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
    p1 <- vis(spectratype(contig_monarch$data[[i]], .quant = "id", .col = "nt"))+
      ggtitle(paste("SAMPLE: ",names(contig_monarch$data)[i],sep = "")) +
      xlab('CCDR3 Nucleotide length')
    p2 <- vis(spectratype(contig_monarch$data[[i]], .quant = "count", .col = "aa+v")) +
      xlab('CCDR3 AA length')
    p <- p1 + p2
    p13plots[[i]] <- p
    names(p13plots)[i] <- names(contig_monarch$data)[i]
  }

  # Diversity estimation
  diversity_estimates <- NULL
  diversity_methods <- c("chao1", "hill", "div", "gini.simp", "inv.simp", "d50","raref")
  p14plots <- NULL

  for(i in 1:length(diversity_methods)){
    p14plots[[i]] <- vis(repDiversity(contig_monarch$data, .method = diversity_methods[i]))
    names(p14plots)[i] <- diversity_methods[i]
  }

  groups <- unique(contig_monarch$meta$CELL_TYPE)

  n <- 10
  p15plots <- NULL
  for(i in 1:length(contig_monarch$data)){
    tc1 <- trackClonotypes(contig_monarch$data, list(i, n), .col = "nt")
    tc2 <- trackClonotypes(contig_monarch$data, list(i, n), .col = "aa")
    p1 <- vis(tc1)+ggtitle(paste(names(contig_monarch$data)[i], ": TOP ",n," CLONOTYPES (CDR3 NT)", sep = ""))+theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=1))
    p2 <- vis(tc2)+ggtitle(paste(names(contig_monarch$data)[i], ": TOP ",n," CLONOTYPES (CDR3 AA)", sep = ""))+theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=1))
    p <- p1/p2
    print(p)
    p15plots[[i]] <- p
    names(p15plots)[i] <- names(contig_monarch$data)[i]
  }

  result_vdjdb <- dbAnnotate(contig_monarch$data, vdjdb, "CDR3.aa", "cdr3")

  # kmer

  p16data <- repOverlap(contig_monarch$data, .col = "v+d+j")

  ov <- scale(p16data)
  ov <- t(scale(t(ov)))
  ov[is.na(ov)] <- 0
  p17plots <- complex_heatmap(ov, col = color_conditions$BlueYellowRed, legendtitle = "Z-Score", col_title = "Repetoire Overlaps (VDJ Regions)")

  print("Processing..")

  contig_types <- c("BCR","TAB","TGD")
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
                             removeNA = TRUE)
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

  tab_presence <- NULL
  tgd_presence <- NULL
  contig_tcr_tab <- NULL
  contig_tcr_tgd <- NULL
  tcr_pheno <- pheno_data[grepl("T_Cell|TCell|T.*Cell|TCR", pheno_data$CELL_TYPE, ignore.case = T) & grepl("contig", pheno_data$DATA_TYPE, ignore.case = T),]
  if(length(tcr_list) > 0){
    contig_tcr <- contig_list[which(names(contig_list) %in% tcr_list)]
    tab_presence <- as.character(unlist(lapply(contig_tcr, function(x){x <- ifelse(length(grep("TRA|TRB",unique(x$chain), ignore.case = T)) > 0, "TRUE")})))
    if(length(grep("TRUE", tab_presence, ignore.case = T)) > 0){
      contig_tcr_tab <- combineTCR(contig_tcr[which(tab_presence == "TRUE")],
                                   samples = tcr_pheno[match(names(contig_tcr[which(tab_presence == "TRUE")]),tcr_pheno$SID),"SID"],
                                   ID = tcr_pheno[match(names(contig_tcr[which(tab_presence == "TRUE")]),tcr_pheno$SID),"PAIR_ID"],
                                   removeNA = TRUE,
                                   cells = "T-AB")
      for(i in 1:length(contig_tcr_tab)){
        contig_tcr_tab[[i]]$SID <- tcr_pheno[match(names(contig_tcr[which(tab_presence == "TRUE")]),tcr_pheno$SID),"SID"][i]
      }
      contig_table$TAB <- contig_tcr_tab
    }else{
      contig_table$TAB <- NULL
    }

    tgd_presence <- as.character(unlist(lapply(contig_tcr, function(x){x <- ifelse(length(grep("TRG|TRD",unique(x$chain), ignore.case = T)) > 1, "TRUE", "FALSE")})))
    if(length(grep("TRUE", tgd_presence, ignore.case = T)) > 0){
      contig_tcr_tgd <- combineTCR(contig_tcr[which(tgd_presence == "TRUE")],
                                   samples = pheno_data[match(names(contig_tcr[which(tgd_presence == "TRUE")]),pheno_data$SAMPLE_ID),"SID"],
                                   # samples = pheno_data[match(names(contig_tcr[which(tgd_presence == "TRUE")]),pheno_data$SAMPLE_ID),"PAIR_ID"],
                                   removeNA = TRUE,
                                   ID = pheno_data[match(names(contig_tcr[which(tgd_presence == "TRUE")]),pheno_data$SAMPLE_ID),"PAIR_ID"],
                                   cells = "T-GD")
      for(i in 1:length(contig_tcr_tgd)){
        contig_tcr_tgd[[i]]$ID <- pheno_data[match(names(contig_tcr[which(tgd_presence == "TRUE")]),pheno_data$SID),"SID"][i]
      }
      contig_table$TGD <- contig_tcr_tgd
    }else{
      contig_table$TGD <- NULL
    }
  }

  # contig_bcr
  # contig_tcr_tab
  # contig_tcr_tgd
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
  name_ref <- data.frame(DATA_TYPE = c("BCR","TAB","TGD"),
                         DESC = c("BCR","ALPHA-BETA TCR","GAMMA-DELTA TCR"))

  for(i in 1:length(contig_table)){

    current_type <- names(contig_table)[i]
    current_name <- name_ref[match(current_type, name_ref$DATA_TYPE),"DESC"]

    if(length(contig_table[[i]]) > 0){
      contig_table[[i]] <- lapply(contig_table[[i]], function(x){
        x <- data.frame(x, contig_pheno[match(unique(x$sample), contig_pheno$SID),])
      })

      number_contig <- quantContig(contig_table[[i]], cloneCall="gene+nt", scale = T, exportTable = T)

      p1 <- abundanceContig(contig_table[[i]], cloneCall = "gene+nt", scale = F) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              plot.title = element_text(size = 20, face = "bold")) +
        guides(fill = guide_legend(ncol = 1))

      p2 <- abundanceContig(contig_table[[i]], cloneCall = "gene+nt", scale = T) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              plot.title = element_text(size = 20, face = "bold")) +
        guides(fill = guide_legend(ncol = 1))

      p18plots[[length(p18plots)+1]] <- (p1|p2)+plot_annotation(title = paste(current_name, " CLONOTYPE ABUNDANCE (VDJC GENE + CDR3 NT)", sep = ""),theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      names(p18plots)[length(p18plots)] <- current_name

      current <- abundanceContig(contig_table[[i]], cloneCall = "aa", exportTable = T)
      p19data[[length(p19data)+1]] <- current[order(current$Abundance, decreasing = T),]
      names(p19data)[length(p19data)] <- current_name

      p1 <- lengthContig(contig_table[[i]], cloneCall="aa", chain = "both")+
        ggtitle("CDR3 AA")+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              axis.text.x = element_text(angle = 0, hjust=1),
              plot.title = element_text(size = 20, face = "bold")) +
        guides(fill = guide_legend(ncol = 1))+
        scale_fill_tableau()
      p1 <- adjust_theme(p1)

      p2 <- lengthContig(contig_table[[i]], cloneCall="nt", chain = "both")+
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
      temp <- compareClonotypes(contig_table[[i]], cloneCall="aa", exportTable = T)
      # temp <- compareClonotypes_debugged(contig_table[[i]], cloneCall = "aa")
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)
      p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="aa", graph = "alluvial")
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]]+
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = ceiling(length(unique(temp$Clonotypes))/12))) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (CDR3 AA)", sep = ""))
      names(p21plots)[k] <- paste("AA_",current_name, sep = "")
      k <- k+1

      temp <- compareClonotypes(contig_table[[i]],cloneCall="gene", exportTable = T)
      # temp <- compareClonotypes_debugged(contig_table[[i]], cloneCall = "gene")
      temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
      temp <- temp[order(temp$Proportion, decreasing = T),]
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)

      p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="gene", graph = "alluvial")
      # plots[[2]] <- plot_alluvial(temp, x = "Sample", y = "Proportion", fill = "Clonotypes",
      #                             group = "Clonotypes", stratum = "Clonotypes",
      #                             alluvium = "Clonotypes", label = "Clonotypes")
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]] +
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = ceiling(length(unique(temp$Clonotypes))/12))) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (VDJC GENE)", sep = ""))
      names(p21plots)[k] <- paste("GENE_",current_name, sep = "")
      k <- k+1

      temp <- compareClonotypes(contig_table[[i]],cloneCall="nt", exportTable = T)
      # temp <- compareClonotypes_debugged(contig_table[[i]], cloneCall = "nt")
      temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
      temp <- temp[order(temp$Proportion, decreasing = T),]
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)

      p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="nt", graph = "alluvial")
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]] +
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (CDR3 NT)", sep = ""))
      names(p21plots)[k] <- paste("NT_",current_name, sep = "")
      k <- k+1

      temp <- compareClonotypes(contig_table[[i]],cloneCall="gene+nt", exportTable = T)
      temp <- temp[!unlist(lapply(strsplit(as.character(temp$Clonotypes), split = "_"), function(x){length(x) == length(grep("NA",x))})),]
      temp <- temp[order(temp$Proportion, decreasing = T),]
      temp <- split(temp, temp$Sample)
      temp <- lapply(temp, function(x){
        x<- x[order(x$Proportion, decreasing = T),]
        x <- x[1:n,]})
      temp <- do.call(rbind.data.frame, temp)

      p21plots[[k]] <- compareClonotypes(contig_table[[i]], clonotypes = unique(as.character(temp$Clonotypes)), cloneCall="gene+nt", graph = "alluvial")
      p21plots[[k]] <- adjust_theme(p21plots[[k]], xangle = 0,legend = "bottom", title_size = 15, legend_title = 12, legend_text = 8)
      p21plots[[k]] <- p21plots[[k]] +
        scale_fill_manual(values = colorRampPalette(color_conditions$general)(length(unique(temp$Clonotypes)))) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust=1),plot.title = element_text(size = 15, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name,"\nTOP CLONOTYPES (VDJC GENE + CDR3 NT)", sep = ""))
      names(p21plots)[k] <- paste("GENE_NT_",current_name, sep = "")
      k <- k+1

      p <- NULL
      p <- vizGenes(contig_table[[i]], gene="V", scale = T)+
        ggtitle(paste(current_name,": V GENE USAGE", sep = ""))+
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              plot.title = element_text(size = 15))
      p22plots[[i]] <- p
      names(p22plots)[i] <- current_name

      p1 <- NULL
      p2 <- NULL
      p1 <- clonalHomeostasis(contig_table[[i]], cloneCall = "gene+nt") +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 20, face = "bold"))+
        scale_x_discrete(labels= names(contig_table[[i]]))
      p2 <- clonalProportion(contig_table[[i]], cloneCall = "gene+nt") +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45, hjust=1))
      p23plots[[i]] <- p1 + p2 + plot_annotation(title = paste(current_name,": VDJC GENE + CDR3 NUCLEOTIDE CLONOTYPE", sep = ""), theme =  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      names(p23plots)[i] <- current_name

      p24plots[[i]] <- clonalOverlap(contig_table[[i]], cloneCall = "gene+nt", method = "morisita") +
        scale_fill_continuous_tableau(palette = "Blue-Teal") +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
              plot.title = element_text(size = 15, face = "bold")) +
        ggtitle(paste(current_name, ": MORISITA INDEX FOR CLONAL OVERLAP (VDJC GENE + CDR3 NT)", sep = ""))
      p24plots[[i]] <- adjust_theme(p24plots[[i]], title_size = 15)
      names(p24plots)[i] <- current_name

      t4 <- clonesizeDistribution(contig_table[[i]], cloneCall = "gene+nt", method="ward.D2", exportTable = T)
      hclust <- hclust(as.dist(t4), method = "ward.D2")
      p25data[[i]] <- as.dendrogram(hclust)
      names(p25data)[i] <- current_name # JENSEN-SHANNON DISTANCE CLUSTERING (VDJC GENE + CDR3 NT)

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "aa", group = "sample")+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (CDR3 AA)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("AA_", current_name, sep = "")

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "gene", group = "sample")+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (VDJC GENE)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("GENE_", current_name, sep = "")

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "nt", group = "sample")+
        theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "right",
              # axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 10, face = "bold")) +
        guides(fill = guide_legend(ncol = 1)) +
        ggtitle(paste(current_name, ": CLONAL DIVERSITY OF SAMPLES (CDR3 NT)", sep = ""))
      names(p26plots)[length(p26plots)] <- paste("NT_", current_name, sep = "")

      p26plots[[length(p26plots)+1]] <- clonalDiversity(contig_table[[i]], cloneCall = "gene+nt", group = "sample")+
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
        data_current[[j]] <- CreateSeuratObject(counts = temp, project = project_name, min.cells = 3, min.features = 200)
      }
      names(data_current)[j] <- current_name
      DefaultAssay(data_current[[j]]) <- "RNA"
      Idents(data_current[[j]]) <- "orig.ident"
      data_current[[j]][["percent.mt"]] <- PercentageFeatureSet(data_current[[j]], pattern = "^MT-")
      data_current[[j]] <- subset(data_current[[j]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
      # data_current[[j]] <- suppressWarnings(SCTransform(data_current[[j]], verbose = FALSE))
      data_current[[j]]@meta.data <- cbind(data_current[[j]]@meta.data, rna_pheno[match(current_name, rna_pheno$SID),])
      # DefaultAssay(data_current[[j]]) <- "RNA"
      data_current[[j]] <- NormalizeData(data_current[[j]])
      data_current[[j]] <- FindVariableFeatures(data_current[[j]])
      print(paste("Complete: ", rna_pheno[j,"SID"], "..", sep = ""))
    }

    if(length(data_current) > 1){
      data_current <- lapply(data_current, function(x) {
        x <- ScaleData(x, verbose=FALSE, features = VariableFeatures(x))
        x <- RunPCA(x, npcs = 30, verbose = FALSE, features = VariableFeatures(x))
      })

      if(toupper(integration_method) == "SEURAT" | is.null(integration_method)){
        red_method <- "pca"
        integrative_features <- SelectIntegrationFeatures(object.list = data_current)
        data_anchors <- FindIntegrationAnchors(object.list = data_current,
                                               reduction = "rpca", anchor.features = integrative_features)
        data <- IntegrateData(anchorset = data_anchors)
        DefaultAssay(data) <- "integrated"
        rm(data_anchors)
      }else if(toupper(integration_method) == "HARMONY"){
        data <- merge(data_current[[1]], data_current[c(2:length(data_current))]) # , add.cell.ids = names(data_current)
        data <- NormalizeData(data, verbose = FALSE) %>%
          FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
          ScaleData(verbose = FALSE) %>%
          RunPCA(pc.genes = data@var.genes, npcs = 30, verbose = FALSE)
        data <- RunHarmony(data, group.by.vars = "BATCH")
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
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    data <- ScaleData(data, verbose = FALSE)
    data <- RunPCA(data, verbose = FALSE)
    data <- RunUMAP(data, reduction = red_method, dims = 1:ifelse(length(data@reductions[[red_method]]) < 30, length(data@reductions[[red_method]]), 30))
    # data <- RunTSNE(data, reduction = red_method, dims = 1:ifelse(length(data@reductions[[red_method]]) < 30, length(data@reductions[[red_method]]), 30), check_duplicates = FALSE)
    data <- FindNeighbors(data, reduction = red_method, dims = 1:ifelse(length(data@reductions[[red_method]]) < 30, length(data@reductions[[red_method]]), 30))
    data <- FindClusters(data, resolution = 0.8)
    cluster_colors <- gen_colors(color_conditions$tenx, length(unique(data$seurat_clusters)))
    names(cluster_colors) <- c(sort(as.numeric(as.character(unique(data$seurat_clusters)))))

    print("Running SingleR..")
    Idents(data) <- "seurat_clusters"
    DefaultAssay(data) <- "RNA"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(data)),
                       clusters =  data$seurat_clusters,
                       ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.main)
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
        data@meta.data[which(row.names(data@meta.data) %in% current[[j]]$barcode),"DATA_TYPE"] <- unique(current[[j]][,which(toupper(colnames(current[[j]])) %in% toupper("cellType"))])[1]
      }
      current_contigs <- c(current_contigs, current)
    }
    data <- combineExpression(current_contigs, data, cloneCall="gene+nt", proportion = TRUE)
    # table(data$cloneType)
    data@meta.data[which(!toupper(data$DATA_TYPE) %in% toupper(c("BCR","TCR","B","T","T-AB","TGD"))),"DATA_TYPE"] <- "OTHERS"

    current_groups <- c("Hyperexpanded (0.1 < X <= 1)",
                        "Large (0.01 < X <= 0.1)",
                        "Medium (0.001 < X <= 0.01)",
                        "Small (1e-04 < X <= 0.001)",
                        "Rare (0 < X <= 1e-04)", NA)
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
    p1 <- occupiedscRepertoire(data, x.axis = "CELL_TYPE", proportion = F)+
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

    p2 <- occupiedscRepertoire(data, x.axis = "CELL_TYPE", proportion = T)+
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
    p36plots$GENE_NT <- clonalDiversity(combined_meta, cloneCall = "gene+nt") +
      ggtitle("CLONAL DIVERISTY: VDJC GENE + CDR3 NUCLEOTIDE") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    p36plots$GENE <- clonalDiversity(combined_meta, cloneCall = "gene") +
      ggtitle("CLONAL DIVERISTY: VDJC GENE") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    p36plots$NT <- clonalDiversity(combined_meta, cloneCall = "nt") +
      ggtitle("CLONAL DIVERISTY: CDR3 NUCLEOTIDE") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    p36plots$AA <- clonalDiversity(combined_meta, cloneCall = "aa") +
      ggtitle("CLONAL DIVERISTY: CDR3 AMINO ACID") +
      guides(fill=guide_legend(title="CELL TYPES")) +
      theme(plot.title = element_text(size=15, face = "bold"))

    Idents(data) <- "seurat_clusters"
    p37plots <- clonalOverlay(data, reduction = "umap",
                              freq.cutpoint = 0, bins = 10, facet = "GROUP") +
      guides(color = "none")+ scale_color_manual(values = cluster_colors)
    p37plots <- adjust_theme(p37plots, strip_size = 25)

    Idents(data) <- "CELL_TYPE"
    p38plots <- clonalOverlay(data, reduction = "umap",
                              freq.cutpoint = 0, bins = 10, facet = "GROUP") +
      guides(color = "none")+ scale_color_manual(values = celltype_colors)
    p38plots <- adjust_theme(p38plots, strip_size = 25)

    p39plots <- clonalNetwork(data,
                              reduction = "umap",
                              identity = "seurat_clusters",
                              filter.clones = NULL,
                              filter.identity = NULL,
                              cloneCall = "aa") + scale_color_manual(values = cluster_colors)
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

  somePDFPath = paste(cdir,"10URSA_PLOT_scIMMUNE_GENE_USAGE_JS-DIVERGENCE_CORRELATION_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  for(i in 1:length(p10plots)){
    plotx <- p10plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"11URSA_PLOT_scIMMUNE_HCLUSTERING_KMEANS_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=12, height=5.3,pointsize=12)
  for(i in 1:length(p11plots)){
    plotx <- p11plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"12URSA_PLOT_scIMMUNE_GENE_WISE_SAMPLE_DISTANCE_PCA_MDS_",project_name,".pdf", sep = "")
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
    kmers <- getKmers(plotx, i)
    kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
    kmers <- kmers[order(kmers[,2], decreasing = T),]
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
      kmers <- getKmers(plotx, i)
      kmers <- kmers[grep(";", kmers$Kmer, ignore.case = T, invert = T),]
      kmers <- kmers[order(kmers[,2], decreasing = T),]
      kmers_aa_stats <- kmer_profile(kmers)
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
  title(paste(project_name,": Repertoire Overlaps (CDR3 AA)", sep = ""), cex = 15)
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

  somePDFPath = paste(cdir,"20URSA_PLOT_scIMMUNE_CDR3_CONTIG_LENGTH_AA_NT_",project_name,".pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=10, pointsize=12)
  for(i in 1:length(p20plots)){
    plotx <- p20plots[[i]]
    print(plotx)
  }
  dev.off()

  somePDFPath = paste(cdir,"21URSA_PLOT_scIMMUNE_ALLUVIAL_CLONOTYPES_PROPORTIONS_CDR3_AA_",project_name,".pdf", sep = "")
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



