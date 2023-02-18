#' @keywords internal

time_ini <- function(){
  x <- format(Sys.time(), format = "%Y%m%d%H%M%S")
  return(x)
}

gen_colors <- function(col = NULL, n, lum = "light"){
  if(is.null(col)){
    set.seed(9)
    col <- randomcoloR::distinctColorPalette(k = n, runTsne = T)
  }else{
    col <- colorRampPalette(col)(n)
  }
  return(col)
}

color_ini <- function(){
  color_conditions <- NULL
  color_conditions$monet <- c("#76A2BB","#BCB180","#E8AC75","#EE9873","#D69896","#F4AFA3","#1166AE","#1C9FD3","#C45A35","#9C9B0B","#E34615","#F6E576","#1A3D40","#7F2125","#1E6652","#1BA0C7","#075B25","#0F8441","#F6C04A","#E1B5B2")
  color_conditions$vanngogh <- c("#EEE6B2","#F9E956","#A66100","#54ABAB","#4062B4","#6D87B9","#E6AE95","#5D1C01","#964000","#010101","#F15E01","#EB9D00","#5B1700","#6E1700","#8A7B01","#DED6AD")
  color_conditions$manycolors <- c(unique(unlist(lapply(ggthemes::ggthemes_data$tableau$`color-palettes`$regular,function(x){x <- x$value}))))
  color_conditions$ggplot <- c("#F8766D", "#EA8331", "#D89000", "#C09B00", "#A3A500", "#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3", "#00BFC4", "#00BAE0", "#00B0F6", "#35A2FF", "#9590FF", "#C77CFF", "#E76BF3", "#FA62DB", "#FF62BC", "#FF6A98")
  color_conditions$colorful <- c("#EC6077", "#9AF270", "#CBD157", "#F0914C", "#6745A4", "#70E695", "#56B6CD", "#5076D7", "#B94AA7")
  color_conditions$general <- c("#5F75A5","#91D376","#D75B58","#F5BFD3","#A8C9EA","#B09C16","#F69F99","#AC79A3","#E89314","#EAD256",
                                "#78706E","#D1A5CA","#F7C277","#569794","#B9B0AC","#99785E","#5FA346","#8DBCB6","#CC7296","#D3B6A5")
  color_conditions$tableau20 <- c("#9edae5","#17becf","#dbdb8d","#bcbd22","#c7c7c7","#7f7f7f","#f7b6d2","#e377c2","#c49c94","#8c564b",
                                  "#c5b0d5","#9467bd","#ff9896","#d62728","#98df8a","#2ca02c","#ffbb78","#ff7f0e","#aec7e8","#1f77b4")
  color_conditions$tenx <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold", "#a65628", "#999999", "#9e2a2b", "grey", "#58bc82", "purple")
  color_conditions$mark <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                             "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                             "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                             "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                             "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                             "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d")
  color_conditions$bright <- c("#6C85B6","#F5B317","#E57369","#94C1CB","#61AE87","#BFC029","#B67694")
  color_conditions$cold <- c("#5F75A5","#91d376","#A8C9EA","#78706E","#569794","#8DBCB6","#5FA346","#B9B0AC")
  color_conditions$warm <- c("#D75B58","#E89314","#F5BFD3","#F69F99","#F7C277","#CC7296","#D3B6A5","#D1A5CA","#99785E","#B09C16","#EAD256","#AC79A3")
  color_conditions$alternate <- c("#5F75A5","#E89314","#91d376","#D75B58","#A8C9EA","#F5BFD3","#78706E","#F69F99","#569794","#F7C277","#8DBCB6","#CC7296","#5FA346","#D3B6A5","#B9B0AC","#D1A5CA")
  color_conditions$gradient <- c("#584C8E","#5079A9","#5EA1AC","#77C1A1","#A2CF9F",
                                 "#CBDE9C","#E9ECA5","#FAE79E","#F7CC82","#EEA56C",
                                 "#E77A58","#D4534E","#B2334A","#8A2140")
  color_conditions$redbluecont <- c("#26456e","#1c5998","#1c73b1","#3a87b7","#67add4","#cacaca","#fc8375","#df513f","#d11719","#bd1316","#9c0824")
  color_conditions$orangebluecont <- c("#26456e","#1c5998","#1c73b1","#3a87b7","#67add4","#7bc8e2","#cacaca","#fdab67","#fd8938","#f06511","#d74401","#a33202","#7b3014")
  color_conditions$trafficlightcont <- c("#9fcd99","#ffdd71","#f26c64")
  color_conditions$orangegradient <- c("white", "cornflowerblue", "yellow", "red")
  color_conditions$flow <- c("blue","cyan","green","yellow","orange","red","red4")
  color_conditions$Rainbow <- rev(rainbow(10))
  color_conditions$OrangeBlue <- c("#26456e","#1c5998","#1c73b1","#3a87b7","#67add4","#7bc8e2","#cacaca","#fdab67", "#fd8938","#f06511","#d74401","#a33202","#7b3014")
  color_conditions$GreenRed <- c("green", "white", "red")
  color_conditions$BlueRed <- c("blue", "#EEEEEE", "red")
  color_conditions$Viridis <- hcl.colors(10)
  color_conditions$Heat <- heat.colors(10)
  color_conditions$Plasma <- hcl.colors(10, palette = "Plasma")
  color_conditions$Zissou <- hcl.colors(10, palette = "Zissou 1")
  color_conditions$RedYellowBlue <- hcl.colors(10, palette = "RdYlBu")
  color_conditions$Spectral <- hcl.colors(10, palette = "Spectral")
  color_conditions$Terrain <- terrain.colors(10)
  color_conditions$CM <- cm.colors(10)
  color_conditions$BlueYellowRed <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
  return(color_conditions)
}

#' @import xlsx
readfile <- function(file_dir){
  if(length(grep("\\.csv$", file_dir, ignore.case = T)) > 0){
    if(ncol(read.csv(file_dir, header = T, sep = ";")) > 1){
      current <- suppressWarnings(read.csv(file_dir, header = T, sep = ";"))
    }else{
      current <- suppressWarnings(read.csv(file_dir, header = T))
    }
  }else if(length(grep("\\.txt$", file_dir, ignore.case = T)) > 0){
    if(ncol(suppressWarnings(read.table(file_dir, header = T))) == 1){
      current <- suppressWarnings(read.table(file_dir, header = T, sep = ";"))
    } else if(ncol(suppressWarnings(read.table(file_dir, header = T, sep = ","))) > 1){
      current <- suppressWarnings(read.table(file_dir, header = T, sep = ","))
    }else{
      current <- suppressWarnings(read.table(file_dir, header = T))
    }
  }else if(length(grep("\\.xls$|\\.xlsx$", file_dir, ignore.case = T)) > 0){
    current <- suppressWarnings(read.xlsx(file_dir, 1, header=TRUE))
  }else{
    current <- suppressWarnings(read.table(file_dir, header=TRUE, sep = "\t"))
  }

  if(toupper(colnames(current)[1]) == "X"){
    row.names(current) <- current[,1]
    current <- current[,which(toupper(colnames(current)) != "X")]
  }
  return(current)
}

pheno_ini <- function(phenodata_dir, pipeline, extra_para = NULL, isDir = TRUE){

  if(isDir == TRUE){
    pheno_data <- readfile(phenodata_dir)
  }else{
    pheno_data <- phenodata_dir
  }

  if(toupper(pipeline) == "SPATIAL"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","GROUP")
  }else if (toupper(pipeline) == "SCRNASEQ"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP") # "CELL_CYCLE", "DATA_TYPE","CELL_TYPE","PAIR_ID","TISSUE_TYPE","CANCER_TYPE"
  } else if(toupper(pipeline) == "SMARTSEQ2"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
  }else if(toupper(pipeline) == "SCCNV"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID")
  }else if(toupper(pipeline) == toupper("scIMMUNE")){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP","DATA_TYPE","CELL_TYPE","PAIR_ID")
  }else if(toupper(pipeline) == "BULK"){
    if(length(grep("Treatment", colnames(pheno_data), ignore.case = TRUE)) > 0){
      colnames(pheno_data) <- c("SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP","TREATMENT")
    }else{
      colnames(pheno_data) <- c("SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
    }
  }else if(toupper(pipeline) == "BULK_TCR"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","CELL_TYPE","GROUP","BATCH")
  }else if (toupper(pipeline) == "FLOW"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","BATCH","GROUP")
  }else if (toupper(pipeline) == "SCATAC"){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","GROUP","REF_GENOME")
  }else if (toupper(pipeline) %in% c("EV","CYTOF")){
    colnames(pheno_data) <- c("FILE","SAMPLE_ID","INDIVIDUAL_ID","BATCH","GROUP")
  }
  if((length(grep("run|BATCH", pheno_data$BATCH, ignore.case = T)) == 0) &
     (length(grep("BATCH", colnames(pheno_data), ignore.case = T)) > 0) &
     (all(grepl("^[0-9]+$", pheno_data$BATCH)))){
    pheno_data$BATCH <- paste("BATCH_", pheno_data$BATCH, sep = "")
  }

  if((length(grep("SAMPLE|^S_|^S.*", pheno_data$SAMPLE_ID, ignore.case = T)) == 0) &
     (length(grep("SAMPLE_ID", colnames(pheno_data), ignore.case = T)) > 0) &
     (all(grepl("^[0-9]+$", pheno_data$SAMPLE_ID)))){
    pheno_data$SAMPLE_ID <- paste("SAMPLE_", pheno_data$SAMPLE_ID, sep = "")
  }

  if((length(grep("INDIVIDUAL_ID", colnames(pheno_data), ignore.case = T)) > 0) &
     (all(grepl("^[0-9]+$", pheno_data$INDIVIDUAL_ID)))){
    pheno_data$INDIVIDUAL_ID <- paste("P", pheno_data$INDIVIDUAL_ID, sep = "")
  }

  for(i in grep("file|ref_genome|RNASEQ_FILE", colnames(pheno_data), ignore.case = T, invert = T)){
    pheno_data[,i] <- gsub("^\\s+|\\s+$|\\.fcs|\\.csv|\\.h5|\\.txt|\\.mtx|\\.tsv|\\.gz","",pheno_data[,i], ignore.case = T)
    pheno_data[,i] <- gsub("\\s+","_",pheno_data[,i], ignore.case = T)
    pheno_data[,i] <- gsub("-","_",pheno_data[,i], ignore.case = T)
    pheno_data[,i] <- toupper(pheno_data[,i])
  }

  if(toupper(pipeline) == "BULK"){
    for(i in grep("SAMPLE_ID", colnames(pheno_data), ignore.case = T, invert = T)){
      pheno_data[,i] <- factor(pheno_data[,i])
    }
  }
  return(pheno_data)
}

ini_colorschemes <- function(assay = "scRNA", data){
  color_conditions <- color_ini()
  ccolor_schemes <- NULL
  if(assay == "scCNV"){
    ploidy_levels <- c("0","1","2","3","4","5","6","7","8","9",">=10")
    color_heatmap <- c("#306192","#548292","white","#7EB0C0","#9BBAA0","#B5C07F","#C8C667","#DEBD42","#D8AC3A","#D96829","#DE3C22")
    names(color_heatmap) <- ploidy_levels
    cell_ploidy <- c("diploid", "non-diploid", "noisy")
    cell_ploidy_colors <- c("#5F75A5", "#B0A062", "#D3B6A5")
    names(cell_ploidy_colors) <- cell_ploidy
    sample_colors <- gen_colors(color_conditions$tenx, length(unique(data$data_cell_stats$Sample)))
    names(sample_colors) <- unique(data$data_cell_stats$Sample)
    cluster_colors <- gen_colors(color_conditions$bright, length(unique(data$umap_coords$Cluster)))
    names(cluster_colors) <- sort(as.numeric(as.character(unique(data$umap_coords$Cluster))))
    ccolor_schemes$color_heatmap <- color_heatmap
    ccolor_schemes$cell_ploidy_colors <- cell_ploidy_colors
    ccolor_schemes$sample_colors <- sample_colors
    ccolor_schemes$cluster_colors <- cluster_colors
    return(ccolor_schemes)
  }
}

prepare_cols <- function(col, n){

  color_conditions <- color_ini()
  if(length(col) > 1 & !is.null(col)){
    col <- colorRampPalette(col)(n)
  } else if(is.null(col)){
    col <- colorRampPalette(color_conditions$colorful)(n)
  }else if(toupper(col) == toupper("default")){
    col <- gg_color_hue(n)
  }

}
