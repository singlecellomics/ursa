#' @keywords internal

plot_ploidy <- function(results, output_dir){
  data <- results$data
  data_cell_stats <- results$cell_stats
  project_name <- results$project
  for(i in 1:length(data)){
  sample_id <- names(data)[i]
  current <- data[[i]]
  cell_stats <- data_cell_stats[which(data_cell_stats$Sample == sample_id),]
  cell_stats$id <- factor(cell_stats$id, levels = levels(current$id))
  print("Printing ploidy by chrom position plot ..")
  p1 <- ggplot(current, aes(range, id, color = Ploidy, fill = Ploidy))+ # [sample(1:nrow(data),size = 100, replace = F),]
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),
          plot.title = element_text(size=25, face = "bold")) +
    # legend.text=element_text(size=18)) + # hjust = 0.5
    # geom_point(alpha = 0.6)+
    geom_tile(size = 0.8)+
    xlab("Chromosome Positions") + ylab("") +
    scale_color_manual(values = results$color_schemes$color_heatmap)+
    scale_fill_manual(values = results$color_schemes$color_heatmap)+
    facet_grid(.~chrom, scales = "free", switch = "x", space = "free_x")+
    ggtitle(paste(sample_id, ", ", signif(length(which(current$Ploidy == "2"))/nrow(current)*100,3), "% diploid"))

  p2 <- ggplot(cell_stats)+
    geom_bar(mapping = aes(x = 1, y = id, fill = Cell_Ploidy), stat = "identity", width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm")) +
    scale_fill_manual(values = results$color_schemes$cell_ploidy_colors)

  legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))
  p <- plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1))

  somePNGPath <- paste(output_dir,"01SCA_PLOT_PLOIDY_INFO_CHROM_POSITION_",sample_id,"_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
  print(p)
  dev.off()

  }

}

plot_ploidy_binary <- function(results, output_dir){
  plotx <- results$binary_cnv
  p1 <- ggplot(plotx, aes(CNV_Event, cell_id, fill = Count))+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),
          plot.title = element_text(size=25, face = "bold"),
          axis.title = element_text(size=20, face = "bold"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    geom_tile()+
    xlab("CNV Events") + ylab("") +
    scale_fill_manual(values = c("black","#eec643"))

  p2 <- ggplot(plotx)+
    geom_bar(mapping = aes(x = 1, y = cell_id, fill = Cluster), stat = "identity", width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    scale_fill_manual(values = results$color_schemes$cluster_colors)

  legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))

  somePNGPath <- paste(output_dir,"5URSA_PLOT_scCNV_BINARY_CNV_EVENTS_BY_CLUSTERS_",results$project,".png", sep = "")
  png(somePNGPath, width = 3000, height =2000, units = "px", res = 300)
  print(plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1)))
  dev.off()

}

plot_selected_chroms <- function(results, output_dir){
  plotx <- results$selected_chrom
  plotx_bar <- results$chrom_bar

  p1 <- ggplot(plotx, aes(range, cell_id, color = Ploidy, fill = Ploidy))+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),
          plot.title = element_text(size=25, face = "bold"),
          axis.title = element_text(size=20, face = "bold"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    geom_tile(size = 0.8)+
    scale_fill_manual(values = results$color_schemes$color_heatmap)+
    scale_color_manual(values = results$color_schemes$color_heatmap)+
    xlab("Chromosome Positions") + ylab("") +
    facet_grid(.~chrom, scales = "free", switch = "x", space = "free_x")

  p2 <- ggplot(plotx_bar)+
    geom_bar(mapping = aes(x = 1, y = cell_id, fill = Cluster), stat = "identity", width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    scale_fill_manual(values = results$color_schemes$cluster_colors)

  legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))

  somePNGPath <- paste(output_dir,"6URSA_PLOT_scCNV_PLOIDY_INFO_TOP_50_CNVs_",results$project,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
  print(plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1)))
  dev.off()
}

adjust_theme <- function(p, xangle = 0,legend = "right", title_size = 20, xsize=20, hejust = 0, vejust = 0, strip_size = 20, legend_title = 20, legend_text = 15){
  p <- p+ theme_classic()+
    # ggtitle(title) +
    # scale_fill_manual(values = col) +
    theme(axis.text.x = element_text(size = xsize, angle = xangle, hjust=hejust,vjust = vejust),
          axis.text.y = element_text(size = 20), legend.position = legend,
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = strip_size),
          strip.background = element_blank(),
          plot.title = element_text(size =title_size, face = "bold", hjust = 0.5))
  return(p)
}

plot_bygroup <- function(plotx, x, y, group, plot_title, col = NULL, annot = TRUE, legend_position = "right", numeric = FALSE,point_size = 4, label_size = 5, legendsize = 20){

  color_conditions <- color_ini()
  n <- length(unique(plotx[,group]))
  if(is.null(names(col)) | n > length(col)){
    # if(n > 1){
    col <- prepare_cols(col, n)
    # }
  }else{
    if(!is.null(names(col))){
      col <- col[which(names(col) %in% unique(plotx[,group]))]
    }else{
      col <- col[1:length(unique(plotx[,group]))]
    }
  }

  if(annot == TRUE){
    if(numeric == FALSE){
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.character(plotx[,group]))))
    }else{
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.numeric(as.character(plotx[,group])))))
    }
    centroids <- create_centroids(plotx, x, y, group)
    centroids$Col <- col

    plotx$Label <- ""
    for(i in 1:nrow(centroids)){
      plotx[nrow(plotx)+1,x] <- centroids[i,"Centroid_X"]
      plotx[nrow(plotx),y] <- centroids[i,"Centroid_Y"]
      plotx[nrow(plotx),group] <- centroids[i,"Cluster"]
      plotx[nrow(plotx),"Label"] <- centroids[i,"Cluster"]
    }

    p <- ggplot(plotx, aes(x = plotx[,x], y = plotx[,y], color = plotx[,group], label = Label)) +
      geom_point(alpha = ifelse(plotx$Label != "", 0, 1), size = point_size) +
      scale_color_manual(values = col) +
      theme_classic() +
      xlab(x) + ylab(y) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 20),
            legend.position = legend_position,
            legend.title = element_blank(),
            legend.text = element_text(size = legendsize),
            legend.key.size = unit(0.8, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(color=guide_legend(override.aes = list(size = 3), title=toupper(group), ncol=ifelse(length(unique(plotx[,group])) > 20, 2, 1)))

    if(numeric == FALSE){
      p <- p + geom_text_repel(max.overlaps = Inf,force=1,
                               point.padding = 0, # additional padding around each point
                               min.segment.length = 0, # draw all line segments
                               max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
                               box.padding = 0.3, size = label_size, colour = "black")
      # p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = label_size, hjust = 0, fontface =1)
    }else{
      p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = label_size, hjust = 0, fontface =1)
    }
  }else{
    p <- ggplot(plotx, aes(x = plotx[,x], y = plotx[,y], color = plotx[,group])) +
      geom_point(alpha = 1, size = point_size) +
      scale_color_manual(values = col) +
      theme_classic() +
      xlab(x) + ylab(y) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 25),
            legend.position = legend_position,
            legend.title = element_blank(),
            legend.text = element_text(size = legendsize),
            legend.key.size = unit(1, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(color=guide_legend(override.aes = list(size = 3), title=toupper(group), ncol=ifelse(length(unique(plotx[,group])) > 20, 2, 1)))
  }

  return(p)

}

create_centroids <- function(plotx, x, y, group){

  centroids <- split(plotx, plotx[,group])
  # centroid_groups <- names(centroids)
  centroids <- lapply(centroids, function(cent){
    cent <- data.frame(Cluster = ifelse(nrow(cent) > 0,as.character(unique(cent[,group])), NA),
                       Centroid_X = ifelse(nrow(cent) > 0, median(cent[,x]), NA),
                       Centroid_Y = ifelse(nrow(cent) > 0, median(cent[,y]), NA))})
  centroids <- do.call(rbind.data.frame, centroids)
  # centroids[,group] <- centroid_groups
  centroids <- centroids[!is.na(centroids[,"Cluster"]),]
  return(centroids)

}

#' @import ComplexHeatmap
complex_heatmap <- function(x, legend_title = NULL, col_title = NULL, row_title = NULL,
                            col = rev(rainbow(10)),nrow_clus = 1, row_annot = NULL, legendtitle = "EXPRESSION"){
  p <- Heatmap(x, name = legend_title, col = col, row_km = nrow_clus, right_annotation = row_annot,
               row_names_gp = gpar(fontsize = 8), row_names_side = "left",
               column_title = col_title,column_names_rot = 45,
               heatmap_legend_param = list(title = legendtitle))
  return(p)
}

violin_plot <- function(current, features, ncol = NULL, col, x_lab, log_status = TRUE){
  p <- VlnPlot(object = current, features = features, ncol = ncol, cols = col, pt.size = 0.05, log = log_status)  & xlab(x_lab)
  return(p)
}

own_facet_scatter <- function(plotx, feature1, feature2, isfacet = T,point_size = 0.9,
                              title, col = NULL, color_by,
                              group_by = NULL, xlabel = feature1, ylabel = feature2,
                              strip_size = 15, legend_pos = "right", ncol = 2){
  # set.seed(2022)
  if(is.null(col)){
  col <- color_ini()$bright
  }
  p <- ggplot(plotx, aes(x = plotx[,feature1], y = plotx[,feature2], color = plotx[,color_by])) +
    geom_point(size = point_size) +
    scale_color_manual(values = gen_colors(col, length(unique(plotx[,color_by])))) +
    # scale_color_manual(values = sample(gen_colors(col, length(unique(plotx[,color_by]))), size = length(unique(plotx[,color_by])), replace = F)) +
    theme_classic() +
    xlab(xlabel)+ylab(ylabel)+
    ggtitle(title) +
    guides(colour = guide_legend(title = color_by, override.aes = list(size=5)))+
    theme(legend.position = legend_pos,
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.8, "cm"),
          plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
          strip.text = element_text(size = strip_size, face = "bold"),
          strip.background = element_blank())
  if(isfacet == T & !is.null(group_by)){
    p <- p + facet_wrap(~plotx[,group_by], ncol = ncol)
  }
  return(p)
}

add_names <- function(current, sample_name,current_ident){
  sample_name <- gsub("_FILTERED_PEAK_BC|_FILTERED_FEATURE_BC","", sample_name, ignore.case = T)
  # sample_name <- ifelse(nchar(sample_name) > 8, gsub(".*\\.(.*)","\\1",sample_name), sample_name)
  current <- SetIdent(current, value = current_ident)
  levels(current@active.ident) <- sample_name
  current@project.name <- sample_name
  current$orig.ident <- sample_name
  return(current)
}

plot_slice <- function(plotx, plot_title, col = c(NULL,"default"), is.facet = FALSE, annot = TRUE, strip_size = 15, pt_size = 2){

  color_conditions <- color_ini()
  n <- length(unique(plotx$Cluster))
  col <- prepare_cols(col, n)

  # plotx$Cluster <- factor(plotx$Cluster, levels = sort(unique(as.character(plotx$Cluster))))

  if(annot == TRUE){
    centroids <- create_centroids(plotx, "X","Y", "Cluster")
    centroids$Col <- col
  }

  if(is.facet == TRUE){
    p <- ggplot(plotx, aes(X, Y, color = Cluster)) +
      # geom_point(shape = 15, aes(colour=Cluster)) +
      # geom_point(shape = 0,colour = "black", stroke = 0.05)+
      geom_point(shape = 16, size = pt_size) +
      coord_flip() + scale_x_reverse() + theme(legend.position = "right") +
      scale_color_manual(values = col) + facet_wrap(~Cluster) +
      theme_map() + guides(color = FALSE, size = FALSE) + ggtitle("") +
      theme(strip.text.x = element_text(face="bold", size = strip_size),
            strip.background = element_rect(colour="black", fill="white"))
  }else{
    p <- ggplot(plotx, aes(X, Y, color = Cluster)) +
      # geom_point(shape = 20, colour = "black", size = 1, stroke = 0.5) +
      geom_point(shape = 16, size = pt_size) +
      coord_flip() + scale_x_reverse() + theme(legend.position = "right") +
      scale_color_manual(values = col) +
      theme_map() + guides(color=guide_legend(title="Clusters", ncol=ifelse(length(unique(plotx$Cluster)) > 8, 2, 1))) +
      ggtitle(plot_title) + theme(plot.title = element_text(size=15, face = "bold", hjust = 0.5), legend.position = "right")
  }

  if(annot == TRUE){
    p <- p + guides(color = FALSE, size = FALSE) +
      annotate(geom="label", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = 3, fill=centroids$Col, col = "white")
  }
  return(p)

}

gen10x_plotx <- function(data, groups = NULL, selected = "ALL", include_meta = FALSE){
  if(toupper(selected[1]) == "ALL"){
    plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                        tSNE_1 = data@reductions$tsne@cell.embeddings[,"tSNE_1"],
                        tSNE_2 = data@reductions$tsne@cell.embeddings[,"tSNE_2"],
                        PC_1 = data@reductions$pca@cell.embeddings[,"PC_1"],
                        PC_2 = data@reductions$pca@cell.embeddings[,"PC_2"])
  }else{
    for(i in 1:length(selected)){
      if(i == 1){
        if(toupper(selected[i]) == "UMAP_SELECTED"){
          plotx <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("_SELECTED","",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }else{
          plotx <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("PCA|PCA_SELECTED","PC",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }
        colnames(plotx) <- paste(toupper(gsub("PCA","PC",selected[i], ignore.case = T)),c("_1","_2"), sep = "")
      }else{
        if(toupper(selected[i]) == "UMAP_SELECTED"){
          temp <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("_SELECTED","",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }else{
          temp <- data.frame(data@reductions[[tolower(selected[i])]]@cell.embeddings[,grep(paste(paste(gsub("PCA|PCA_SELECTED","PC",selected[i], ignore.case = T),c("_1$","_2$"), sep = ""), collapse = "|"), colnames(data@reductions[[tolower(selected[i])]]@cell.embeddings), ignore.case = T)])
        }
        colnames(temp) <- paste(toupper(gsub("PCA","PC",selected[i], ignore.case = T)),c("_1","_2"), sep = "")
        plotx <- cbind(plotx, temp)
      }
    }
    colnames(plotx) <- toupper(gsub("PCA_([0-9]+)$","PC_\\1",colnames(plotx), ignore.case = T))
    colnames(plotx) <- toupper(gsub("tSNE","tSNE",colnames(plotx), ignore.case = T))
  }

  if(!is.null(groups)){
    for(i in 1:length(groups)){
      plotx[,groups[i]] <- data@meta.data[,groups[i]]
    }
  }

  if(include_meta == TRUE){
    plotx <- cbind(data@meta.data, plotx)
  }

  return(plotx)

}

own_violin <- function(plotx, x = "SAMPLE_ID", feature, plotx_title, col = NULL, title.size = 20, angle = 0, hjust=NULL,vjust = NULL){
  p <- ggplot(plotx, aes_string(x=x, y=feature, fill = x)) +
    geom_violin(trim=TRUE) + scale_fill_manual(values = col)+
    theme_classic()+
    # scale_y_continuous(trans='log10') +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = angle, size = ifelse(nchar(as.character(unique(plotx$SAMPLE_ID))) > 20, 12,15), hjust=hjust,vjust = vjust),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = title.size, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 15, margin=margin(10,0,0,0),face = "bold")) +
    xlab("SAMPLE_ID") + ylab("") + ggtitle(plotx_title)
  if(max(plotx[,feature]) > 0){
    p <- p + geom_jitter(size = 0.05)
  }

  return(p)
}

own_feature <- function(plotx, feature1, feature2, title_name, col, xlabel, ylabel){
  p <- ggplot(plotx, aes(plotx[,feature1], plotx[,feature2], color = SAMPLE_ID))+
    geom_point()+theme_classic()+
    scale_color_manual(values = gen_colors(col, length(unique(plotx[,"SAMPLE_ID"])))) +
    stat_cor() +
    xlab(xlabel)+ylab(ylabel)+ # ggtitle(title_name)+
    theme(legend.position = "none", axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 15, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 15, margin=margin(0,10,0,0)),
          plot.margin=unit(c(1,1,1,1),"cm"))
  # plot.title = element_text(size = title.size, face = "bold", hjust = 0.5))
  return(p)
}

own_2d_scatter <- function(current_data, reduction_method, split_var, plot_title, cols = NULL){
  if(is.null(cols)){
    cols <- color_ini()
    cols <- cols$bright
  }
  p <- DimPlot(current_data, reduction = reduction_method,split.by = split_var, cols = cols) +
    ggtitle(paste(plot_title, sep = "")) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.2, "cm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 15, margin=margin(0,10,0,0)))
  return(p)
}

deanalysis <- function(data, current_clusters, plot_title,group=NULL,
                       de_analysis = "findallmarkers", n=10, cluster_name = "seurat_clusters"){

  DefaultAssay(data) <- "RNA"
  Idents(data) <- cluster_name
  out <- NULL
  current_data_markers <- NULL
  p_thresh <- 0.05
  # fc_threshold <- 1.25
  current_clusters <- sort(as.numeric(as.character(unique(current_clusters))), decreasing = F)
  top1 <- NULL
  topn <- NULL

  if(toupper(de_analysis) == toupper("findallmarkers")){
    de_type <- "TOP_DEGENES_IN_CLUSTERS"
    de_name <- ""
    current_data_markers <- FindAllMarkers(data, min.pct = 0.5, logfc.threshold = 0.25)
    current_data_markers <- current_data_markers[current_data_markers$p_val_adj < p_thresh,]
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]

    wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
    # temp <- current_data_markers[current_data_markers[,wt] > fc_threshold,]
    top1 <- choose_topn(current_data_markers, wt = wt, group = "cluster", n = 1, max = T)
    # wt <- colnames(current_data_markers)[grep("p_val_adj", colnames(current_data_markers), ignore.case = T)]
    # top1 <- rbind(top1,choose_topn(current_data_markers, wt = wt, group = "cluster", n = 1, max = F))
    top1 <- top1[which(!is.na(top1$gene)),]
    top1 <- unique(top1)

    wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
    topn <- choose_topn(current_data_markers, wt = wt, group = "cluster", n = n, max = T)
    # wt <- colnames(current_data_markers)[grep("p_val_adj", colnames(current_data_markers), ignore.case = T)]
    # topn <- rbind(topn,choose_topn(current_data_markers, wt = wt, group = "cluster", n = n, max = F))
    # topn <- unique(topn)

    current_data_markers <- current_data_markers[order(current_data_markers$cluster, decreasing = F),]
    top1 <- top1[order(top1$cluster, decreasing = F),]
    top1 <- unique(top1)
    topn <- topn[order(topn$cluster, decreasing = F),]
    topn <- unique(topn)

  }else if(toupper(de_analysis) == toupper("findconservedmarkers")){
    de_type <- "TOP_CONSERVED_GENES_GROUPS_IN_CLUSTERS"
    de_name <- "CONSERVED MARKERS: "

    for(i in 1:length(current_clusters)){
      if(min(table(data@meta.data[,group][Idents(data) == current_clusters[i]])) > 3){
        print(i)
        current <- FindConservedMarkers(data, ident.1 = current_clusters[i], grouping.var = group, verbose = FALSE)
        current$gene <- row.names(current)
        current$cluster <- current_clusters[i]
        current <- current[current$minimump_p_val < p_thresh,]
        if(nrow(current) > 0){
          current_data_markers <- rbind(current_data_markers, current)
          current_groups <- colnames(current)[grep("avg_log.*FC", colnames(current), ignore.case = T)]
          top1 <- rbind(top1,current[unique(apply(current[,current_groups],2,function(x){which.max(x)})),])
          for(m in 1:length(current_groups)){
            topn <- rbind(topn,current[order(current[,current_groups[m]], decreasing = T),][1:n,])
          }
          # current_groups <- colnames(current)[grep("_p_val_adj", colnames(current), ignore.case = T)]
          # top1 <- rbind(top1,current[unique(apply(current[,current_groups],2,function(x){which.min(x)})),])
          # for(m in 1:length(current_groups)){
          #   topn <- rbind(topn,current[order(current[,current_groups[m]], decreasing = F),][1:n,])
          # }
          top1 <- unique(top1)
          topn <- unique(topn)
        }
      }
    }
    current_data_markers <- current_data_markers[order(current_data_markers$minimump_p_val, decreasing = F),]
    top1 <- top1[order(top1$cluster, decreasing = F),]
    top1 <- unique(top1)
    topn <- topn[order(topn$cluster, decreasing = F),]
    topn <- unique(topn)

  }else if(toupper(de_analysis) == toupper("finddegroup")){
    de_type <- "TOP_DE_GENES_GROUPS_WITHIN_EACH_CLUSTER"
    de_name <- "DE BETWEEN GROUPS: "

    for(i in 1:length(current_clusters)){
      current <-  subset(data, idents = current_clusters[i])
      if(min(table(current$GROUP)) > 3){
        Idents(current) <- group
        current_groups <- unique(current@meta.data[,group])
        for(l in 1:length(current_groups)){
          for(m in 1:length(current_groups)){
            if(l != m){
              current_result <- FindMarkers(current, ident.1 = current_groups[l], ident.2 = current_groups[m])
              current_result <- current_result[current_result$p_val_adj < p_thresh,]
              current_result <- current_result[order(current_result$p_val_adj, decreasing = F),]
              current_result$gene <- row.names(current_result)
              current_result$cluster <- current_clusters[i]
              current_result$group1 <- current_groups[l]
              current_result$group2 <- current_groups[m]

              wt <- colnames(current_result)[grep("log.*FC", colnames(current_result), ignore.case = T)]
              top1 <- rbind(top1, current_result[which.max(current_result[,wt]),])
              # wt <- colnames(current_result)[grep("p_val_adj", colnames(current_result), ignore.case = T)]
              # top1 <- rbind(top1, current_result[which.min(current_result[,wt]),])

              wt <- colnames(current_result)[grep("log.*FC", colnames(current_result), ignore.case = T)]
              current_result <- current_result[order(current_result[,wt], decreasing = T),]
              topn <- rbind(topn, current_result[1:n,])
              # wt <- colnames(current_result)[grep("p_val_adj", colnames(current_result), ignore.case = T)]
              # current_result <- current_result[order(current_result[,wt], decreasing = F),]
              # topn <- rbind(topn, current_result[1:n,])

              top1 <- unique(top1)
              topn <- unique(topn)
              current_data_markers <- rbind(current_result, current_data_markers)
            }
          }
        }
      }
    }
    current_data_markers <- current_data_markers[order(current_data_markers$p_val_adj, decreasing = F),]
    top1 <- top1[order(top1$cluster, decreasing = F),]
    top1 <- unique(top1)
    topn <- topn[order(topn$cluster, decreasing = F),]
    topn <- unique(topn)
  }

  # findallmarkers
  # findconservedmarkers
  # finddegroup

  out$current_data_markers <- current_data_markers
  p <- NULL

  DefaultAssay(data) <- "RNA"
  for(k in 1:length(top1$gene)){
    if(toupper(de_analysis) %in% c(toupper("finddegroup"), toupper("findconservedmarkers"))){
      if(toupper(de_analysis) == toupper("finddegroup")){
        extra_info <- paste(top1$group1[k], " VS ", top1$group2[k], sep = "")
      }else{
        extra_info <- ""
      }
      p[[k]] <- FeaturePlot(data, features = top1$gene[k], split.by = group,
                            label = T, pt.size = 0.5, label.size = 6, max.cutoff = 'q95')+
        plot_annotation(title = paste(extra_info,"\nCLUSTER ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""),
                        theme =
                          theme(axis.text.x = element_text(size = 20),
                                axis.text.y = element_text(size = 20),
                                axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                legend.title = element_blank(),
                                legend.text = element_text(size = 15),
                                plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
      if(toupper(de_analysis) == toupper("finddegroup")){
        names(p)[k] <- paste("CLUSTER ",top1$cluster[k],": TOP GENE IN ", top1$group1[k],sep = "")
      }else{
        names(p)[k] <- top1$cluster[k]
      }
    }else{
      extra_info <- ""
      if(k == 1){
        p <- FeaturePlot(data, features = top1$gene[k], split.by = group,
                         label = T, pt.size = 0.5, label.size = 6, max.cutoff = 'q95')+
          plot_annotation(title = paste(extra_info,"\nCLUSTER ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""),
                          theme =
                            theme(axis.text.x = element_text(size = 20),
                                  axis.text.y = element_text(size = 20),
                                  axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                  axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 15),
                                  plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
      }else{
        p1 <- FeaturePlot(data, features = top1$gene[k], split.by = group,
                          label = T, pt.size = 0.5, label.size = 6, max.cutoff = 'q95')+
          plot_annotation(title = paste(extra_info,"\nCLUSTER ", top1$cluster[k], ", GENE: ",top1$gene[k],sep = ""),
                          theme =
                            theme(axis.text.x = element_text(size = 20),
                                  axis.text.y = element_text(size = 20),
                                  axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
                                  axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 15),
                                  plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))
        p <- p+p1
      }
    }
  }

  out$de_type <- de_type
  out$featureplot <- p
  out$top1 <- top1
  out$topn <- topn
  out$de_name <- de_name

  return(out)
}

choose_topn <- function(current, wt, group, n, max = T){
  current <- current[order(current[,wt], decreasing = ifelse(max == T, T, F)),]
  temp <- split(current,current[,group])
  temp <- lapply(temp, function(x){x <- x[1:n,]})
  temp <- do.call(rbind.data.frame, temp)
  return(temp)
}

group_medianexpr <- function(current_data_markers, data, ref_group = "seurat_clusters", group = "seurat_clusters", cell_type = F,n=10){

  wt <- colnames(current_data_markers)[grep("log.*FC", colnames(current_data_markers), ignore.case = T)]
  if(cell_type == TRUE){
    current_data_markers$CELL_TYPE <- data@meta.data[match(current_data_markers$cluster, data@meta.data[,ref_group]),"CELL_TYPE"]
    top_markers <- current_data_markers %>% group_by(CELL_TYPE) %>% top_n(n = n, eval(parse(text = wt)))
  }else{
    top_markers <- current_data_markers %>% group_by(cluster) %>% top_n(n = n, eval(parse(text = wt)))
  }
  top_markers <- top_markers[order(top_markers[,grep("log.*FC", colnames(current_data_markers), ignore.case = T)], decreasing = T),]
  top_markers$gene <- factor(top_markers$gene, levels = c(as.character(unlist(unique(top_markers[order(top_markers[,grep("log.*FC", colnames(current_data_markers), ignore.case = T)], decreasing = T),"gene"])))))
  # top_markers[,group] <- factor(top_markers[,group], levels = sort(as.numeric(as.character(unique(top_markers[,group])), decreasing = F)))

  top_expr <- data.frame(t(as.matrix(GetAssayData(data))[which(row.names(data) %in% unique(top_markers$gene)),]))

  if(cell_type == FALSE){
    top_expr$population <- data@meta.data[match(row.names(top_expr), row.names(data@meta.data)),group]
    clusters <- sort(as.numeric(as.character(unique(top_expr$population))))
  }else{
    top_expr$population <- data@meta.data[match(row.names(top_expr), row.names(data@meta.data)),"CELL_TYPE"]
    clusters <- sort(as.character(unique(top_expr$population)))
  }

  median_expr <- NULL
  cell_number <- NULL
  expr <- NULL

  for(i in 1:length(clusters)){
    current <- top_expr[which(top_expr$population == clusters[i]),grep("population",colnames(top_expr), ignore.case = T, invert = T)]
    cell_number <- c(cell_number,nrow(current))
    for (k in 1:ncol(current)){
      expr <- c(expr,median(current[current[,k]>0,k]))
    }
    median_expr <- rbind(median_expr, expr)
    expr <- NULL
  }

  median_expr <- data.frame(t(median_expr))
  row.names(median_expr) <- colnames(top_expr)[grep("population",colnames(top_expr), ignore.case = T, invert = T)]
  colnames(median_expr) <- clusters
  # median_expr <- median_expr[,colSums(median_expr) > 0]

  plot_median <- scale(median_expr)
  plot_median <- t(scale(t(plot_median)))
  out <- NULL
  out$plot_median <- plot_median
  out$top_markers <- top_markers
  out$top_expr <- top_expr
  # out$current_data_markers
  return(out)
}

beforeafter_dimplot <- function(data1, data2, dim1, dim2, group, subtitle1, subtitle2, maintitle, titlesize, col, legend_position = "bottom"){
  p1 <- plot_bygroup(data1, x = dim1, y = dim2, group = group, plot_title = subtitle1,
                     col = col, annot = FALSE, legend_position = legend_position, point_size = 1)
  p2 <- plot_bygroup(data2,  x = dim1, y = dim2, group = group, plot_title = subtitle2,
                     col = col, annot = FALSE, legend_position = legend_position, point_size = 1)

  p <- p1+p2+
    plot_annotation(title = maintitle, theme = theme(plot.title = element_text(size = titlesize, face = "bold", hjust = 0.5)))

  return(p)
}

plot_pseudo <- function(data, reduction, group, label_size, plot_title, col, n_col,
                        cell_size = 2, traj_size = 1.5){

  p <- plot_cells(data,
                  reduction_method = reduction,
                  color_cells_by = group,
                  group_label_size = label_size,
                  graph_label_size = 6,
                  cell_size = cell_size,
                  cell_stroke = I(2/2),
                  alpha = 0.9,
                  trajectory_graph_segment_size = traj_size,
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE)+
    scale_color_manual(values = gen_colors(col,n_col))+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 15),
          # legend.position = "right",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.2, "cm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

  return(p)
}

# get_earliest_principal_node() is from Monocle3
get_earliest_principal_node <- function(mono3data, group, group_element){
  cell_ids <- which(colData(mono3data)[,group] == group_element)
  # mono3data <- preprocess_cds(mono3data, num_dim = 100)
  closest_vertex <-
    mono3data@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(mono3data),])
  root_pr_nodes <-
    igraph::V(principal_graph(mono3data)[["UMAP"]])$name[as.numeric(names
                                                                    (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

adjust_plot <- function(p,col,n_col,plot_title="", fill = F){
  if(fill == F){
    p <- p + scale_color_manual(values = gen_colors(col,n_col))
  }
  p <- p +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 15),
          # legend.position = "right",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.2, "cm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)))+
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  return(p)
}

#' @import plotly
plot_3d <- function(d1, d2, d3, color_group, plot_title, n, lt,
                    current_text,t1,t2,t3,col = NULL){
  if(is.null(col)){
    col <- color_ini()$general
  }
  p <- plot_ly(x=d1,y=d2,z=d3,type="scatter3d", mode="markers",size = 0.2,
               colors = gen_colors(col, n),
               color=color_group,
               hoverinfo = "text",
               hovertext = current_text)%>%
    layout(scene = list(xaxis = list(title = t1),
                        yaxis = list(title = t2),
                        zaxis = list(title = t3)),
           title =paste("<b>",plot_title,"</b>",sep = ""),
           legend = list(title = list(text = lt)))
  return(p)

}
