library(ggplot2)
library(patchwork)
dr_tsne <- function(data = data, tsneSeed = 42) {
  data <- as.matrix(data)
  if (is.numeric(tsneSeed)) set.seed(tsneSeed)
  mapped <- Rtsne::Rtsne(data, initial_dims = ncol(data), check_duplicates = FALSE, pca = TRUE, num_threads = parallel::detectCores()-1)$Y
  if (!is.null(mapped)) {
    colnames(mapped) <- paste("tsne", 1:2, sep = "_")
    rownames(mapped) <- row.names(data)
  }
  return(mapped)
}

ROC_plot <- function (i, test, case_problist) {
  res <- try(suppressMessages(pROC::roc(test[[i]]$condition, case_problist[[i]], auc = T)),
             silent = T)
  if (class(res) == "try-error") return(NULL)
  plot(res, main = paste("Cluster", i), col = '#e73667', lwd = 3,
       print.auc = TRUE, print.auc.col = 'black',
       auc.polygon = TRUE, max.auc.polygon = TRUE, auc.polygon.col = "#fef8f8",
       identity.lwd = 2, identity.col = 'grey',
       print.thres = TRUE)
}

mergedimred <- function(subdata_cluster_DEG) {
  vars <- c("filename", "condition", "cluster")
  if (is.null(subdata_cluster_DEG)) return(NULL)
  data <- dplyr::select(subdata_cluster_DEG, -dplyr::one_of(vars))
  data_dimred <- try(dr_tsne(data = data, tsneSeed = 42), silent = T)
  if (any(class(data_dimred) == "try-error")) return(NULL)
  dplyr::left_join(tibble::rownames_to_column(subdata_cluster_DEG),
                   tibble::rownames_to_column(as.data.frame(data_dimred)),
                   by = "rowname", copy = TRUE)
}

purity_plot <- function(i, subdata_cluster_dimred, sub_cluster_label) {
  data <- subdata_cluster_dimred[[i]]
  # if (is.null(data)) return(NULL)
  cluster_label <- as.factor(sub_cluster_label[[i]])
  data_condition <- data.frame(dplyr::select(data, dplyr::one_of(c("tsne_1", "tsne_2", "condition"))),
                               cluster = cluster_label)
  
  data_ell <- data.frame()
  colorset <- c('#fec600', '#810061')
  for (g in 1:length(unique(data_condition$cluster))) {
    data <- dplyr::filter(data_condition, cluster == g)
    elli <- as.data.frame(ellipse::ellipse(cor(data$tsne_1, data$tsne_2),
                                           scale = c(sd(data$tsne_1), sd(data$tsne_2)),
                                           centre =c(mean(data$tsne_1), mean(data$tsne_2))))
    colnames(elli) <- c("tsne_1", "tsne_2")
    elli$cluster <- as.factor(g)
    elli$color <- colorset[elli$cluster]
    data_ell <- rbind(data_ell, elli)
  }
  
  p_condition <- ggplot(data = data_condition, aes(tsne_1, tsne_2, color = condition)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(name = "Condition", values = c("#00ced1", "#F87168")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  
  p_cluster <- ggplot(data = data_condition, aes(tsne_1, tsne_2, color = cluster)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(name = "Cluster", values = c("#fec600", "#810061")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  # p_condition_cluster <- ggplot() +
  #   geom_point(data = data_condition,
  #              aes(tsne_1, tsne_2, color = condition),
  #              alpha = 0.5) +
  #   geom_path(data = data_ell,
  #             aes(tsne_1, tsne_2, color = cluster),
  #             size = 1.2, linetype = 2, inherit.aes = F, alpha = 0.8) + 
  #   scale_color_manual(name = "Condition & Cluster",
  #                      breaks = c(levels(data_condition$condition), "1", "2"),
  #                      labels = c(levels(data_condition$condition), "Cluster 1", "Cluster 2"),
  #                      values = c("#00ced1", "#F87168", "#fec600", "#810061"),
  #                      guide = guide_legend(nrow = 2)) +
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         legend.position = "top")
  
  # return(ggpubr::ggarrange(p_condition, p_cluster, p_condition_cluster, ncol = 3, labels = paste("cluster", i)))
  return(ggpubr::ggarrange(p_condition, p_cluster, ncol = 2, labels = paste("cluster", i)))
}


heatmap_plot <- function (i, subdata_cluster_DEG) {
  data <- subdata_cluster_DEG[[i]]
  # if (is.null(data) || (length(data) == 4)) return(NULL)
  data_pure <- dplyr::select(data, -dplyr::one_of(c("filename", "condition", "cluster")))
  data_heatmap <- aggregate(data_pure, by = list(data$filename, data$condition), mean)
  rownames(data_heatmap) <- data_heatmap[, 1]
  data_heatmap <- data_heatmap[, -c(1,2)]
  
  
  data_heatmap <- t(data_heatmap)
  pheatmap::pheatmap(
    data_heatmap,
    color = gplots::colorpanel(100, "#c10022", "#ffffbd", "#2c3298"),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste("Cluster", i),
    fontsize_row = 8,
    fontsize_col = 8,
    show_rownames = TRUE,
    show_colnames = TRUE,
    legend = FALSE
  )
  
}

# venn_plot <- function (i, CConsistency_marker) {
#   data <- CConsistency_marker[[i]]
#   venn::venn(data,
#        zcolor = c("#00ced1", "#F87168", "#fec600"), opacity = 0.5,
#        ilcs = 1.5, sncs = 2,
#        box = F, ggplot = T)
# }

venn_plot <- function (i, CConsistency_marker) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  data <- CConsistency_marker[[i]]
  return(
  VennDiagram::venn.diagram(
    data, filename = NULL, category.names = c("Set 1" , "Set 2" , "Set 3"), main = paste("Cluster", i),
    col = "black", fill = c("#00ced1", "#F87168", "#fec600")
  ))
}


volcanovalue <- function(subdata_cluster = subdata_cluster, control = "Biopsy", case = "PBMC") {
  # pvalue
  data1 <- data_ttest(subdata_cluster)
  selected_markers <- colnames(data1)[-c(which(colnames(data1) == "filename"), which(colnames(data1) == "condition"), which(colnames(data1) == "cluster"))]
  pvalue_res <- try(as.data.frame(-log10(sapply(selected_markers, pvalue_method, data = data1))), silent = T)
  if (class(pvalue_res) == "try-error") return(NULL)
  # FC
  data1$condition <- factor(data1$condition, levels = c(control, case))
  data2 <- dplyr::group_by(data1, condition) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
  data3 <- data2[, -c(which(colnames(data2) == "condition"), which(colnames(data2) == "cluster"))]
  FC_res <- t(log2(data3[2,]/data3[1,])) # case/control
  FC_res <- data.frame(FC_res[!is.infinite(FC_res),])
  # combine pvalue & FC
  res <- dplyr::inner_join(tibble::rownames_to_column(pvalue_res), tibble::rownames_to_column(FC_res), by = "rowname")
  res <- na.omit(res)
  colnames(res) <- c("marker", "pvalue", "log2FC")
  return(res)
}


volcano_plot <- function (i, volcanodata) {
  data <- volcanodata[[i]]
  # if (is.null(data)) return(NULL)
  
  threshold <- ifelse(data$pvalue > -log10(0.05) & abs(data$log2FC) >= 1,
                      ifelse(data$log2FC >= 1 , "Up", "Down"),
                      "Not")
  
  data <- cbind(data, threshold)
  ggplot2::ggplot(data, aes(x = log2FC, y = pvalue, colour = threshold)) +
    xlab("log2FC") + ylab("-log10pvalue") +
    labs(title = paste("Cluster", i)) +
    geom_point(size = 2, alpha = 1) +
    ylim(floor(range(data$pvalue)[1]), ceiling(range(data$pvalue)[2])) +
    xlim(floor(range(data$log2FC)[1]), ceiling(range(data$log2FC)[2])) +
    scale_color_manual(breaks = c("Down", "Not", "Up"), values = c("blue", "grey", "red"), name = NULL) +
    geom_vline(xintercept = c(-1, 1), lty = 2, colour = "#000000") +
    geom_hline(yintercept = c(-log10(0.05)), lty = 2, colour = "#000000") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
}


cluster_plot <- function (data_with_cluster) {
  data_dimred <- mergedimred(data_with_cluster)
  
  data_dimred$cluster <- as.factor(data_dimred$cluster)
  cluster_num <- length(unique(data_dimred$cluster))
  col_legend_row <- ceiling(cluster_num/15)
  
  edata <- data_dimred[, c("tsne_1", "tsne_2", "cluster")]
  colnames(edata) <- c("x", "y", "z")
  center <- aggregate(cbind(x, y) ~ z, data = edata, median)
  
  plot1 <- ggplot() +
    geom_point(data = data_dimred, aes(x = tsne_1, y = tsne_2, colour = cluster), size = 0.5) +
    geom_text(data = center, aes(x = x, y = y, label = z)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)))
  plot2 <- plot1 + facet_wrap(~ condition, nrow = 1) + ggtitle("Grid plot divided by condition")
  
  # return(ggpubr::ggarrange(plot1, plot2, ncol = 2, common.legend = TRUE, legend = "bottom", widths = c(1, 2)))
  return(ggpubr::ggarrange(plot1, plot2, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(2, 1)))
}


recall_plot <- function(data_with_cluster, known_marker) {
  data_known_marker <- data_with_cluster[colnames(data_with_cluster) %in% c(known_marker, "filename", "condition")]
  data_known_marker2 <- dplyr::group_by(data_known_marker, filename, condition) %>%
    dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
  data_known_marker3 <- reshape2::melt(data_known_marker2, id = c("filename", "condition"))
  colnames(data_known_marker3) <- c("filename", "condition", "protein", "expression")
  
  ggplot(data = data_known_marker3, mapping = aes(x = protein, y = expression, color = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, size = 1.1, na.rm = T) +
    scale_color_manual(values = c("#810061", "#fec600")) +
    ggpubr::stat_compare_means(aes(label = ..p.signif..), method = "t.test") +
    ggtitle("Differential expression of known markers between two cnditions") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(size = 10, angle = 10, vjust = 0.7),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10)
    )
}


densityplot <- function (data_with_condition, markername) {
  data_with_condition2 <- data_with_condition %>%
    dplyr::filter_at(vars("filename"),all_vars(.%in% unique(data_with_condition$filename))) %>%
    dplyr::select(one_of(c("filename",markername)))
  
  data_with_condition3 <- tidyr::gather_(data_with_condition2, key_col = "Markers", value = "exprs", markername)
  
  cols <- c(RColorBrewer::brewer.pal(12, "Paired"),
            RColorBrewer::brewer.pal(5, "Accent"),
            RColorBrewer::brewer.pal(3, "Dark2"),
            RColorBrewer::brewer.pal(8, "Pastel1"),
            RColorBrewer::brewer.pal(9, "Set1"),
            RColorBrewer::brewer.pal(12, "Set3"))
  
  p1 <- ggplot(data_with_condition3) +
    ggridges::geom_density_ridges(aes(x = exprs, y = filename, fill = filename)) +
    scale_fill_manual(values = cols[1:length(unique(data_with_condition3$filename))]) +
    # facet_wrap(as.formula("filename ~ Markers"), scales = "free_y", ncol = 1) +
    # facet_wrap(~ filename, scales = "free_y", ncol = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top")
  return(p1)
}

pointdensity_plot <- function(data, x_marker, y_marker) {
  
  data1 <- data[c(x_marker, y_marker)]
  colnames(data1) <- c("x", "y")
  
  p1 <- ggplot(data = data1, aes(x = x, y = y)) +
    ggpointdensity::geom_pointdensity() +
    scale_colour_viridis_c() +
    theme_classic() +
    xlab(x_marker) + ylab(y_marker) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 10))
  return(p1)
}

merge_multi <- function(multi_flowFrame,condition_info) {
  
  exprsL <- list()
  label <- list()
  condition <- list()
  
  #
  multi_flowFrame <- Filter(function(x) length(x@exprs) > 0, multi_flowFrame)
  
  for (i in 1:length(multi_flowFrame)) {
    
    fcsFile <- multi_flowFrame[i]
    index <- sapply(strsplit(rownames(fcsFile[[1]]@exprs), split = "_"), function(x) x[length(x)])
    condition[[i]] <- rep(condition_info[i], length(index))
    
    
    name <- names(fcsFile)
    fcs <- fcsFile[[1]]
    
    pd <- fcs@parameters@data # abstract the protein name and its characterization
    
    size_channels <- grep("FSC|SSC", colnames(fcs@exprs), ignore.case = TRUE)
    
    marker_id <- seq_along(colnames(fcs@exprs))
    
    exprs <- fcs@exprs[, marker_id, drop = FALSE] # data after removed the "Time" and "Event"
    
    if (length(size_channels) > 0) { # judge whether it is FC data
      if (any(size_channels %in% marker_id)) {
        used_size_channel <- size_channels[size_channels %in% marker_id]
        used_size_channel_id <- match(used_size_channel, marker_id)
        exprs[, used_size_channel_id] <- apply(exprs[, used_size_channel_id, drop = FALSE], 2,
                                               function(x) scaleData(x, range = c(0,4.5))) # scale the "FSC" and "SSC" data
      }
    }
    col_names <- paste0(pd$desc, "(", pd$name, ")")
    colnames(exprs) <- col_names[marker_id]
    row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
    exprsL[[i]] <- exprs
  }
  
  merged_data <- do.call(rbind, exprsL)
  merged_lable <- unlist(label)
  merged_condition <- unlist(condition)
  
  return(list(data = merged_data,
              condition = merged_condition))
}

get_markers <- function(marker_path) {
  
  # check 'marker_path'
  if(grepl(".csv$", marker_path)){
    markers <- read.csv(file = marker_path, header = T)
  } else if(grepl(".xlsx$", marker_path)){
    markers <- openxlsx::read.xlsx(marker_path)
  } else if(sum(dim(marker_path))==0){
    stop("No CSV or XLSX file found, please check the parameter of 'marker_path'.")
  }
  
  markers$posGene <- gsub(" ", "", markers$posGene)
  markers$negGene <- gsub(" ", "", markers$negGene)
  
  # correct gene symbols from the given    (up-genes)
  markers$posGene <- sapply(1:nrow(markers), function(i){
    markers_all <- gsub(" ", "", unlist(strsplit(markers$posGene[i], ",")))
    markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all <- sort(markers_all)
    
    return(paste0(markers_all, collapse=","))
  })
  
  # correct gene symbols from the given DB (down-genes)
  markers$negGene <- sapply(1:nrow(markers), function(i){
    markers_all <- gsub(" ", "", unlist(strsplit(markers$negGene[i],",")))
    markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all <- sort(markers_all)
    
    return(paste0(markers_all, collapse = ","))
    
  })
  
  markers$posGene <- gsub("///",",",markers$posGene)
  markers$posGene <- gsub(" ","",markers$posGene)
  markers$negGene <- gsub("///",",",markers$negGene)
  markers$negGene <- gsub(" ","",markers$negGene)
  
  gs <- lapply(1:nrow(markers), function(j) gsub(" ","",unlist(strsplit(toString(markers$posGene[j]),","))))
  names(gs) <- markers$cellName
  gs2 <- lapply(1:nrow(markers), function(j) gsub(" ","",unlist(strsplit(toString(markers$negGene[j]),","))))
  names(gs2) <- markers$cellName
  
  return(list(gs_positive = gs, gs_negative = gs2))
}

cell_annotation <- function(data, marker_path, colsToUse = NULL,
                            resolution = "cell") {
  
  if (resolution == "cluster") {
    set.seed(123)
    data_cluster <- suppressWarnings(FlowSOM::FlowSOM(as.matrix(data), nClus = 50,  colsToUse = colsToUse))
    node_lable <- data_cluster[["map"]][["mapping"]][,1]
    cluster_label <- data_cluster[["metaclustering"]][node_lable]
    if (any(table(cluster_label) == 0)) {
      cluster_label <- droplevels(cluster_label)
    }

  }
  
  colnames(data) <- gsub(pattern = "\\(.*\\)", "", x = colnames(data))
  data <- t(data)
  
  markers <- get_markers(marker_path = marker_path)
  
  # assign cell types
  es.max <- anno_score(input_data = data, scaled = TRUE, gs = markers$gs_positive, gs2 = markers$gs_negative)
  
  # predicted label
  if (resolution == "cell") {
    pred_label <- apply(es.max, 2, function(col) {
      rownames(es.max)[which.max(col)]
    })
  } else if (resolution == "cluster") {
    # Extract top cell types for each cluster
    cL_resutls <- do.call("rbind", lapply(unique(cluster_label), function(cl){
      if (is.array(es.max[, cluster_label == cl])) {
        es.max.cl <- sort(rowSums(es.max[, cluster_label == cl]), decreasing = T)
      } else {
        es.max.cl <- sort(es.max[, cluster_label == cl], decreasing = T)
      }
      
      return(head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, 
                             ncells = sum(cluster_label == cl)), 10))
    }))
    anno_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    anno_scores$type[as.numeric(as.character(anno_scores$scores)) < anno_scores$ncells/4] <- "Unknown"
    anno_scores <- anno_scores[order(anno_scores$cluster),]
    pred_label <- anno_scores$type[cluster_label]
  }
  return(pred_label)
}

anno_score <- function(input_data, scaled = T, gs, gs2 = NULL, to_uppercase = T){
  
  # check input matrix
  if(!is.matrix(input_data)){
    warning("input_data doesn't seem to be a matrix")
  } else {
    if(sum(dim(input_data))==0){
      warning("The dimension of input input_data matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat <- sort(table(unlist(gs)), decreasing = T)
  marker_sensitivity <- data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), 
                                                                              to = c(0,1), 
                                                                              from = c(length(gs),1)),
                                   gene_ = names(marker_stat), stringsAsFactors = F)
  
  # convert gene names to Uppercase
  if(to_uppercase){
    rownames(input_data) = toupper(rownames(input_data));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){
    GeneIndToKeep = rownames(input_data) %in% as.character(gs[[d_]])
    rownames(input_data)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(input_data) %in% as.character(gs2[[d_]]); rownames(input_data)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(scaled) Z <- t(scale(t(input_data))) else Z <- input_data
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  return(es.max)
}

