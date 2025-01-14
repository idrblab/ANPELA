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
  if (class(data_dimred) == "try-error") return(NULL)
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
  data_heatmap <- data_heatmap[, -1]
  
  metabolomics::HeatMap(data_heatmap, colramp = gplots::colorpanel(100, "#c10022", "#ffffbd", "#2c3298"), dendrogram = "both",
                        margins = c(10, 10), key = FALSE, cexRow = 0.8, main = paste("cluster", i))
  
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
  if (is.null(data)) return(NULL)
  grid::grid.draw(VennDiagram::venn.diagram(
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
  ggplot(data, aes(x = log2FC, y = pvalue, colour = threshold)) +
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
  return(ggpubr::ggarrange(plot1, plot2, ncol = 1, common.legend = TRUE, legend = "bottom", widths = c(2, 1)))
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
