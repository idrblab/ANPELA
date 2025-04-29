dr_multi <- function(TIres, D){
  colors <- colorRampPalette(c("#eef4ed","#97c8c5","#4661a5","#183f7f"))(length(unique(D$timepoint)))
  point_colors <- colors[as.factor(D$timepoint)]
  if (!is.null(TIres$linInd)) {
    current_curve <- TIres$linInd
    sds <- slingshot::SlingshotDataSet(TIres$crv1)
    curve_weights <- sds@curves[[current_curve]]$w
    curve_points <- which(curve_weights > 0)
    point_colors[-curve_points] <- "grey80"
  }
  dr.plot <- plot(TIres$dimRed, col = point_colors,
       pch=16, cex = 1.2,
       asp = 1,axes = FALSE,xlab = "", ylab = "")
  if (!is.null(TIres$linInd)) {
    lines(slingshot::SlingshotDataSet(TIres$crv1), linInd=TIres$linInd, lwd=6, col="grey40")
  } else lines(slingshot::SlingshotDataSet(TIres$crv1), lwd=6, col="grey40")
  return(dr.plot)
}

dr <- function(TIres, D, sub = F, seed = 3, cell.subset = 0.2){

  # load the data
  input_matrix <- cbind(D$expr, D$timepoint, TIres$pseudotime, TIres$trajectory)
  colnames(input_matrix)[(length(input_matrix)-3):(length(input_matrix)-2)] <- c("D.timepoint", "pseudot")

  # load the dimensionality reduction method results
  dr_matrix <- TIres$dimRed
  colnames(dr_matrix) <- c("V1","V2")

  # calculate plot parameters
  npar <- which(colnames(input_matrix) == "D.timepoint") - 1
  ntimepoints <- length(unique(input_matrix$D.timepoint))

  # extract data to plot
  # dat <- input_matrix[, 1:npar]
  Timepoints <- input_matrix[, npar + 1]
  Pseudotime <- input_matrix[, npar + 2]
  trajectory <- input_matrix[, c(npar + 3, npar + 4)]


  if (sub) {
    set.seed(seed)
    ncells <- nrow(input_matrix)
    nsub.cells <- round(cell.subset*ncells)
    subIndex <- sample(1:ncells, nsub.cells)
    subsample <- 1:ncells %in% subIndex
    dr_matrix <- cbind(dr_matrix, drsubsample = subsample)
    Timepoints <- cbind(Timepoints, tmsubsample = subsample)
  } else {
    ncells <- nrow(input_matrix)
    subsample <- rep(1,ncells)
    dr_matrix <- cbind(dr_matrix, drsubsample = subsample)
    Timepoints <- cbind(Timepoints, tmsubsample = subsample)
  }

  # sort the plot values in ascending pseudotime order
  ix <- order(Pseudotime)
  # dat <- dat[ix,]
  Timepoints <- Timepoints[ix,]
  # Pseudotime <- Pseudotime[ix]
  # if (sum(is.na(Pseudotime)) > 0){
  #   print("The algorithm returned a branching trajectory. Cannot visualize!")
  #   next
  # }

  trajectory <- trajectory[ix,]
  colnames(trajectory) <- c("V1","V2")

  dr_matrix <- dr_matrix[ix,]

  cols <- c(RColorBrewer::brewer.pal(12, "Paired"),
            RColorBrewer::brewer.pal(5, "Accent"),
            RColorBrewer::brewer.pal(3, "Dark2"),
            RColorBrewer::brewer.pal(8, "Pastel1"),
            RColorBrewer::brewer.pal(9, "Set1"),
            RColorBrewer::brewer.pal(12, "Set3"))

  # create the reduced manifold figure
  dr.df <- subset(data.frame(dr_matrix, Timepoints), drsubsample == 1)

  ratio.display <- 1/1
  ratio.values <- (max(dr.df$V1)-min(dr.df$V1))/(max(dr.df$V2)-min(dr.df$V2))


  dr.plot <- ggplot() +
    geom_point(data = dr.df, aes(x = V1,y = V2, color = as.factor(Timepoints)), alpha = 0.5, size = 1) +
    scale_color_manual(values = cols, name = "time") +
    geom_point(data = trajectory, aes(x = V1, y = V2), size = 0.5, color = "#666666") +
    xlab("reduced dimension 1") + ylab("reduced dimension 2") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 10),
          legend.position = "bottom"
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1, title = "experimental time")) +
    coord_fixed(ratio.values / ratio.display)


  return(dr.plot)
}



abund_pt_single <- function(to_plot, prot, dat) {
  cols <- c(RColorBrewer::brewer.pal(12, "Paired"),
            RColorBrewer::brewer.pal(5, "Accent"),
            RColorBrewer::brewer.pal(3, "Dark2"),
            RColorBrewer::brewer.pal(8, "Pastel1"),
            RColorBrewer::brewer.pal(9, "Set1"),
            RColorBrewer::brewer.pal(12, "Set3"))
  y_min <- min(dat)
  y_max <- max(dat)

  p4 <- ggplot(to_plot[[prot]], aes(pst, y)) +
    geom_point(aes(color = t), alpha = 0.5, size = 0.8) +
    scale_color_manual(values = cols, name = "time") +
    ylim(y_min, y_max) + xlim(0, 1) +
    facet_wrap(~t, ncol = 4) +
    ggtitle(prot) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          legend.title = element_text(size = 10),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1, title = "experimental time"))

  # create the final subplot figure
  # library(cowplot)
  y_label <- cowplot::ggdraw() + cowplot::draw_label("protein abundance", size = 10, angle = 90, hjust = 0.5)
  x_label <- cowplot::ggdraw() + cowplot::draw_label("standardized pseudotime", size = 10, hjust = 0.5)
  p5 <- cowplot::plot_grid(y_label, p4, rel_widths = c(0.1, 1))
  p6 <- cowplot::plot_grid(p5, x_label, rel_heights = c(1, 0.05), nrow = 2)

  # add the legend
  leg <- cowplot::get_legend(p4 + theme(legend.position = 'bottom'))
  p7 <- cowplot::plot_grid(leg, p6, ncol = 1, rel_heights = c(0.05, 1), axis = "l")

  return(p7)
}


robustness_new_plot <- function(input_matrix, finalMatrix, method = method) {
  plist <- list()
  time_point <- initial_order <- subset_order <- as.data.frame(matrix(nrow = nrow(finalMatrix[[1]]$pseudotime), ncol = length(finalMatrix)))
  
  # # Compare the initial ordering with those returned by the subsets
  for (n in 1:length(finalMatrix)){
    if (finalMatrix[[n]]$lineages == 1){ # "Linear"
      subset_order[,n] <- finalMatrix[[n]]$pseudotime
    } else {
      subset_order[,n] <- apply(finalMatrix[[n]]$pseudotime,1, function(row) mean(row,na.rm = T))
    }
    time_point[,n] <- input_matrix$D.timepoint[finalMatrix[[n]]$subIndex]
    initial_order[,n] <- input_matrix$pseudot[finalMatrix[[n]]$subIndex]
    
    
    df <- data.frame(a = initial_order[,n], b = subset_order[,n], c = time_point[,n])
    df <- na.omit(df)
    df$a <-  (df$a - min(df$a)) / (max(df$a) - min(df$a))
    df$b <-  (df$b - min(df$b)) / (max(df$b) - min(df$b))
    
    initial_ix <- order(df$a)
    df <- df[initial_ix,]
    df$a <- shift_start(exp.time = df$c, pseudotime = df$a, ps.ordered = T, circular = F, randomize_t0 = T)
    
    subset_ix <- order(df$b)
    df <- df[subset_ix,]
    df$b <- shift_start(exp.time = df$c, pseudotime = df$b, ps.ordered = T, circular = F, randomize_t0 = T)
    
    if (method == "Spearman rank correlation") {
      metric_value <- paste("Spearman rank correlation:", round(cor(df$a, df$b, use = "pairwise.complete.obs", method = "spearman"), 2))
    }
    if (method == "Kendall rank correlation") {
      metric_value <- paste("Kendall rank correlation:", round(cor(df$a, df$b, use = "pairwise.complete.obs", method = "kendall"), 2))
    }
    
    cols <- c(RColorBrewer::brewer.pal(12, "Paired"), 
              RColorBrewer::brewer.pal(5, "Accent"), 
              RColorBrewer::brewer.pal(3, "Dark2"),
              RColorBrewer::brewer.pal(8, "Pastel1"),
              RColorBrewer::brewer.pal(9, "Set1"), 
              RColorBrewer::brewer.pal(12, "Set3"))
    
    
    plist[[n]] <- ggplot(df) + 
      geom_point(size = 1, aes(x = a, y = b, color = as.factor(c))) +
      scale_color_manual(values = cols, name = "time") +
      geom_smooth(aes(x = a, y = b), method = "lm") + 
      ylim(0, 1) + xlim(0, 1) + 
      xlab("initial pseudotime") + ylab("pseudotime of subset") +
      # annotate("text", x = 0.3, y = 0.7, label = metric_value, size = 5) +
      ggtitle(paste("Subset", n, "  ", metric_value)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 10)) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1, title = "experimental time")) +
      coord_fixed(1)
    
    
    
  }
  
  plist1 <- ggpubr::ggarrange(plotlist = plist, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
  return(plist1)
}


roughness_plot <- function(TIres, D) {
  # load the data
  input_matrix <- cbind(D$expr, D$timepoint, TIres$pseudotime)
  colnames(input_matrix)[(length(input_matrix)-1):(length(input_matrix))] <- c("D.timepoint", "pseudot")

  # calculate plot parameters
  npar <- which(colnames(input_matrix) == "D.timepoint") - 1 # 14

  # extract data to plot
  dat <- input_matrix[, 1:npar]
  Pseudotime <- input_matrix[, npar + 2]

  # # Find axis limits => you need common limits of asinh or log for y axis
  # y_min <- min(dat)
  # y_max <- max(dat)

  # sort the plot values in ascending pseudotime order
  ix <- order(Pseudotime)
  dat0 <- dat[ix,]
  set.seed(123)
  dat1 <- data.frame(t(sample(data.frame(t(dat0)))))

  test0 <- abs(sapply(dat0, diff))
  test1 <- abs(sapply(dat1, diff))

  dfm0 <- reshape2::melt(test0)
  dfm1 <- reshape2::melt(test1)
  final <- dplyr::full_join(dfm0, dfm1, by = c("Var1", "Var2"))
  final <- reshape2::melt(final, id.vars = c("Var1", "Var2"))

  cols <- c(RColorBrewer::brewer.pal(12, "Paired"),
            RColorBrewer::brewer.pal(5, "Accent"),
            RColorBrewer::brewer.pal(3, "Dark2"),
            RColorBrewer::brewer.pal(8, "Pastel1"),
            RColorBrewer::brewer.pal(9, "Set1"),
            RColorBrewer::brewer.pal(12, "Set3"))



  param.plot <- list()
  # repeat for each cell parameter
  for (j in 1:npar){
    data <- subset(final, Var2 == names(dat)[j])

    plot <- ggplot(data, aes(x = variable, y = value, color = variable)) + geom_boxplot()
    outlier_data <- layer_data(plot)['outliers']

    for(i in levels(data$variable)){
      data[data$variable == i, "value"] <- replace(data[data$variable == i, "value"],
                                                   data[data$variable == i, "value"] %in% outlier_data[which(levels(data$variable) %in% i), ][[1]],
                                                   NA)
    }

    mediandata0 <- dplyr::group_by(data, variable) %>%
      dplyr::summarise(median(value, na.rm = T))
    mediandata <- as.numeric(as.matrix(mediandata0[,2]))


    param.plot[[j]] <- ggplot(data, aes(x = variable, y = value, color = variable)) +
      geom_violin() +
      geom_boxplot(outlier.shape = NA, width = 0.1, na.rm = T) +
      geom_hline(yintercept = mediandata, colour = c("#810061", "#fec600"), linetype = "dashed") +
      scale_color_manual(values = c("#810061", "#fec600")) +
      scale_x_discrete(breaks = c("value.x", "value.y"), label = c("TI", "naive")) +
      ggsignif::geom_signif(comparisons = list(c("value.x", "value.y")), map_signif_level = F,
                            test = t.test, test.args = list(alternative = "less", paired = F),
                            vjust = 1.5, textsize = 4) +
      ggtitle(colnames(dat)[j]) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            # axis.text.x = element_text(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 10),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 10))
  }

  # draw the parameters
  p0 <- cowplot::plot_grid(plotlist = param.plot, ncol = 4)

  return(p0)
}


abund_pt_plot <- function(TIres, D) {
  # load the data
  input_matrix <- cbind(D$expr, D$timepoint, TIres$pseudotime)
  colnames(input_matrix)[(length(input_matrix)-1):(length(input_matrix))] <- c("D.timepoint", "pseudot")

  # calculate plot parameters
  npar <- which(colnames(input_matrix) == "D.timepoint") - 1 # 14
  ntimepoints <- length(unique(input_matrix$D.timepoint)) # 8

  # extract data to plot
  dat <- input_matrix[, 1:npar]
  Timepoints <- input_matrix[, npar + 1]
  Pseudotime <- input_matrix[, npar + 2]
  valid_idx <- !is.na(Pseudotime)
  Pseudotime <- Pseudotime[valid_idx]
  Timepoints <- Timepoints[valid_idx]
  dat <- dat[valid_idx, ]

  # Find axis limits => you need common limits of asinh or log for y axis
  y_min <- min(dat)
  y_max <- max(dat)

  # sort the plot values in ascending pseudotime order
  # ix <- order(Pseudotime)
  # dat <- dat[ix,]
  # Timepoints <- Timepoints[ix]
  # Pseudotime <- Pseudotime[ix]
  # if (sum(is.na(Pseudotime)) > 0){
  #   print("The algorithm returned a branching trajectory. Cannot visualize!")
  #   return(NULL)
  # }
  # Pseudotime <- na.omit(Pseudotime)

  df <- data.frame(a = Pseudotime, b = Timepoints, c = 1:length(Pseudotime))
  initial_ix <- order(df$a)
  df <- df[initial_ix,]
  df$a <- shift_start(exp.time = df$b, pseudotime = df$a, ps.ordered = T, circular = F, randomize_t0 = T)
  initial_ix2 <- order(df$c)
  df <- df[initial_ix2,]
  Pseudotime <- df$a

  # rescale pseudotime to [0,1] for the figures between algorithms to be comparable
  Pseudotime_zscore <- (Pseudotime - min(Pseudotime, na.rm = T))/(max(Pseudotime, na.rm = T) - min(Pseudotime, na.rm = T))

  # reposition the pseudotime values so that cells in experimental time 0 are close to the origin
  # Pseudotime_zscore <- shift_start(exp.time = Timepoints, pseudotime = Pseudotime_zscore, ps.ordered = T, circular = T, randomize_t0 = T)

  # add the spline fit
  spline.df <- calc_spline(dat, Pseudotime_zscore)
  dfm <- reshape2::melt(spline.df, id.vars = 'time')

  cols <- c(RColorBrewer::brewer.pal(12, "Paired"),
            RColorBrewer::brewer.pal(5, "Accent"),
            RColorBrewer::brewer.pal(3, "Dark2"),
            RColorBrewer::brewer.pal(8, "Pastel1"),
            RColorBrewer::brewer.pal(9, "Set1"),
            RColorBrewer::brewer.pal(12, "Set3"))


  param.plot <- list()
  to_plot <- list()

  # repeat for each cell parameter
  for (j in 1:npar){
    plotvar <- dat[,j]
    to_plot[[j]] <- data.frame(pst = Pseudotime_zscore, y = plotvar, t = as.factor(Timepoints))
    names(to_plot)[[j]] <- colnames(dat)[j]
    # make the plot
    param.plot[[j]] <- ggplot(to_plot[[j]], aes(pst, y)) +
      geom_point(aes(color = t), alpha = 0.5, size = 0.5) +
      scale_color_manual(values = cols, name = "time") +
      ylim(y_min,y_max) + xlim(0,1) +
      geom_line(data = subset(dfm, variable == names(dat)[j]), aes(x = time, y = value), size = 1) +
      ggtitle(colnames(dat)[j]) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 10),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 10)) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1, title = "experimental time"))
  }

  # draw the parameters
  p0 <- cowplot::plot_grid(plotlist = param.plot, ncol = 4)

  # create the final subplot figure
  # library(cowplot)
  y_label <- cowplot::ggdraw() + cowplot::draw_label("protein abundance", size = 10, angle = 90, hjust = 0.5)
  x_label <- cowplot::ggdraw() + cowplot::draw_label("standardized pseudotime", size = 10, hjust = 0.5)
  p1 <- cowplot::plot_grid(y_label, p0, rel_widths = c(0.05, 1))
  p2 <- cowplot::plot_grid(p1, x_label, rel_heights = c(1, 0.05), nrow = 2)

  # add the legend
  leg <- cowplot::get_legend(param.plot[[1]] + theme(legend.position = 'bottom'))
  p3 <- cowplot::plot_grid(p2, leg, ncol = 1, rel_heights = c(1,0.05), axis = "l")

  return(list(p3 = p3, spline.df = spline.df, to_plot = to_plot, dat = dat))
}


diagonal_plot_pre <- function(data, pro_seq, i) {
  ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    geom_segment(aes(xend = c(tail(x, n = -1), NA),
                     yend = c(tail(y, n = -1), NA)),
                 arrow = arrow(length = unit(0.4, "cm"),
                               type = "closed")) +
    scale_y_continuous(breaks = seq(1, length(pro_seq), 1),
                       label = pro_seq) +
    scale_x_continuous(breaks = seq(1, nrow(data), 1),
                       label = data$namex) +
    xlab("Order of the respective proteins in known signal transduction cascades") +
    ylab("Pseudo-temporal ordering of proteins at peak activation") +
    ggtitle(paste("Pathway", i)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.4))
}


diagonal_plot <- function(data, pro_seq) {
  p0 <- list()
  for (i in 1:length(data)) {
    p0[[i]] <- ggplot(data[[i]], aes(x = x, y = y)) +
      geom_point() +
      geom_segment(aes(xend = c(tail(x, n = -1), NA),
                       yend = c(tail(y, n = -1), NA)),
                   arrow = arrow(length = unit(0.4, "cm"),
                                 type = "closed")) +
      scale_y_continuous(breaks = seq(1, length(pro_seq), 1),
                         label = pro_seq) +
      scale_x_continuous(breaks = seq(1, nrow(data[[i]]), 1),
                         label = data[[i]]$namex) +
      xlab("Order of the respective proteins in known signal transduction cascades") +
      ylab("Pseudo-temporal ordering of proteins at peak activation") +
      ggtitle(paste("Pathway", i)) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, hjust = 0.4))
  }
  p1 <- cowplot::plot_grid(plotlist = p0, ncol = 2)
  return(p1)
}

arrow_plot <- function(data, pro_seq) {
  p0 <- list()
  for (i in 1:length(data)) {
    p0[[i]] <- ggplot(data[[i]], aes(x = x, y = y)) +
      geom_point() +
      geom_segment(aes(x=0, y=0, xend = c(tail(x, n = -1), NA), yend = c(tail(y, n = -1), NA)),
                   arrow = arrow(length = unit(0.4, "cm"),
                                 type = "closed")) +
      scale_y_continuous(breaks = c(1, 0, -1),
                         label = c("Consistent", as.character(data[[i]]$namex[1]), "Inconsistent")) +
      scale_x_continuous(breaks = seq(1, nrow(data[[i]]), 1),
                         label = data[[i]]$namex) +
      xlab("Order of the respective proteins in known signal transduction cascades") +
      ylab(paste("Consistency of protein activation after", as.character(data[[i]]$namex[1]))) +
      ggtitle(paste("Pathway", i)) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, hjust = 0.4))
  }
  p1 <- cowplot::plot_grid(plotlist = p0, ncol = 2)
  return(p1)
}


