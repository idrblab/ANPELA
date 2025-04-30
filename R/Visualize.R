#' @title Visualize
#' @description Visualize() loads processed SCP data and generates a suite of visualizations designed to evaluate the quality and characteristics of either CSI or PTI results, sav-ing these plots for review.
#' @param studytype Character, the type of your study, including ‘CSI (Cell Subpopulation Identification)’ and ‘PTI (Pseudotime Trajectory Inference)’.
#' @param respath Character, the absolute path of the folder storing the resulting ‘info_saved.RData’ file and the ‘process_res’ folder of the function ‘Process’, ‘FCprocess’ or ‘MCprocess’ when the ‘save_processed_res’ parameter in these functions is set to ‘one_folder’.
#' @param save_processed_res Character, the format of the data processing output files. ‘one_folder’ denotes that successfully processed results will be saved as separate RData files in the ‘process_res’ folder. ‘one_RData’ denotes that all processed results will be saved as one RData file in the ‘process_res’ folder.
#' @param workflow Character, the combinations of data processing methods specified by users according to their research interests.
#' @param savepath Character, the absolute path of the folder which will store the results.
#' @param studyname Character, the filename of the PNG files generated that includes the value of this required character string parameter, serving as a user-specified identifier for the analysis run.
#' @param plot_metric Character, the metric(s) which will be visualized and saved.
#' @param clusteringM Character, the method of clustering the processed data prior to performance assessment, including ‘FlowSOM’ and ‘PhenoGraph’.
#' @param ncluster Integer, the number of clusters for meta clustering in FlowSOM.
#' @param ntop Integer, the number of the most differentially expressed markers that are truncated for calculating the CWrel value.
#' @param DEP Character, the absolute file path of the CSV file including the differentially expressed proteins used as the prior knowledge for the fourth criterion.
#'   <br>It is a table of one column without the column name, each table cell includes one protein typically in the format of "channel description (channel name)", for example: "CD20(FITC.A)".
#' @param TIM Character, the method of trajectory inference for the processed data prior to performance assessment, consisted of tra-jectory reconstruction and data space representation, including ‘scorpius_distSpear’, ‘scorpius_distPear’, ‘scorpius_distEucl’, ‘scorpius_distManh’, ‘slingshot_tSNE’, ‘slingshot_FLOWMAP’, ‘slingshot_PCA’, ‘slingshot_diffMaps’.
#' @param clustering.var Character, the vector naming channels to be used to calculate distances/differences between cells for clustering (if re-quested) and edge-drawing steps.
#' @param pathwayhierarchy Character, the absolute file path of the pathway hierarchy file.
#' @param Cc_metric Character, the assessing metric under Criterion Cc for the ‘PTI’ study type, including ‘Spearman correlation’ and ‘Kendall Rank Correlation’.
#' @return PNG files visualizing the assessment standards.
#' @export
#'
#' @examples
#' \donttest{
#' }

Visualize <- function(
    studytype = c("CSI","PTI"),
    respath,  save_processed_res ="one_folder", workflow, savepath = "./ANPELA_res", studyname = NULL,
    plot_metric = c("Ca_metric","Cb_metric","Cc_metric","Cd_metric"),
    #CSI
    clusteringM = c("FlowSOM"), ncluster = 8,
    ntop = NULL, DEP = NULL,
    #PTI
    TIM = "scorpius_distSpear",
    clustering.var = NULL, pathwayhierarchy = NULL,
    Cc_metric = "Spearman rank correlation"
){

  #load data&metadata
  datapath <- list.files(paste0(respath, "/process_res/"), pattern = "\\.RData$", full.names = T)
  if (length(datapath)==0){
    stop("The parameter of 'respath' is incorrect. The 'datafile' cannot be loaded.")
  }
  if (save_processed_res == "one_RData") {
    assign("data",load(datapath))
    data <- get(data)
    metadata <- data[["metadata"]]
    colsToUse <- data[["index_TIclass"]]
  } else if (save_processed_res == "one_folder") {
    info_saved <- try(load(paste0(respath, "/info_saved.RData")), silent = T)
    if (class(info_saved) == "try-error") {
      stop("The parameter of 'respath' is incorrect. The 'info_saved.RData' cannot be loaded.")
    } else if (info_saved == "info_saved"){
      load(paste0(respath, "/info_saved.RData"))
    }
    metadata <- info_saved[["metadata"]]
    colsToUse <- info_saved[["index_TIclass"]]
  }
  condition_info  <- metadata[["condition"]]


  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }

  errorplot <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", x = 0.5, y = 0.5,
                      label = paste("During the plotting process,\n",
                                    "the program had some problems and will not be able to display the output."),
                      cex = 5, color = "black") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
  ################CSI Visualize#####################
  if (studytype == "CSI") {
    try(source("./CSI/1readfcs.R"))
    try(source("./CSI/2cluster.R"))
    try(source("./CSI/3criteria.R"))
    try(source("./CSI/4plot.R"))
    #arg_check
    # clusteringM
    if (missing(clusteringM)) {
      clusteringM <- "FlowSOM"
    } else {
      clusteringM <- match.arg(clusteringM)
    }


    # ncluster
    if (clusteringM == "FlowSOM") {
      if (ncluster < 2 || ncluster %% 1 != 0) {
        stop("The number of clusters for meta clustering in FlowSOM should be positive integer number and greater than 1.")
      }
    }

    # ntop
    if (missing(ntop)) {
      ntop <- floor(length(colsToUse)/2)
    } else if (ntop >= length(colsToUse) || ntop %% 1 != 0 || ntop <= 0) {
      stop("The value of 'ntop' is incorrect. It should be positive whole number and less than the number of your selected markers.")
    }



    # known_marker
    if (is.null(DEP)) {
      cat("*************************************************************************", "\n")
      cat("The standardized marker names are listed below. \n")
      cat("*************************************************************************", "\n")
      if (save_processed_res == "one_RData") {
        cat(paste0(colnames(data$AP2_pro1_frame_classTI[vapply(data$AP2_pro1_frame_classTI, Negate(is.null), NA)][[1]][[1]]@exprs), collapse = "\n"), "\n")
      } else if (save_processed_res == "one_folder") {
        load(datapath[1])
        cat(paste0(colnames(res[[1]]@exprs), collapse = "\n"), "\n")
        rm(res)
      }
      cat("*************************************************************************", "\n")
      cat("Please enter the marker names which are separated by comma on a single line.
        \nFor example, CD103(La139Di), CD11b(Nd144Di), CD8a(Nd146Di), CD7(Sm147Di)")
      DEP <- readline("Now, you can select the known biomarker(s) differentially expressed between two conditions (if NULL, please press the Enter key directly):")
    } else if (file.exists(DEP)) {
      DEP_data <- read.csv(DEP, header = F, stringsAsFactors = FALSE)
      DEP <-paste(as.character(unlist(DEP_data)), collapse = ",")
    } else if (!DEP == "" & !file.exists(DEP)){
      stop("The filepath of parameter 'DEP' is incorrect. Please input the absolute filepath of the csv file.")
    }

    for ( w in 1:length(workflow)){
      dataset_name <- workflow[w]
      print(paste0("The visualization of data processed by workflow '", dataset_name, "' begins."))
      if (save_processed_res == "one_folder"){
        load(grep(dataset_name, datapath, value = T))
      } else if (save_processed_res == "one_RData") {
        res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      }
      if (class(res) == "try-error") {
        print(paste0("Can not load data for dataset: ", dataset_name))
        next
      }

      AP2_processed_D_class <- res
      rm(res)


      # AP2_processed_data_class0
      data1 <- try(as.data.frame(readfcs_multi(fcsFiles = AP2_processed_D_class, mergeMethod = "all")), silent = T)
      if (class(data1) == "try-error") {
        print(paste0("Can not read data for dataset: ", dataset_name))
        next
      }
      colnames(data1) <- sub("<", "_", colnames(data1))
      colnames(data1) <- sub(">", "_", colnames(data1))
      rm(AP2_processed_D_class)

      condition <- metadata
      data_with_filename <- cbind(data1, "filename" = sub("_[0-9]*$","", row.names(data1)))
      condition_label <- as.character(condition[match(data_with_filename[,"filename"], condition$filename),2])
      data_with_condition <- cbind(data_with_filename, "condition" = condition_label)
      AP2_processed_data_class0 <- data_with_condition
      rm(data1, condition, data_with_filename, condition_label, data_with_condition)


      # AP2_processed_data_class
      AP2_processed_data_class <- AP2_processed_data_class0[, c(colsToUse, "filename", "condition")]
      rm(AP2_processed_data_class0)


      # cluster_label
      cluster_label <- try(data_cluster(data = AP2_processed_data_class[,1:(dim(AP2_processed_data_class)[2]-2)],
                                        method = clusteringM, Phenograph_k = Phenograph_k, FlowSOM_k = ncluster,
                                        FlowSeed = 40), silent = T)

      if (class(cluster_label) == "try-error") {
        print(paste0("Can not cluster data for dataset: ", dataset_name))
        next
      }


      # data_with_cluster
      data_with_cluster <- cbind(AP2_processed_data_class, "cluster" = cluster_label)


      # test_KNN
      subdata_cluster <- list()
      subdata_cluster_DEG <- list()
      train <- list()
      test <- list()
      KNN_res <- list()

      for (j in 1:length(unique(cluster_label))) {

        subdata_cluster[[j]] <- subset(as.data.frame(data_with_cluster, stringsAsFactors = F), cluster == j)
        subdata_cluster[[j]]$condition <- as.factor(subdata_cluster[[j]]$condition)

        non_DEG <- try(feature_selection(subdata_cluster[[j]]), silent = T)
        if ((class(non_DEG) == "try-error") || (length(non_DEG) == length(subdata_cluster[[j]]) - 3)) next # The second case is that the number of non_DEG is equal to the number of markers in the original data, that is, subdata_cluster_DEG will have no data

        subdata_cluster_DEG[[j]] <- select(subdata_cluster[[j]], -one_of(non_DEG))

        data_condition1 <- subset(subdata_cluster_DEG[[j]], condition == unique(subdata_cluster_DEG[[j]]$condition)[1])
        data_condition2 <- subset(subdata_cluster_DEG[[j]], condition == unique(subdata_cluster_DEG[[j]]$condition)[2])

        set.seed(123)
        tr_ind1 <- sample(1:nrow(data_condition1), nrow(data_condition1)*0.7)
        tr_ind2 <- sample(1:nrow(data_condition2), nrow(data_condition2)*0.7)

        train[[j]] <- rbind(data_condition1[tr_ind1,], data_condition2[tr_ind2,])
        test[[j]] <- rbind(data_condition1[-tr_ind1,], data_condition2[-tr_ind2,])


        KNN_res[[j]] <- try(FNN::knn(train = train[[j]][, 1:(ncol(train[[j]]) - 3), drop = FALSE],
                                     test = test[[j]][, 1:(ncol(test[[j]]) - 3), drop = FALSE],
                                     cl = train[[j]]$condition, k = 5, prob = TRUE), silent = TRUE)
      }

      test_KNN <- list(test = test, KNN_res = KNN_res, subdata_cluster_DEG = subdata_cluster_DEG, subdata_cluster = subdata_cluster)
      rm(subdata_cluster, subdata_cluster_DEG, train, test, KNN_res)


      # sub_cluster_label
      sub_cluster_label <- lapply(test_KNN$subdata_cluster_DEG, sub_cluster, FlowSeed = 40)
      # return(sub_cluster_label)

      gc()

      #CSI_Ca_Plot####
      if ("Ca_metric" %in% plot_metric) {
        resAUC <- try(AUC(test = test_KNN$test, KNN_res = test_KNN$KNN_res), silent = T)
        case_problist <- resAUC$case_problist

        num_rows <- round(length(which(!sapply(test_KNN$test, is.null)))/2 + 1e-10)
        if (length(which(!sapply(test_KNN$test, is.null))) == 1) num_rows <- 2
        grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Ca.png"),bg = "white",
                       width = 10, height = 5*num_rows, res=300, units ="in")
        if(length(which(!sapply(test_KNN$test, is.null))) > 1) par(mfrow=c(num_rows, 2))
        CSI_Ca_plot <- list()
        for (j in which(!sapply(test_KNN$test, is.null))) {
          CSI_Ca_plot[[j]] <- ROC_plot(i = j, test = test_KNN$test,case_problist = case_problist )
          if (any(class(CSI_Ca_plot[[j]]) == "try-error")) {
            print(errorplot)
          }
        }
        dev.off()
        print("The plot illustrating the CSI Ca metric is stored in the 'savpath' folder")
        rm(resAUC, case_problist, CSI_Ca_plot)
      }


      #CSI_Cb_Plot####
      if ("Cb_metric" %in% plot_metric){
        CSI_Cb_plot <- try(cluster_plot(data_with_cluster = data_with_cluster), silent = T)
        if (any(class(CSI_Cb_plot) == "try-error")) {
          CSI_Cb_plot <- errorplot
        }
        grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Cb.png"),bg = "white",
                       width = 10, height = 15, res=300, units ="in")
        print(CSI_Cb_plot)
        dev.off()
        print("The plot illustrating the CSI Cb metric is stored in the 'savpath' folder")
        rm(CSI_Cb_plot)
      }

      #CSI_Cc_plot####
      if ("Cc_metric" %in% plot_metric) {
        #Venn
        CS_preres <- try(CS_pre(test_KNN$subdata_cluster), silent = T)
        resCS <- try(CWfun(CS_preres = CS_preres, top = ntop), silent = T)

        CSI_Cc_plot_Venn <- list()
        a <- as.numeric(which(!sapply(resCS$DEGlist, is.null)))
        if (length(a) == 0) {
          message("There is no marker in the list of differentially expressed genes.
                   The program will not be able to display the Cc_metric output.")
        } else {
          num_rows <- round(length(a)/2 + 1e-10)
          if (length(a) == 1) num_rows <- 2
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Cc_Venn.png"),bg = "white",
                         width = 10, height = 5*num_rows, res=300, units ="in")
          for (j in a) {
            CSI_Cc_plot_Venn[[j]] <- try(venn_plot(i = j, CConsistency_marker = resCS$DEGlist), silent = T)
            if (length(CSI_Cc_plot_Venn) != 0 && any(class(CSI_Cc_plot_Venn[[j]]) == "try-error")) {
              CSI_Cc_plot_Venn[[j]] <- errorplot
            }
          }
          if (length(a) > 1){
            gridExtra::grid.arrange(grobs = CSI_Cc_plot_Venn, nrow = num_rows, ncol = 2)
          } else {
            gridExtra::grid.arrange(grobs = CSI_Cc_plot_Venn, nrow = 1, ncol = 1)
          }
          dev.off()
          rm(CS_preres, resCS, CSI_Cc_plot_Venn, a)
        }
        #Volcano
        volcanodata <- lapply(test_KNN$subdata_cluster, volcanovalue, control = as.character(unique(metadata$condition)[1]),
                              case = as.character(unique(metadata$condition)[2]))

        CSI_Cc_plot_Volcano <- list()
        non_null_indices <- as.numeric(which(!sapply(volcanodata, is.null)))
        combined_plot <- NULL
        for (j in non_null_indices) {
          CSI_Cc_plot_Volcano[[j]] <- try(volcano_plot(i = j, volcanodata = volcanodata), silent = T)
          if (any(class(CSI_Cc_plot_Volcano[[j]]) == "try-error")) {
            CSI_Cc_plot_Volcano[[j]] <- errorplot
          }

          if (is.null(combined_plot)) {
            combined_plot <- CSI_Cc_plot_Volcano[[j]]
          } else {
            combined_plot <- combined_plot + CSI_Cc_plot_Volcano[[j]]
          }
        }

        # 根据非空元素数量设置布局
        if (length(non_null_indices) > 1) {
          num_rows <- round(length(non_null_indices) / 2+ 1e-10)
          combined_plot <- combined_plot + patchwork::plot_layout(nrow = num_rows, ncol = 2)
        } else {
          num_rows <- 2
        }

        # 打开 PNG 设备并保存组合后的图形
        grDevices::png(paste0(savepath, "/", studyname, "_", dataset_name, "_CSI_Cc_Volcano.png"),
                       bg = "white",
                       width = 12,
                       height = 5 * num_rows,
                       res = 300,
                       units = "in")
        print(combined_plot)
        dev.off()
        print("The plot illustrating the CSI Cc metric is stored in the 'savpath' folder")
        rm(volcanodata, CSI_Cc_plot_Volcano, combined_plot)
      }

      #CSI_Cd_plot####
      if("Cd_metric" %in% plot_metric){
        if (DEP == "" || is.null(DEP)) {
          print("The parameter 'DEP' is NULL. Please input the known biomarker(s) differentially expressed between two conditions.")
          next
        } else {
          known_marker <- unlist(strsplit(DEP, "\\s*,\\s*"))
          if (length(known_marker) != 0) {
            if (length(known_marker) > 8) {
              widthE_plot <- ceiling(length(known_marker)) * 1.4
            } else widthE_plot <- 12

            CSI_Cd_plot <- try(recall_plot(data_with_cluster = data_with_cluster, known_marker = known_marker), silent = T)
            if (any(class(CSI_Cd_plot) == "try-error")) {
              CSI_Cd_plot <- errorplot
            }
            grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Cd.png"),bg = "white",
                           width = widthE_plot, height = 7, res=300, units ="in")
            print(CSI_Cd_plot)
            dev.off()
            print("The plot illustrating the CSI Cd metric is stored in the 'savpath' folder")
            rm(CSI_Cd_plot)
          }
        }
      }
      rm(data_with_cluster, test_KNN)
      gc()
      print(paste0("The visualization of data processed through workflow: ", dataset_name, " is finished"))
    }
    ###########PTI Visualize#####################
  } else if (studytype == "PTI") {
    try(source("./PTI/load_data2.R"))
    try(source("./PTI/TI_method.R"))

    try(source("./PTI/Bio_con_4.R"))
    try(source("./PTI/time_metric_3.R"))
    try(source("./PTI/Robustness_4.R"))
    try(source("./PTI/Rough_3.R"))

    try(source("./PTI/shift_start.R"))
    try(source("./PTI/calc_spline.R"))
    try(source("./PTI/cycle_pseudotime.R"))
    try(source("./PTI/reverse_pseudotime.R"))
    try(source("./PTI/check_pairs.R"))

    try(source("./PTI/plot.R"))
    try(source("./PTI/ANPELA_FLOWMAP.R"))
    try(source("./PTI/ANPELA_FLOWMAP-function.R"))


    # TIM
    if (missing(TIM)) {
      TIM <- "scorpius_distSpear"
    }
    # Cc_metric
    if (missing(Cc_metric)) {
      Cc_metric <- "Spearman rank correlation"
    } else {
      Cc_metric <- match.arg(Cc_metric)
    }
    #index
    index <- stringr::str_replace_all(colsToUse, "\\(.*", "")

    for (w in 1:length(workflow)){
      dataset_name <- workflow[w]
      print(paste0("The visualization of data processed by workflow '", dataset_name, "' begins."))

      #load data
      if (save_processed_res == "one_folder"){
        load(grep(dataset_name, datapath, value = T))
      } else if (save_processed_res == "one_RData") {
        res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      }

      AP2_processed_D_TI <- try(load.Data(res, index = index, measurement.time = as.matrix(metadata$timepoint), TIM = TIM), silent = T)
      if (class(res) == "try-error") {
        print(paste0("Can not load data for dataset: ", dataset_name))
        next
      }

      TIres <- try(TI(D = AP2_processed_D_TI, method = TIM,
                      dataset_name = dataset_name, clustering.var = clustering.var), silent = T)

      if (class(TIres) == "try-error") {
        print(paste0( "Error in the pseudo-time trajectory inference for dataset: ", dataset_name))
        next
      }


      #PTI_Ca_Plot####
      #time
      if ("Ca_metric" %in% plot_metric) {
        grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Ca_time.png"),
                       bg = "white",width = 10, height = 10, res=300, units ="in")
        if (!grepl("slingshot", TIM)){
          PTI_Ca_Plot_time <- try(dr(TIres, AP2_processed_D_TI), silent = T)
          if (any(class(PTI_Ca_Plot_time) == "try-error")) {
            print(errorplot)
          } else {
            print(PTI_Ca_Plot_time)
          }
        } else {
          PTI_Ca_Plot_time <- try(dr_multi(TIres, AP2_processed_D_TI), silent = T)
        }
        if (any(class(PTI_Ca_Plot_time) == "try-error")) {
          print(errorplot)
        }
        dev.off()

        rm(PTI_Ca_Plot_time)
      }
      #abunf_pt
      result <- TIres
      if (is.list(result$trajectory)) {
        colnames(result$pseudotime) <- paste0("Trajectory_", seq_len(ncol(result$pseudotime)))
        names(result$trajectory) <- paste0("Trajectory_", seq_len(ncol(result$pseudotime)))
        colnames(result$crv1) <- paste0("Trajectory_", seq_len(ncol(result$pseudotime)))
      }
      for(traj_idx in 1:TIres$lineages){
        if (TIres[["cr_method"]] == "run_slingshot"){
          result$pseudotime <- TIres$pseudotime[, traj_idx]
          result$trajectory <- TIres$trajectory[[paste0("Trajectory_", traj_idx)]]
          result$linInd <- as.numeric(traj_idx)
        }
        bio_meaning <- try(suppressWarnings(abund_pt_plot(TIres = result, D = AP2_processed_D_TI)), silent = T)

        if ("Ca_metric" %in% plot_metric) {
          heightA2_plot <- ceiling(length(unique(metadata$timepoint))/4) * 3.5
          widthA2_plot <- 12
          PTI_Ca_Plot_prot <- list()
          combined_plot <- NULL
          for (i in 1:length(index)) {
            PTI_Ca_Plot_prot[[i]] <- try(abund_pt_single(to_plot = bio_meaning$to_plot, prot = index[i], dat = bio_meaning$dat), silent = T)
            if (any(class(PTI_Ca_Plot_prot[[i]]) == "try-error")) {
              PTI_Ca_Plot_prot[[i]] <- errorplot
              print(i)
            }
            if (is.null(combined_plot)) {
              combined_plot <- PTI_Ca_Plot_prot[[i]]
            } else {
              combined_plot <- combined_plot + PTI_Ca_Plot_prot[[i]]
            }
          }
          if (length(index) > 1) {
            num_rows <- round(length(index) / 2+ 1e-10)
            combined_plot <- combined_plot + patchwork::plot_layout(nrow = num_rows, ncol = 2)
          } else {
            num_rows <- 2
          }
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Ca_prot_trajectory",traj_idx,".png"),
                         bg = "white",width = widthA2_plot*2, height = heightA2_plot*num_rows, res=300, units ="in")
          print(combined_plot)
          dev.off()
          print(paste0("The plot illustrating the PTI Ca metric of trajectory",traj_idx, " is stored in the 'savpath' folder"))
          rm(PTI_Ca_Plot_prot, combined_plot)
        }

        # PTI_Cb_plot #####
        if ("Cb_metric" %in% plot_metric) {
          heightB_plot <- ceiling(length(index)/4) * 2.75
          widthB_plot <- 12
          PTI_Cb_plot <- try(roughness_plot(TIres = result, D = AP2_processed_D_TI), silent = T)
          if (any(class(PTI_Cb_plot) == "try-error")) {
            PTI_Cb_plot <- errorplot
          }
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Cb_trajectory",traj_idx,".png"),
                         bg = "white",width = widthB_plot, height = heightB_plot, res=300, units ="in")
          print(PTI_Cb_plot)
          dev.off()
          print(paste0("The plot illustrating the PTI Cb metric of trajectory",traj_idx, " is stored in the 'savpath' folder"))
          rm(PTI_Cb_Plot)
        }

        # PTI_Cc_plot ####
        if ("Cc_metric" %in% plot_metric) {
          Rob <- try(Robustness(TIres = result, D = AP2_processed_D_TI, nruns = 4, cell.subset = 0.8,
                                clustering.var = clustering.var, dataset_name =dataset_name), silent = T)

          PTI_Cc_plot <- BBmisc::suppressAll(try(robustness_new_plot(input_matrix = Rob$input_matrix,
                                                                     finalMatrix = Rob$finalMatrix, method =Cc_metric), silent = T))
          if (any(class(PTI_Cc_plot) == "try-error")) {
            PTI_Cc_plot <- errorplot
          }

          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Cc_trajectory",traj_idx,".png"),
                         bg = "white",width = 10, height = 10, res=300, units ="in")
          print(PTI_Cc_plot)
          dev.off()
          print(paste0("The plot illustrating the PTI Cc metric of trajectory",traj_idx, " is stored in the 'savpath' folder"))
          rm(PTI_Cc_plot,Rob)
        }

        # PTI_Cd_plot ####
        if ("Cd_metric" %in% plot_metric) {
          if (!is.null(pathwayhierarchy) && file.exists(pathwayhierarchy)) {
            heightD1_plot <- ceiling(length(index)/4) * 2.75

            PTI_Cd_plot <- try(bio_meaning$p3, silent = T)
            if (any(class(PTI_Cd_plot) == "try-error")) {
              PTI_Cd_plot <- errorplot
            }

            grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Cd_trajectory",traj_idx,".png"),
                           bg = "white",width = 12, height = heightD1_plot, res=300, units ="in")
            print(PTI_Cd_plot)
            dev.off()
            print(paste0("The plot illustrating the PTI Cd metric of trajectory",traj_idx, " is stored in the 'savpath' folder"))
            rm(PTI_Cd_plot)
          } else {
            print("Please input the filepath of the csv file for parameter 'pathwayhierarchy'.")
          }
        }
        rm(bio_meaning)
      }

      rm(AP2_processed_D_TI, TIres)
      gc()
      print(paste0("The visualization of data processed through workflow: ", dataset_name, " is finished"))
    }
  }
}
