
#' @title CSI (Cell Subpopulation Identification) Assess
#' @description CSIassess() assesses processing performance of all workflows which are used while running the function ""FCprocess" or "MCprocess" based on comprehensive criteria (each with distinct underlying theories) from the perspective of CSI studies.
#'
#' @param name Character, the filename of the RData file in the "assess_res" folder which will store the assessment results.
#' @param data Character, the R object resulting from the function "Process", "FCprocess" or "MCprocess", or obtained by loading from the resulting RData file of these funcitons when the `save_processed_res` parameter in these functions is set to "one_RData".
#' @param respath Character, the absolute path of the folder storing the resulting "info_saved.RData" file and the "process_res" folder of the function "Process", "FCprocess" or "MCprocess" when the `save_processed_res` parameter in these functions is set to "one_folder".
#' @param clusteringM Character, the method of clustering the processed data prior to performance assessment, including "FlowSOM" and "PhenoGraph".
#'   <br>**FlowSOM**: a widely used cluster clustering algorithm designed for high-dimensional cytometry data (the number of clusters needs to be specified). ANPELA uses the function “SOM” in R package "FlowSOM" to implement the algorithm.
#'   <br>**PhenoGraph**: a well-designed clustering algorithm developed to define phenotypes in high-dimensional single-cell data (the number of clusters need not be specified). ANPELA uses the function “Rphenograph” in R package "cytofkit".
#' @param ncluster Integer, the number of clusters for meta clustering in FlowSOM.
#'   <br>Only needed when the argument of "clusteringM" is selected as "FlowSOM".
#' @param Phenograph_k Character, the number of nearest neighbours used in PhenoGraph clustering method.
#'   <br>Only needed when the argument of "clusteringM" is selected as "PhenoGraph".
#' @param Ca_metric Character, the assessing metric under Criterion Ca for the "CSI" study type, including "AUC and "F1 score".
#'   <br>**AUC**: a metric that measures the ability of a model to distinguish between classes by calculating the area under the Receiver Operating Characteristic (ROC) curve, where a higher AUC indicates better model performance in terms of classification.
#'   <br>**F1 score**: a metric that evaluates a model’s accuracy by calculating the harmonic mean of precision and recall, where a higher F1 Score indicates a better balance between the model’s ability to identify true positives and avoid false positives.
#' @param Cb_metric Character, the assessing metric under Criterion Cb for the "CSI" study type, including "Silhouette coefficient (SC)", "Xie-Beni index (XB)", "Calinski-Harabasz index (CH)", "Davies-Bouldin index (DB)", "purity" and "Rand index (RI)".
#'   <br>**Silhouette coefficient (SC)**: a metric that measures how similar a sample is to its own cluster compared to other clusters by calculating the difference between the mean intra-cluster distance and the mean nearest-cluster distance.
#'   <br>**Xie-Beni index (XB)**: a metric that assesses the compactness and separation of clusters by calculating the ratio of the sum of squared distances of samples to their cluster centroids over the minimum distance squared between any two cluster centroids.
#'   <br>**Calinski-Harabasz index (CH)**: a metric that evaluates clustering quality by assessing the ratio of between-cluster dispersion to within-cluster dispersion, aiming to maximize inter-cluster distances while minimizing intra-cluster distances.
#'   <br>**Davies-Bouldin index (DB)**: a metric that evaluates clustering results by calculating the average similarity ratio between each cluster and its most similar other cluster, where lower values indicate better clustering performance.
#'   <br>**purity**: a metric that measures the extent to which each cluster contains data points primarily from a single class by identifying the dominant class in each cluster and calculating the proportion of these dominant class counts to the total number of samples.
#'   <br>**Rand index (RI)**: a metric that evaluates the similarity between clustering results and ground truth labels by calculating the ratio of correctly classified sample pairs (both true positives and true negatives) to the total number of pairs.
#' @param Cc_metric Character, the assessing metric under Criterion Cc for the "CSI" study type, including "relative weighted consistency (CWrel)" and "consistency score (CS)".
#'   <br>**relative weighted consistency (CWrel)**: a metric that evaluates the consistency of clustering with respect to a reference partition by weighting the consistency of data point assignments based on their relative importance, with a higher CWrel indicating a stronger alignment with the reference clustering structure.
#'   <br>**consistency score (CS)**: a metric that assesses the stability of clustering results by measuring how consistently data points are assigned to the same clusters across different runs or variations of the clustering process, with a higher CS indicating more stable and reliable clustering.
#' @param ntop Integer, the number of the most differentially expressed markers that are truncated for calculating the CWrel value.
#'   <br>Only needed when the argument of "Cc_metric" is selected as "relative weighted consistency (CWrel)". This value must be less than the number of your selected markers.
#' @param DEP Character, the absolute filepath of the CSV file including the differentially expressed proteins used as the prior knowledge for the fourth criterion.
#'   <br>It is a table of one column without the column name, each table cell includes one protein typically in the format of "channel description (channel name)", for example: "CD20(FITC.A)".
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param marker_path Character, the absolute file path of the CSV or XLSX file containing the markers for cell type annotation, and detailed format requirements can be found on the website https://github.com/idrblab/ANPELA.
#' @param known_celltype_path Character, the absolute file path of the CSV file containing the gold-standard cell type annotation results, with the first column being cell IDs and the second column being cell types, and these cell types should correspond pre-cisely to those in the ‘marker_path’ file.
#' @param save_processed_res Character, the format of the data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
#' @param savepath Character, the absolute path of the folder which will store the assessment results.
#' @return The **assess_res** folder stores the assessment output file named `name`**_assess.RData**, which contains 2 lists, "table" and "table2", providing the raw scores for different assessment criteria and performance assessment levels categorized by thresholds, respectively.
#'   <br>In addition, the file **log.txt** is also generated simultaneously, recording the processing details.
#' @export
#'
#' @examples
#' \donttest{
#' }


CSIassess <- function(name = "result", data, respath,
                      clusteringM = c("FlowSOM", "PhenoGraph","Mclust"),
                      Phenograph_k = 30,
                      ncluster = 8,
                      Ca_metric = "AUC",
                      Cb_metric = "Silhouette coefficient (SC)",
                      Cc_metric = "relative weighted consistency (CWrel)",
                      ntop = NULL,
                      DEP = NULL,
                      marker_path = NULL, known_celltype_path = NULL,
                      save_processed_res = "one_folder",
                      savepath = "./ANPELA_res",
                      cores = floor(parallel::detectCores()/2), ...) {

  # data & info_saved & process_res
  if (save_processed_res == "one_RData") {
    if (missing(data)) {
      stop("The parameter of 'data' is missing.")
    }else if (any(sapply(data, is.null))) {
      stop("The parameter of 'data' is incorrect. Please check the data.")
    }
  } else if (save_processed_res == "one_folder") {
    info_saved <- try(load(paste0(respath, "/info_saved.RData")), silent = T)
    datapath <- list.files(paste0(respath, "/process_res/"), pattern = "\\.RData$", full.names = T)
    if (class(info_saved) == "try-error") {
      stop("The parameter of 'respath' is incorrect. The 'info_saved.RData' cannot be loaded.")
    } else if (info_saved == "info_saved"){
      load(paste0(respath, "/info_saved.RData"))
    }
    if (length(datapath)==0){
      stop("The parameter of 'respath' is incorrect. The 'datafile' cannot be loaded.")
    }
  }


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


  # Ca_metric
  if (missing(Ca_metric)) {
    Ca_metric <- "AUC"
  } else {
    Ca_metric <- match.arg(Ca_metric)
  }


  # Cb_metric
  if (missing(Cb_metric)) {
    Cb_metric <- "Silhouette coefficient (SC)"
  } else {
    Cb_metric <- match.arg(Cb_metric)
  }


  # Cc_metric
  if (missing(Cc_metric)) {
    Cc_metric <- "relative weighted consistency (CWrel)"
  } else {
    Cc_metric <- match.arg(Cc_metric)
  }


  # ntop
  if (Cc_metric == "relative weighted consistency (CWrel)") {
    if (save_processed_res == "one_RData") {
      if (missing(ntop)) {
        ntop <- floor(length(data$index_TIclass)/2)
      } else if (ntop >= length(data$index_TIclass) || ntop %% 1 != 0 || ntop <= 0) {
        stop("The value of 'ntop' is incorrect. It should be positive whole number and less than the number of your selected markers.")
      }
    } else if (save_processed_res == "one_folder") {
      if (missing(ntop)) {
        ntop <- floor(length(info_saved$index_TIclass)/2)
      } else if (ntop >= length(info_saved$index_TIclass) || ntop %% 1 != 0 || ntop <= 0) {
        stop("The value of 'ntop' is incorrect. It should be positive whole number and less than the number of your selected markers.")
      }
    }
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

  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }

  cat("The program is running. Please wait patiently.")



  if (save_processed_res == "one_RData") {
    # parallel start
    opts <- list(progress = function(n) setTxtProgressBar(txtProgressBar(min = 1, max = length(data$AP2_pro1_frame_classTI), style = 3), n))
    cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/CSIassess_log.txt"))
    doSNOW::registerDoSNOW(cl)
    time = proc.time()


    table <- foreach::foreach(i = 1:length(data$AP2_pro1_frame_classTI), .options.snow = opts,
                              .packages = c("Rphenograph", "dplyr","mclust"), .combine = rbind) %dopar% {

                                try(source("./CSI/1readfcs.R"))
                                try(source("./CSI/2cluster.R"))
                                try(source("./CSI/3criteria.R"))
                                try(source("./CSI/4plot.R"))


                                # AP2_processed_D_class
                                res <- data$AP2_pro1_frame_classTI[[i]]
                                if (is.null(res)) {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                  return(res)
                                }
                                names(res) <- limma::removeExt(basename(data$dataFileNames), sep=".")
                                AP2_processed_D_class <- res
                                rm(res)


                                # AP2_processed_data_class0
                                data1 <- try(as.data.frame(readfcs_multi(fcsFiles = AP2_processed_D_class, mergeMethod = "all")), silent = T)
                                if (class(data1) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                  return(res)
                                }
                                # make.names(colnames(data1))
                                colnames(data1) <- sub("<", "_", colnames(data1))
                                colnames(data1) <- sub(">", "_", colnames(data1))
                                rm(AP2_processed_D_class)

                                condition <- data$metadata
                                data_with_filename <- cbind(data1, "filename" = sub("_[0-9]*$","", row.names(data1)))
                                condition_label <- as.character(condition[match(data_with_filename[,"filename"], condition$filename),2])
                                data_with_condition <- cbind(data_with_filename, "condition" = condition_label)
                                AP2_processed_data_class0 <- data_with_condition
                                rm(data1, condition, data_with_filename, condition_label, data_with_condition)


                                # AP2_processed_data_class
                                AP2_processed_data_class <- AP2_processed_data_class0[, c(data$index_TIclass, "filename", "condition")]
                                rm(AP2_processed_data_class0)


                                # cluster_label
                                cluster_label <- try(data_cluster(data = AP2_processed_data_class[,1:(dim(AP2_processed_data_class)[2]-2)], method = clusteringM, Phenograph_k = Phenograph_k, FlowSOM_k = ncluster, FlowSeed = 40), silent = T)
                                if (class(cluster_label) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                  return(res)
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

                                if (Ca_metric == "AUC") {
                                  # Criterion A Accuracy-AUC
                                  resAUC <- try(AUC(test = test_KNN$test, KNN_res = test_KNN$KNN_res), silent = T)
                                  Cauc <- try(round(sum(resAUC$auc, na.rm = TRUE)/length(test_KNN[["subdata_cluster"]]), 5), silent = T)
                                  if (class(Cauc) != "numeric") {
                                    Cauc <- NA
                                  }
                                  Ca <- Cauc
                                  rm(resAUC, Cauc)
                                } else if (Ca_metric == "F1 score") {
                                  # Criterion A Accuracy-F1 score
                                  resF1 <- try(F1_score(test = test_KNN$test, KNN_res = test_KNN$KNN_res, label = AP2_processed_data_class$condition), silent = T)
                                  CF1_score <- try(round(max(sapply(resF1, mean, na.rm = TRUE)), 5), silent = T)
                                  if (class(CF1_score) != "numeric") {
                                    CF1_score <- NA
                                  }
                                  Ca <- CF1_score

                                  rm(resF1, CF1_score)
                                }

                                if (Cb_metric == "Silhouette coefficient (SC)") {
                                  # Criterion B Silhouette coefficient (SC)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Silhouette")[[1]], silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  } else {
                                    Cb <- round((Cb + 1)/2, 5)
                                  }
                                } else if (Cb_metric == "Xie-Beni index (XB)") {
                                  # Criterion B Xie-Beni index (XB)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Xie_Beni")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "Calinski-Harabasz index (CH)") {
                                  # Criterion B Calinski-Harabasz index (CH)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Calinski_Harabasz")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "Davies-Bouldin index (DB)") {
                                  # Criterion B Davies-Bouldin index (DB)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Davies_Bouldin")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "purity") {
                                  # Criterion B Precision-purity
                                  respurity <- try(Purity(data = test_KNN$subdata_cluster_DEG, sub_cluster_label = sub_cluster_label, FlowSeed = 40), silent = T)
                                  Cpurity <- try(round(mean(respurity, na.rm = TRUE), 5), silent = T)
                                  if (class(Cpurity) != "numeric") {
                                    Cpurity <- NA
                                  }
                                  Cb <- Cpurity
                                  rm(respurity, Cpurity)
                                } else if (Cb_metric == "Rand index (RI)") {
                                  # Criterion B Precision-RI
                                  resCRI <- try(RI(data = test_KNN$subdata_cluster_DEG, sub_cluster_label = sub_cluster_label, FlowSeed = 40), silent = T)
                                  CRI <- try(round(mean(resCRI, na.rm = TRUE), 5), silent = T)
                                  if (class(CRI) != "numeric") {
                                    CRI <- NA
                                  }
                                  Cb <- CRI
                                  rm(resCRI, CRI)
                                }
                                rm(sub_cluster_label, AP2_processed_data_class, cluster_label, datab)


                                # Criterion C Robustness-CS_pre
                                CS_preres <- try(CS_pre(subdata_cluster = test_KNN$subdata_cluster), silent = T)
                                # Robustness-CS, CW
                                if (Cc_metric == "consistency score (CS)") {
                                  resCS <- try(CSfun(CS_preres = CS_preres), silent = T)
                                } else if (Cc_metric == "relative weighted consistency (CWrel)") {
                                  resCS <- try(CWfun(CS_preres = CS_preres, top = ntop), silent = T)
                                }
                                CS <- try(round(mean(resCS$consistency, na.rm = TRUE), 5), silent = T)
                                if (class(CS) != "numeric") {
                                  CS <- NA
                                }
                                Cc <- CS
                                rm(test_KNN, CS_preres, resCS, CS)


                                # Criterion D Biological Meaning
                                if ((is.null(marker_path)||is.null(known_celltype_path))& (DEP == "" || is.null(DEP))) {
                                  Cd <- NA
                                } else if (!is.null(marker_path) && !is.null(known_celltype_path)) {
                                  CRecall <- try(round(AP2_Recall(data_with_cluster = data_with_cluster,
                                                                  marker_path = marker_path, known_celltype_path = known_celltype_path), 5), silent = T)
                                  if (class(CRecall) != "numeric") {
                                    CRecall <- NA
                                  }
                                  Cd <- CRecall
                                } else if (!DEP == "" && !is.null(DEP)) {
                                  known_marker <- unlist(strsplit(DEP, "\\s*,\\s*"))
                                  CRecall <- try(round(AP2_Recall(data_with_cluster = data_with_cluster, known_marker = known_marker), 5), silent = T)
                                  if (class(CRecall) != "numeric") {
                                    CRecall <- NA
                                  }
                                  Cd <- CRecall
                                  rm(known_marker, CRecall)
                                }
                                rm(data_with_cluster)
                                gc()

                                res <- data.frame(Ca, Cb, Cc, Cd)
                                rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                return(res)
                              }


    parallel::stopCluster(cl)
    print(proc.time()-time)
    # parallel end
  } else if (save_processed_res == "one_folder") {
    # parallel start
    opts <- list(progress = function(n) setTxtProgressBar(txtProgressBar(min = 1, max = length(datapath), style = 3), n))#
    cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/CSIassess_log.txt"))
    doSNOW::registerDoSNOW(cl)
    time = proc.time()


    table <- foreach::foreach(i = 1:length(datapath), .options.snow = opts,
                              .packages = c("Rphenograph", "dplyr","mclust"), .combine = rbind) %dopar% {

                                try(source("./CSI/1readfcs.R"))
                                try(source("./CSI/2cluster.R"))
                                try(source("./CSI/3criteria.R"))
                                try(source("./CSI/4plot.R"))

                                # AP2_processed_D_class
                                load(datapath[i])

                                if (is.null(res)) {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]
                                  return(res)
                                }
                                names(res) <- limma::removeExt(basename(info_saved$dataFileNames), sep=".")
                                AP2_processed_D_class <- res
                                rm(res)


                                # AP2_processed_data_class0
                                data1 <- try(as.data.frame(readfcs_multi(fcsFiles = AP2_processed_D_class, mergeMethod = "all")), silent = T)
                                if (class(data1) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]
                                  return(res)
                                }
                                # make.names(colnames(data1))
                                colnames(data1) <- sub("<", "_", colnames(data1))
                                colnames(data1) <- sub(">", "_", colnames(data1))
                                rm(AP2_processed_D_class)

                                condition <- info_saved$metadata
                                data_with_filename <- cbind(data1, "filename" = sub("_[0-9]*$","", row.names(data1)))
                                condition_label <- as.character(condition[match(data_with_filename[,"filename"], condition$filename),2])
                                data_with_condition <- cbind(data_with_filename, "condition" = condition_label)
                                AP2_processed_data_class0 <- data_with_condition
                                rm(data1, condition, data_with_filename, condition_label, data_with_condition)


                                # AP2_processed_data_class
                                AP2_processed_data_class <- AP2_processed_data_class0[, c(info_saved$index_TIclass, "filename", "condition")]
                                rm(AP2_processed_data_class0)


                                # cluster_label
                                cluster_label <- try(data_cluster(data = AP2_processed_data_class[,1:(dim(AP2_processed_data_class)[2]-2)],
                                                                  method = clusteringM, Phenograph_k = Phenograph_k,
                                                                  FlowSOM_k = ncluster, FlowSeed = 40), silent = T)
                                if (class(cluster_label) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]
                                  return(res)
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

                                gc()

                                if (Ca_metric == "AUC") {
                                  # Criterion A Accuracy-AUC
                                  resAUC <- try(AUC(test = test_KNN$test, KNN_res = test_KNN$KNN_res), silent = T)
                                  Cauc <- try(round(sum(resAUC$auc, na.rm = TRUE)/length(test_KNN[["subdata_cluster"]]), 5), silent = T)
                                  if (class(Cauc) != "numeric") {
                                    Cauc <- NA
                                  }
                                  Ca <- Cauc
                                  rm(resAUC, Cauc)
                                } else if (Ca_metric == "F1 score") {
                                  # Criterion A Accuracy-F1 score
                                  resF1 <- try(F1_score(test = test_KNN$test, KNN_res = test_KNN$KNN_res, label = AP2_processed_data_class$condition), silent = T)
                                  CF1_score <- try(round(max(sapply(resF1, mean, na.rm = TRUE)), 5), silent = T)
                                  if (class(CF1_score) != "numeric") {
                                    CF1_score <- NA
                                  }
                                  Ca <- CF1_score

                                  rm(resF1, CF1_score)
                                }

                                if (Cb_metric == "Silhouette coefficient (SC)") {
                                  # Criterion B Silhouette coefficient (SC)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Silhouette")[[1]], silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  } else {
                                    Cb <- round((Cb + 1)/2, 5)
                                  }
                                } else if (Cb_metric == "Xie-Beni index (XB)") {
                                  # Criterion B Xie-Beni index (XB)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Xie_Beni")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "Calinski-Harabasz index (CH)") {
                                  # Criterion B Calinski-Harabasz index (CH)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Calinski_Harabasz")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "Davies-Bouldin index (DB)") {
                                  # Criterion B Davies-Bouldin index (DB)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Davies_Bouldin")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "purity") {
                                  # Criterion B Precision-purity
                                  respurity <- try(Purity(data = test_KNN$subdata_cluster_DEG, sub_cluster_label = sub_cluster_label, FlowSeed = 40), silent = T)
                                  Cpurity <- try(round(mean(respurity, na.rm = TRUE), 5), silent = T)
                                  if (class(Cpurity) != "numeric") {
                                    Cpurity <- NA
                                  }
                                  Cb <- Cpurity
                                  rm(respurity, Cpurity)
                                } else if (Cb_metric == "Rand index (RI)") {
                                  # Criterion B Precision-RI
                                  resCRI <- try(RI(data = test_KNN$subdata_cluster_DEG, sub_cluster_label = sub_cluster_label, FlowSeed = 40), silent = T)
                                  CRI <- try(round(mean(resCRI, na.rm = TRUE), 5), silent = T)
                                  if (class(CRI) != "numeric") {
                                    CRI <- NA
                                  }
                                  Cb <- CRI
                                  rm(resCRI, CRI)
                                }
                                rm(sub_cluster_label, AP2_processed_data_class, cluster_label, datab)


                                # Criterion C Robustness-CS_pre
                                CS_preres <- try(CS_pre(subdata_cluster = test_KNN$subdata_cluster), silent = T)
                                # Robustness-CS, CW
                                if (Cc_metric == "consistency score (CS)") {
                                  resCS <- try(CSfun(CS_preres = CS_preres), silent = T)
                                } else if (Cc_metric == "relative weighted consistency (CWrel)") {
                                  resCS <- try(CWfun(CS_preres = CS_preres, top = ntop), silent = T)
                                }
                                CS <- try(round(mean(resCS$consistency, na.rm = TRUE), 5), silent = T)
                                if (class(CS) != "numeric") {
                                  CS <- NA
                                }
                                Cc <- CS
                                rm(test_KNN, CS_preres, resCS, CS)


                                # Criterion D Biological Meaning
                                if ((is.null(marker_path)||is.null(known_celltype_path))& (DEP == "" || is.null(DEP))) {
                                  Cd <- NA
                                } else if (!is.null(marker_path) && !is.null(known_celltype_path)) {
                                  CRecall <- try(round(AP2_Recall(data_with_cluster = data_with_cluster,
                                                                  marker_path = marker_path, known_celltype_path = known_celltype_path), 5))#, silent = T)
                                  if (class(CRecall) != "numeric") {
                                    CRecall <- NA
                                  }
                                  Cd <- CRecall
                                } else if (!DEP == "" && !is.null(DEP)) {
                                  known_marker <- unlist(strsplit(DEP, "\\s*,\\s*"))
                                  CRecall <- try(round(AP2_Recall(data_with_cluster = data_with_cluster, known_marker = known_marker), 5), silent = T)
                                  if (class(CRecall) != "numeric") {
                                    CRecall <- NA
                                  }
                                  Cd <- CRecall
                                  rm(known_marker, CRecall)
                                }
                                rm(data_with_cluster)
                                gc()

                                res <- data.frame(Ca, Cb, Cc, Cd)
                                rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]
                                return(res)
                              }


    parallel::stopCluster(cl)
    print(proc.time()-time)
    # parallel end
  }


  colnames(table) <- c("Accuracy", "Tightness", "Robustness", "Correspondence")
  table2 <- table
  table2["Accuracy"][table2["Accuracy"] > 0.7] <- 10
  table2["Accuracy"][table2["Accuracy"] <= 0.7] <- 4

  table2["Tightness"][table2["Tightness"] > 0.5] <- 10
  table2["Tightness"][table2["Tightness"] <= 0.5] <- 4

  table2["Robustness"][table2["Robustness"] > 0.35] <- 10
  table2["Robustness"][table2["Robustness"] <= 0.35] <- 4

  if (!is.null(marker_path) && !is.null(known_celltype_path)) {
    table2["Correspondence"][table2["Correspondence"] > 0.7] <- 10
    table2["Correspondence"][table2["Correspondence"] <= 0.5] <- 4
  } else {
    table2["Correspondence"][table2["Correspondence"] > 0.5] <- 10
    table2["Correspondence"][table2["Correspondence"] <= 0.5] <- 4
  }

  assess_res <- list(table = table, table2 = table2)
  if (!dir.exists(paste0(savepath, "/assess_res"))) {
    dir.create(paste0(savepath, "/assess_res"), recursive = T)
  }
  save(assess_res, file = paste0(savepath, "/assess_res/", name, "_assess.RData"))
  return(assess_res)
}

