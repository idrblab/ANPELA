
#' @title PTI (Pseudotime Trajectory Inference) Assess
#' @description PTIassess() assesses processing performance of all workflows which are used while running the function "FCprocess" or "MCprocess" based on comprehensive criteria (each with distinct underlying theories) from the perspective of PTI studies.
#'
#' @param name Character, the filename of the RData file in the "assess_res" folder which will store the assessment results.
#' @param data Character, the R object resulting from the function "Process", "FCprocess" or "MCprocess", or obtained by loading from the resulting RData file of these funcitons when the `save_processed_res` parameter in these functions is set to "one_RData".
#' @param respath Character, the absolute path of the folder storing the resulting "info_saved.RData" file and the "process_res" folder of the function "Process", "FCprocess" or "MCprocess" when the `save_processed_res` parameter in these functions is set to "one_folder".
#' @param TIM Character, the method of trajectory inference for the processed data prior to performance assessment, consisted of tra-jectory reconstruction and data space representation, including ‘scorpius_distSpear’, ‘scorpius_distPear’, ‘scorpius_distEucl’, ‘scorpius_distManh’, ‘slingshot_tSNE’, ‘slingshot_FLOWMAP’, ‘slingshot_PCA’, ‘slingshot_diffMaps’.
#' @param Cc_metric Character, the assessing metric under Criterion Cc for the "PTI" study type, including "Spearman correlation" and "Kendall Rank Correlation".
#'   <br>**Spearman correlation**: a metric that measures the monotonic relationship between two ranked variables by measuring how well the relationship between the variables can be described by a monotonic function, with values closer to ±1 indicating a stronger monotonic association.*
#'   <br>**Kendall Rank Correlation**: a metric that evaluates the similarity between two rankings by calculating the proportion of concordant and discordant pairs, with higher values indicating a stronger agreement in the relative ordering between the variables.
#' @param pathwayhierarchy Character, the absolute filepath of the pathway hierarchy file.
#' @param clustering.var Character, the vector naming channels to be used to calculate distances/differences between cells for clustering (if re-quested) and edge-drawing steps.
#'   <br>Only needed when ‘slingshot_FLOWMAP’ is included in the parameter of ‘TIM’.
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param save_processed_res Character, the format of the data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
#' @param savepath Character, the absolute path of the folder which will store the assessment results.
#' @return The **assess_res** folder stores the assessment output file named `name`**_assess.RData**, which contains 2 lists, "table" and "table2", providing the raw scores for different assessment criteria and performance assessment levels categorized by thresholds, respectively.
#'   <br>In addition, the file **log.txt** is also generated simultaneously, recording the processing details.
#' @export
#'
#' @examples
#' \donttest{
#' }


PTIassess <- function(name = "result", data, respath,
                      TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh",
                              "slingshot_tSNE", "slingshot_FLOWMAP", "slingshot_PCA", "slingshot_diffMaps"),
                      Cc_metric = c("Spearman rank correlation", "Kendall rank correlation"),
                      pathwayhierarchy = NULL,
                      save_processed_res = "one_folder",
                      savepath = "./ANPELA_res",
                      clustering.var = NULL,

                      cores = floor(parallel::detectCores()/2),
                      ...) {
  # data
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


  # TIM
  if (missing(TIM)) {
    TIM <- "scorpius_distSpear"
  } else {
    TIM <- match.arg(TIM)
  }


  # Cc_metric
  if (missing(Cc_metric)) {
    Cc_metric <- "Spearman rank correlation"
  } else {
    Cc_metric <- match.arg(Cc_metric)
  }

  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }

  cat("The program is running. Please wait patiently.")

  if (save_processed_res == "one_RData") {
    # parallel start
    opts <- list(progress = function(n) setTxtProgressBar(txtProgressBar(min = 1, max = length(data$AP2_pro1_frame_classTI), style = 3), n))
    cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/PTIassess_log.txt"))

    doSNOW::registerDoSNOW(cl)
    time = proc.time()

    #length(data$AP2_pro1_frame_classTI)
    table <- foreach::foreach(i = 1:length(data$AP2_pro1_frame_classTI), .options.snow = opts,
                              .packages = c("flowCore", "foreach", "dplyr", "igraph", "mclust"), .combine = rbind) %dopar% {
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

                                # AP2_processed_D_TI
                                index <- stringr::str_replace_all(data$index_TIclass, "\\(.*", "")
                                dataset_name <- names(data$AP2_pro1_frame_classTI)[i]#flowmap
                                res <- data$AP2_pro1_frame_classTI[[i]]
                                if ("condition" %in% colnames(data$metadata)){
                                  res <- try(load.Data(res, index = index, measurement.condition = as.matrix(data$metadata$condition),
                                                       measurement.time = as.matrix(data$metadata$timepoint), TIM = TIM), silent = T)
                                } else {
                                  res <- try(load.Data(res, index = index, measurement.time = as.matrix(data$metadata$timepoint), TIM = TIM), silent = T)
                                }

                                if (class(res) == "try-error") {
                                  res <- data.frame(Tm = NA, R_pvalue = NA, Rob = NA, BC = NA)
                                  rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                  return(res)
                                }

                                AP2_processed_D_TI <- res
                                rm(index, res)

                                # TIres
                                TIres <- try(TI(AP2_processed_D_TI, method = TIM,
                                                dataset_name = dataset_name, clustering.var = clustering.var)) #flowmap
                                print(paste0("TIres","_",i))
                                if(class(TIres) == "try-error") {
                                  res <- data.frame(Tm = NA, R_pvalue = NA, Rob = NA, BC = NA)
                                  rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                  return(res)
                                }



                                # Criterion A Consistency-Time
                                Tm <- try(time_metric(TIres, AP2_processed_D_TI, nruns = 100), silent = T)
                                if (class(Tm) != "numeric") {
                                  Tm <- NA
                                } else {
                                  Tm <- round(Tm, 5)
                                }


                                # Criterion B Roughness-Roughness
                                R <- try(Rough(TIres, AP2_processed_D_TI), silent = T)
                                R_pvalue <- try(R$p.value, silent = T)
                                if (class(R_pvalue) != "numeric") {
                                  R_pvalue <- NA
                                } else {
                                  if (R_pvalue > 0.05) {
                                    R_pvalue <- 0
                                  } else {
                                    R_pvalue <- round(1 - 20*R_pvalue, 5)
                                  }
                                }
                                rm(R)


                                # Criterion C Robustness-Robustness
                                Rob0 <- try(Robustness(TIres, AP2_processed_D_TI, nruns = 4, cell.subset = 0.8,clustering.var = clustering.var), silent = T)
                                if (class(Rob0) != "list") {
                                  Rob <- NA
                                } else {
                                  Rob <- switch (Cc_metric,
                                                 "Spearman rank correlation" = round(Rob0$Robustness_result[1], 5),
                                                 "Kendall rank correlation" = round(Rob0$Robustness_result[2], 5)
                                  )
                                }

                                rm(Rob0)


                                # Criterion D Biological Meaning-Biological consistency
                                if (!is.null(pathwayhierarchy) && file.exists(pathwayhierarchy)) {
                                  BC <- try(suppressWarnings(Bio_con(AP2_processed_D_TI, Pathway_Hierarchy_file = pathwayhierarchy, nruns = 3, dr_method = TIres$dr_method, TIres = TIres)), silent = T)
                                  if (class(BC) != "numeric") {
                                    BC <- NA
                                  } else {
                                    BC <- round(BC, 5)
                                  }
                                } else if (is.null(pathwayhierarchy) && "condition" %in% colnames(data$metadata)){
                                  BC <- try(suppressWarnings(Bio_con(AP2_processed_D_TI,  nruns = 3, dr_method = TIres$dr_method, TIres = TIres)), silent = T)
                                  if (class(BC) != "numeric") {
                                    BC <- NA
                                  } else {
                                    BC <- round(BC, 5)
                                  }
                                } else BC <- NA


                                rm(AP2_processed_D_TI, TIres)
                                gc()

                                res <- data.frame(Tm, R_pvalue, Rob, BC)
                                rownames(res) <- names(data$AP2_pro1_frame_classTI)[i]
                                print(paste0("res","_",i))
                                return(res)
                              }

    parallel::stopCluster(cl)
    print(proc.time()-time)
    # parallel end
  }else if (save_processed_res == "one_folder") {
    # parallel start
    opts <- list(progress = function(n) setTxtProgressBar(txtProgressBar(min = 1, max = length(datapath),  style = 3), n))
    cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/PTIassess_log.txt"))

    doSNOW::registerDoSNOW(cl)

    time = proc.time()


    table <- foreach::foreach(i = 1:length(datapath), .options.snow = opts,
                              .packages = c("flowCore", "foreach", "dplyr", "igraph", "mclust"), .combine = rbind) %dopar% {

                                # # 用一个全局变量防止重复 sink
                                # if (!exists(".log_started", envir = .GlobalEnv)) {
                                #   assign(".log_started", TRUE, envir = .GlobalEnv)
                                #
                                #   pid <- Sys.getpid()
                                #   log_file <- paste0(savepath, "/log_worker_", pid, ".txt")
                                #
                                #   # 尝试打开 sink
                                #   try({
                                #     sink(log_file, split = TRUE)
                                #     sink(log_file, type = "message", append = TRUE)
                                #   }, silent = TRUE)
                                # }
                                #
                                # # 日志输出
                                # pid <- Sys.getpid()
                                # cat(sprintf("[%s] Worker PID %d started task %d\n", Sys.time(), pid, i))

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
                                # AP2_processed_D_TI
                                index <- stringr::str_replace_all(info_saved$index_TIclass, "\\(.*", "")
                                load(datapath[i])
                                dataset_name <-limma::removeExt(basename(datapath), sep=".")[i]#flowmap
                                if ("condition" %in% colnames(info_saved$metadata)){
                                  res <- try(load.Data(res, index = index, measurement.condition = as.matrix(info_saved$metadata$condition),
                                                       measurement.time = as.matrix(info_saved$metadata$timepoint), TIM = TIM), silent = T)
                                } else {
                                  res <- try(load.Data(res, index = index, measurement.time = as.matrix(info_saved$metadata$timepoint), TIM = TIM), silent = T)
                                }

                                if (class(res) == "try-error") {
                                  res <- data.frame(Tm = NA, R_pvalue = NA, Rob = NA, BC = NA)
                                  rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]
                                  return(res)
                                }

                                AP2_processed_D_TI <- res
                                rm(index, res)

                                # TIres
                                TIres <- try(TI(D = AP2_processed_D_TI, method = TIM,  dataset_name = dataset_name,clustering.var =clustering.var))
                                print(paste0("TIres","_",i))
                                if(class(TIres) == "try-error") {
                                  res <- data.frame(Tm = NA, R_pvalue = NA, Rob = NA, BC = NA)
                                  rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]
                                  print(paste0("TIres","_",i,"error:",TIres))
                                  return(res)
                                }



                                # Criterion A Consistency-Time
                                Tm <- try(time_metric(TIres, AP2_processed_D_TI, nruns = 100), silent = T)
                                if (class(Tm) != "numeric") {
                                  Tm <- NA
                                } else {
                                  Tm <- round(Tm, 5)
                                }
                                print(paste0("TIres","_",i,"_Ca:",Tm))

                                # Criterion B Roughness-Roughness
                                R <- try(Rough(TIres, AP2_processed_D_TI), silent = T)
                                R_pvalue <- try(R$p.value, silent = T)
                                if (class(R_pvalue) != "numeric") {
                                  R_pvalue <- NA
                                } else {
                                  if (R_pvalue > 0.05) {
                                    R_pvalue <- 0
                                  } else {
                                    R_pvalue <- round(1 - 20*R_pvalue, 5)
                                  }
                                }
                                rm(R)
                                print(paste0("TIres","_",i,"_Cb:",R_pvalue))

                                # Criterion C Robustness-Robustness
                                Rob0 <- try(Robustness(TIres, AP2_processed_D_TI, nruns = 4, cell.subset = 0.8, clustering.var = clustering.var,
                                                       dataset_name = dataset_name), silent = T)
                                if (class(Rob0) != "list") {
                                  Rob <- NA
                                } else {
                                  Rob <- switch (Cc_metric,
                                                 "Spearman rank correlation" = round(Rob0$Robustness_result[1], 5),
                                                 "Kendall rank correlation" = round(Rob0$Robustness_result[2], 5)
                                  )
                                }

                                rm(Rob0)
                                print(paste0("TIres","_",i,"_Cc:",Rob))

                                # Criterion D Biological Meaning-Biological consistency
                                if (!is.null(pathwayhierarchy) && file.exists(pathwayhierarchy)) {
                                  BC <- try(suppressWarnings(Bio_con(AP2_processed_D_TI, Pathway_Hierarchy_file = pathwayhierarchy, nruns = 3, dr_method = TIres$dr_method, TIres = TIres)), silent = T)
                                  if (class(BC) != "numeric") {
                                    BC <- NA
                                  } else {
                                    BC <- round(BC, 5)
                                  }
                                } else if(is.null(pathwayhierarchy) && "condition" %in% colnames(info_saved$metadata)){
                                  BC <- try(suppressWarnings(Bio_con(AP2_processed_D_TI,  nruns = 3, dr_method = TIres$dr_method, TIres = TIres)), silent = T)
                                  if (class(BC) != "numeric") {
                                    BC <- NA
                                  } else {
                                    BC <- round(BC, 5)
                                  }
                                } else BC <- NA
                                print(paste0("TIres","_",i,"_Cd:",BC))

                                rm(AP2_processed_D_TI, TIres)
                                gc()

                                res <- data.frame(Tm, R_pvalue, Rob, BC)
                                rownames(res) <- limma::removeExt(basename(datapath), sep=".")[i]

                                return(res)
                              }

    # # 关闭 sink（每个 worker 退出时关闭）
    # parallel::clusterEvalQ(cl, {
    #   try(sink(), silent = TRUE)
    #   try(sink(type = "message"), silent = TRUE)
    #   NULL
    # })
    #
    #
    parallel::stopCluster(cl)
    print(proc.time()-time)
    # parallel end
  }

  colnames(table) <- c("Conformance", "Smoothness", "Robustness", "Correspondence")
  table2 <- table
  table2["Conformance"][table2["Conformance"] > 0.6] <- 10
  table2["Conformance"][table2["Conformance"] <= 0.6] <- 4

  table2["Smoothness"][table2["Smoothness"] > 0.8] <- 10
  table2["Smoothness"][table2["Smoothness"] <= 0.8] <- 4

  table2["Robustness"][table2["Robustness"] > 0.5] <- 10
  table2["Robustness"][table2["Robustness"] <= 0.5] <- 4

  table2["Correspondence"][table2["Correspondence"] > 0.6] <- 10
  table2["Correspondence"][table2["Correspondence"] <= 0.6] <- 4

  assess_res <- list(table = table, table2 = table2)
  if (!dir.exists(paste0(savepath, "/assess_res"))) {
    dir.create(paste0(savepath, "/assess_res"), recursive = T)
  }
  save(assess_res, file = paste0(savepath, "/assess_res/", name, "_assess.RData"))
  return(assess_res)
}
