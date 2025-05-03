#' @title Get CSI (Cell Subpopulation Identification) Study Result(s)
#' @description Get_CSIres() performs cell annotation and t-SNE dimensionality reduction on processed SCP data, optionally visualizing the resulting cell populations and saved as a PNG file.
#' @param respath Character, the absolute path of the folder storing the resulting ‘info_saved.RData’ file and the ‘process_res’ folder of the function ‘Process’, ‘FCprocess’ or ‘MCprocess’ when the ‘save_processed_res’ parameter in these functions is set to ‘one_folder’.
#' @param save_processed_res Character, the format of the data processing output files. ‘one_folder’ denotes that successfully processed results will be saved as separate RData files in the ‘process_res’ folder. ‘one_RData’ denotes that all processed results will be saved as one RData file in the ‘process_res’ folder.
#' @param workflow Character, the combinations of data processing methods specified by users according to their research interests.
#'   <br>It is a vector includes one or more method combinations, typically in the format of ‘compensation method name_ transformation method name_ normalization method name_ signal clean method name ‘, for example: c(‘None_Biexponential Transformation_None_None’,’CytoSpill_FlowVS Transformation_None_FlowCut’).
#' @param savepath Character, the absolute path of the folder which will store the results.
#' @param marker_path Character, the absolute file path of the CSV or XLSX file containing the markers for cell type annotation, and detailed format requirements can be found on the website https://github.com/idrblab/ANPELA.
#' @param color Character, the vector of color specifications used to visually distinguish different cell annotation categories within the generated t-SNE plot.
#' @param plot Logical, the logical flag determining whether the t-SNE visualization plot based on the analysis results should be generated and saved as a PNG file.
#' @return An RData file named "workflow_CSIres.RData", recording the CSI results. Optionally a PNG file in addition.
#' @export
#'
#' @examples
#' \donttest{
#' }

Get_CSIres <- function(
    respath, save_processed_res ="one_folder", workflow, savepath = "./ANPELA_res",
    marker_path,color =c("#F39B7FFF","#3C5488FF","#7E6148FF","#B09C85FF","#8491B4FF","#4DBBD5FF","#00A087FF",
                         "#91D1C2FF","#E64B35FF","grey80",RColorBrewer::brewer.pal(12, "Set3")),
    plot = c(T,F)){

  #try(source("./CSI/1readfcs.R"))
  try(source("./CSI/4plot.R"), silent = T)
  library(ggplot2)
  CSIres_list <- list()

  #load data
  datapath <- list.files(paste0(respath, "/process_res/"), pattern = "\\.RData$", full.names = T)
  if (length(datapath)==0){
    stop("The parameter of 'respath' is incorrect. The 'datafile' cannot be loaded.")
  }
  if (save_processed_res == "one_RData") {
    assign("data",load(datapath))
    data <- get(data)
  } else if (save_processed_res == "one_folder") {
    info_saved <- try(load(paste0(respath, "/info_saved.RData")), silent = T)
    if (class(info_saved) == "try-error") {
      stop("The parameter of 'respath' is incorrect. The 'info_saved.RData' cannot be loaded.")
    } else if (info_saved == "info_saved"){
      load(paste0(respath, "/info_saved.RData"))
    }
  }
  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }

  for ( i in 1:length(workflow)){

    dataset_name <- workflow[i]
    if (save_processed_res == "one_folder"){
      load(grep(dataset_name, datapath, value = T))
      condition_info  <- info_saved[["metadata"]][["condition"]]
      colsToUse <- info_saved[["index_TIclass"]]
    } else if (save_processed_res == "one_RData") {
      index <- stringr::str_replace_all(data$index_TIclass, "\\(.*", "")
      res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      condition_info  <- data[["metadata"]][["condition"]]
      colsToUse <- data[["index_TIclass"]]
    }
    if (class(res) == "try-error") {
      print(paste0("Can not load data for dataset: ", dataset_name))
      next
    }

    merge_res <- merge_multi(multi_flowFrame = res,
                             condition_info = condition_info)
    data_j <- try(as.data.frame(merge_res$data, silent = T))
    condition <- merge_res$condition
    data_rd_j <- data_j[,colsToUse] # 仅使用 MCquan_res$index_TIclass 中的列
    data_rd_j$condition <- condition
    rm(merge_res,res)

    set.seed(123) # t-SNE降维的随机种子
    tsne_result <- BBmisc::suppressAll(Rtsne::Rtsne(data_rd_j[,!names(data_rd_j) %in% c("condition")],
                                                    perplexity = 40, dims = 2, verbose = TRUE,
                                                    max_iter = 1000, check_duplicates = F, num_threads= 0)) #10

    tsne_df <- data.frame(
      tsne_x = tsne_result$Y[,1],
      tsne_y = tsne_result$Y[,2],
      condition = as.factor(data_rd_j$condition)
    )

    anno_res <- cell_annotation(data = data_rd_j[,!names(data_rd_j) %in% c("condition")],
                                marker_path = marker_path,
                                colsToUse = colsToUse)
    tsne_df$anno <- anno_res

    CSIres_list[[dataset_name]] <- tsne_df
    save(CSIres_list, file = paste0(savepath, "/", dataset_name,"_CSIres.RData"))
    if (plot){
      p <- ggplot(tsne_df, aes(x = tsne_x, y = tsne_y, color =  .data[["anno"]])) +
        geom_point(alpha = 0.6, size = 0.8) +
        scale_color_manual(values = color,
                           breaks = unique(tsne_df[["anno"]]))+
        theme(
          legend.position = "right",
          panel.background = element_blank(),
          plot.background = element_blank(),  # 设置背景为透明
          panel.grid.major = element_blank(),  # 删除主网格线
          panel.grid.minor = element_blank(),  # 删除次要网格线
          axis.text = element_blank(),         # 删除坐标轴的文字
          axis.ticks = element_blank(),        # 删除坐标轴的刻度线
          axis.line = element_blank(),         # 删除坐标轴的线条
          panel.border = element_blank(),       # 删除边框
          axis.title = element_blank()         # 删除坐标轴的标题
        )

      ggsave(paste0(savepath, "/", dataset_name,"_CSI.png"), p, width = 10, height = 10, dpi = 300, bg="transparent")
    }
  }
}



