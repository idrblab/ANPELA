#' @title Get PTI (Pseudotime Trajectory Inference) Study Result(s)
#' @description Get_PTIres() applies selected trajectory inference algorithms to processed SCP data, optionally visualizing the resulting cellular trajectories and saved as a PNG file.
#' @param respath Character, the absolute path of the folder storing the resulting ‘info_saved.RData’ file and the ‘process_res’ folder of the function ‘Process’, ‘FCprocess’ or ‘MCprocess’ when the ‘save_processed_res’ parameter in these functions is set to ‘one_folder’.
#' @param save_processed_res Character, the format of the data processing output files. ‘one_folder’ denotes that successfully processed results will be saved as separate RData files in the ‘process_res’ folder. ‘one_RData’ denotes that all processed results will be saved as one RData file in the ‘process_res’ folder.
#' @param workflow Character, the combinations of data processing methods specified by users according to their research interests.
#'    <br>It is a vector includes one or more method combinations, typically in the format of ‘compensation method name_ transformation method name_ normalization method name_ signal clean method name ‘, for example: c(‘None_Biexponential Transformation_None_None’,’CytoSpill_FlowVS Transformation_None_FlowCut’).
#' @param savepath Character, the absolute path of the folder which will store the results.
#' @param TIM Character, the method of trajectory inference for the processed data prior to performance assessment, consisted of tra-jectory reconstruction and data space representation, including ‘scorpius_distSpear’, ‘scorpius_distPear’, ‘scorpius_distEucl’, ‘scorpius_distManh’, ‘slingshot_tSNE’, ‘slingshot_FLOWMAP’, ‘slingshot_PCA’, ‘slingshot_diffMaps’.
#' @param clustering.var Character, the vector naming channels to be used to calculate distances/differences between cells for clustering (if re-quested) and edge-drawing steps.
#' @param plot Logical, the logical flag determining whether the visualization plot representing the results of the selected trajectory inference method should be generated and saved as a PNG file.
#' @return An RData file named "workflow_PTIres.RData", recording the PTI results. Optionally a PNG file in addition.
#' @export
#'
#' @examples
#' \donttest{
#' }

Get_PTIres <- function(
    respath, save_processed_res ="one_folder", workflow, savepath,
    TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh",
            "slingshot_tSNE","slingshot_FLOWMAP","slingshot_PCA", "slingshot_diffMaps"),
    clustering.var = NULL,plot =c(T, F)
){
  try(source("./PTI/load_data2.R"))
  try(source("./PTI/TI_method.R"))
  try(source("./PTI/ANPELA_FLOWMAP.R"))
  try(source("./PTI/ANPELA_FLOWMAP-function.R"))
  try(source("./PTI/plot.R"))

  PTIres_list <- list()

  # data & info_saved
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

  # TIM
  if (missing(TIM)) {
    TIM <- "scorpius_distSpear"
  } else {
    TIM <- match.arg(TIM)
  }

  for ( i in 1:length(workflow)){
    # AP2_processed_D_TI
    dataset_name <- workflow[i]
    if (save_processed_res == "one_folder"){
      index <- stringr::str_replace_all(info_saved$index_TIclass, "\\(.*", "")
      load(grep(dataset_name, datapath, value = T))
      if ("condition" %in% colnames(info_saved$metadata)){
        res <- try(load.Data(res, index = index, measurement.condition = as.matrix(info_saved$metadata$condition),
                             measurement.time = as.matrix(info_saved$metadata$timepoint), TIM = TIM), silent = T)
      } else {
        res <- try(load.Data(res, index = index, measurement.time = as.matrix(info_saved$metadata$timepoint), TIM = TIM), silent = T)
      }
    } else if (save_processed_res == "one_RData") {
      index <- stringr::str_replace_all(data$index_TIclass, "\\(.*", "")
      res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      if ("condition" %in% colnames(data$metadata)){
        res <- try(load.Data(res, index = index, measurement.condition = as.matrix(data$metadata$condition),
                             measurement.time = as.matrix(data$metadata$timepoint), TIM = TIM), silent = T)
      } else {
        res <- try(load.Data(res, index = index, measurement.time = as.matrix(data$metadata$timepoint), TIM = TIM), silent = T)
      }
    }

    if (class(res) == "try-error") {
      print(paste0("Can not load data for dataset: ", dataset_name))
      next
    }

    AP2_processed_D_TI <- res
    rm(index, res)

    # TIres
    TIres <- try(TI(D = AP2_processed_D_TI, method = TIM, dataset_name = dataset_name,clustering.var = clustering.var))

    PTIres_list[[dataset_name]] <- TIres
    save(PTIres_list, file = paste0(savepath, "/", dataset_name, "_PTIres.RData"))
    if (plot){
      errorplot <- ggplot() +
        annotate(geom = "text", x = 0.5, y = 0.5,
                 label = paste("During the plotting process,\n",
                               "the program had some problems and will not be able to display the output."),
                 cex = 5, color = "black") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.border = element_blank())

      if(grepl("scorpius", TIM)){
        TI_plot <- try(dr(TIres, AP2_processed_D_TI), silent = T)
        if (any(class(TI_plot) == "try-error")) {
          TI_plot <- errorplot
        }
        grDevices::png(paste0(savepath, "/", dataset_name,"_PTI.png"),bg = "white",width = 10, height = 10,res=300, units ="in")
        print(TI_plot)
        dev.off()
      } else if (grepl("slingshot", TIM)){
        if (grepl("FLOWMAP", TIM)){
          timepoint <- TIres[["timepoint"]]
        } else {
          timepoint <- AP2_processed_D_TI$timepoint
        }
        colors <- colorRampPalette(c("#eef4ed","#97c8c5","#4661a5","#183f7f"))(length(unique(timepoint)))
        grDevices::png(paste0(savepath, "/", dataset_name,"_PTI.png"),bg = "white",width = 10, height = 10,res=300, units ="in")
        plot(TIres[["dimRed"]], col = colors[as.factor(timepoint)],
             pch=16, cex = 1.5,#点的形状和大小
             asp = 1,axes = T,xlab = "reduced dimension 1", ylab = "reduced dimension 2")
        lines(slingshot::SlingshotDataSet(TIres[["crv1"]]), lwd=6, col="grey40")
        legend("topright",
               legend = levels(as.factor(timepoint)),
               col = colors,
               inset=0.8,
               pch = 16)
        if(class(TIres) == "try-error"){
          print(errorplot)
        }
        dev.off()
      }
    }
  }
}
