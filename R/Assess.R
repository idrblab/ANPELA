#'
#' @title Assess
#' @description
#' A short description...
#'
#' @param name Character, the filename of all the files in the "assess_res" folder which will store the assessment results.
#' @param data Character, the resulting RData file of the function "Process", "FCprocess" or "MCprocess" when the "save_processed_res" parameter in these functions is set to "one_RData".
#' @param respath Character, the absolute path of the folder storing the resulting RData files of the function "Process", "FCprocess" or "MCprocess" when the "save_processed_res" parameter in these functions is set to "one_folder".
#' @param studytype Character, the type of your study, including "CSI (Cell Subpopulation Identification)" and "PTI (Pseudotime Trajectory Inference)".
#' @param clusteringM Character, the method of clustering the processed data prior to performance assessment, including "FlowSOM" and "PhenoGraph".
#'   <br>**FlowSOM**: a widely used cluster clustering algorithm designed for high-dimensional cytometry data (the number of clusters needs to be specified). ANPELA uses the function “SOM” in R package "FlowSOM" to implement the algorithm.
#'   <br>**PhenoGraph**: a well-designed clustering algorithm developed to define phenotypes in high-dimensional single-cell data (the number of clusters need not be specified). ANPELA uses the function “Rphenograph” in R package "cytofkit".
#' @param ncluster Integer, the number of clusters for meta clustering in FlowSOM.
#'   <br>Only needed when the argument of "clusteringM" is selected as "FlowSOM".
#' @param Phenograph_k Character, the number of nearest neighbours used in PhenoGraph clustering method.
#'   <br>Only needed when the argument of "clusteringM" is selected as "PhenoGraph".
#' @param DEP Character, the differentially expressed proteins used as the prior knowledge for the fourth criterion.
#' @param TIM Character, the method of trajectory inference for the processed data prior to performance assessment, consisted of trajectory reconstruction and data space representation, including "scorpius_distSpear", "scorpius_distPear", "scorpius_distManh","slingshot_tSNE", "prinCurves_tSNE", "slingshot_PCA", "slingshot_diffMaps" and "prinCurves_diffMaps".
#' @param pathwayhierarchy Character, the absolute filepath of the pathway hierarchy file.
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param save_processed_res Character, the form of data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
#' @param savepath Character, the absolute path of the folder which will store files of the assessment and processed results.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' }

Assess <- function(
    name, data = NULL, respath =NULL,
    studytype = c("CSI", "PTI"),
    clusteringM = "FlowSOM",
    ncluster = 8,
    Phenograph_k = 30,
    DEP = NULL,

    TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh", "slingshot_FLOWMAP", "slingshot_tSNE",
            "prinCurves_tSNE", "slingshot_PCA", "slingshot_diffMaps", "prinCurves_diffMaps"),
    pathwayhierarchy = NULL,
    cores = floor(parallel::detectCores()/2),
    save_processed_res = c("no", "one_folder", "one_RData"),
    savepath = "./"
){
  if(studytype == "CSI"){
    assess_res <- CSIassess(data = data, respath = respath, clusteringM = clusteringM,
                            Phenograph_k = Phenograph_k, ncluster = ncluster,
                            DEP = DEP, cores = cores,
                            save_processed_res = save_processed_res)


  } else if (studytype == "PTI"){
    assess_res <- PTIassess(data = data, respath = respath, TIM = TIM, pathwayhierarchy = pathwayhierarchy,
                            cores = cores,save_processed_res = save_processed_res)
  }
  if (!dir.exists(paste0(savepath, "/assess_res"))) {
    dir.create(paste0(savepath, "/assess_res"), recursive = T)
  }
  save(assess_res, file = paste0(savepath, "/assess_res/", name, "_assess.RData"))
  Ranking(assess_res, name = name, savepath = paste0(savepath, "/assess_res"))
  return(assess_res)
}
