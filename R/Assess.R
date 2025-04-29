#'
#' @title Assess
#' @description Assess() assesses processing performance of all workflows which are used while running the function "Process", "FCprocess" or "MCprocess" based on comprehensive criteria (each with distinct underlying theories) from the perspective of CSI or PTI studies.
#' @param name Character, the filename of the RData file in the "assess_res" folder which will store the assessment results.
#' @param data Character, the R object resulting from the function "Process", "FCprocess" or "MCprocess", or obtained by loading from the resulting RData file of these funcitons when the `save_processed_res` parameter in these functions is set to "one_RData".
#' @param respath Character, the absolute path of the folder storing the resulting "info_saved.RData" file and the "process_res" folder of the function "Process", "FCprocess" or "MCprocess" when the `save_processed_res` parameter in these functions is set to "one_folder".
#' @param studytype Character, the type of your study, including "CSI (Cell Subpopulation Identification)" and "PTI (Pseudotime Trajectory Inference)".
#' @param clusteringM Character, the method of clustering the processed data prior to performance assessment, including "FlowSOM" and "PhenoGraph".
#'   <br>**FlowSOM**: a widely used cluster clustering algorithm designed for high-dimensional cytometry data (the number of clusters needs to be specified). ANPELA uses the function “SOM” in R package "FlowSOM" to implement the algorithm.
#'   <br>**PhenoGraph**: a well-designed clustering algorithm developed to define phenotypes in high-dimensional single-cell data (the number of clusters need not be specified). ANPELA uses the function “Rphenograph” in R package "cytofkit".
#' @param ncluster Integer, the number of clusters for meta clustering in FlowSOM.
#'   <br>Only needed when the argument of "clusteringM" is selected as "FlowSOM".
#' @param Phenograph_k Character, the number of nearest neighbours used in PhenoGraph clustering method.
#'   <br>Only needed when the argument of "clusteringM" is selected as "PhenoGraph".
#' @param DEP Character, the absolute filepath of the CSV file including the differentially expressed proteins used as the prior knowledge for the fourth criterion.
#'   <br>It is a table of one column without the column name, each table cell includes one protein typically in the format of "channel description (channel name)", for example: "CD20(FITC.A)".
#' @param TIM Character, the method of trajectory inference for the processed data prior to performance assessment, consisted of tra-jectory reconstruction and data space representation, including ‘scorpius_distSpear’, ‘scorpius_distPear’, ‘scorpius_distEucl’, ‘scorpius_distManh’, ‘slingshot_tSNE’, ‘slingshot_FLOWMAP’, ‘slingshot_PCA’, ‘slingshot_diffMaps’.
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

Assess <- function(
    name = "result", data = NULL, respath =NULL,
    studytype = c("CSI", "PTI"),
    clusteringM = "FlowSOM",
    ncluster = 8,
    Phenograph_k = 30,
    DEP = NULL,

    TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh",
            "slingshot_tSNE", "slingshot_FLOWMAP", "slingshot_PCA", "slingshot_diffMaps"),
    pathwayhierarchy = NULL,
    clustering.var = NULL,
    cores = floor(parallel::detectCores()/2),
    save_processed_res = "one_folder",
    savepath = "./ANPELA_res"
){
  if (studytype == "CSI"){
    assess_res <- CSIassess(name = name, data = data, respath = respath, clusteringM = clusteringM,
                            Phenograph_k = Phenograph_k, ncluster = ncluster,
                            DEP = DEP, cores = cores,
                            save_processed_res = save_processed_res, savepath = savepath)


  } else if (studytype == "PTI"){
    assess_res <- PTIassess(name = name, data = data, respath = respath, TIM = TIM,
                            pathwayhierarchy = pathwayhierarchy, cores = cores,clustering.var = clustering.var,
                            save_processed_res = save_processed_res, savepath = savepath)
  }

  return(assess_res)
}
