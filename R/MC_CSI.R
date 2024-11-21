
#' @title Mass Cytometry Cell Subpopulation Identification
#' @description MC_CSI enables high-throughput processing of SCP data acquired from MC using up to ~675 available workflows and subsequently assesses the performance of all processing workflows based on comprehensive criteria from the perspective of CSI studies, functioning similarly to a combination of MCprocess() and CSIassess(), while utilizing less memory during execution.
#'
#' @param name Character, the filename of all the files in the "assess_res" folder which will store assessment results.
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#' @param mergeM Character, the method of merging multiple FCS files. When multiple FCS files are selected, cells can be combined using one of the four different methods including "Fixed", "Ceil", "All" and "Min".
#'   <br>**Fixed**: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each FCS file and combined for analysis.
#'   <br>**Ceil**: up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each FCS file and combined for analysis.
#'   <br>**All**: all cells from each FCS file are combined for analysis.
#'   <br>**Min**: The minimum number of cells among all the selected FCS files are sampled from each FCS file and combined for analysis.
#' @param fixedNum Integer, the fixed number of cells to be extracted from each FCS file.
#' @param compensationM Character, the method(s) of compensation for mass cytometry data including "CATALYST", "CytoSpill" and "None". Compensation refers to the processing step of removing unwanted spillover resulting from signal crosstalk and spectral overlap across detection channels.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param transformationM Character, the method(s) of transformation for mass cytometry data including "Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",  "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",  "LnTransform", "Log Transformation", "Logicle Transformation", "QuadraticTransform", "ScaleTransform", "TruncateTransform" and "None". Transformation refers to the processing step of adjusting the data with a heavily skewed distribution to a normal distribution.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param normalizationM Character, the method(s) of normalization for mass cytometry data, including "Bead-based Normalization", "GaussNorm", "WarpSet" "ZScore" and "None". Normalization refers to the processing step of eliminating signal decay and technical variability across all files and batches over long-term data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param signalcleanM Character, the method(s) of signal clean for mass cytometry data, including "FlowAI",  "FlowCut" and "None". Signal cleaning refers to the processing step of identifying and removing abrupt signal shifts and changes that derive from (i) abrupt changes in the flow rate, (ii) clogs within the capillary tubes, (iii) temporary disruptions in cytometer fluidics, and (iv) unstable data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param single_pos_fcs Character, the absolute filepath of the FCS file containing stained samples and control antibody-capture beads/pooled single-stained beads.
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param single_pos_mass Integer, the masses corresponding to barcode channels.
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param beads_mass Integer, the masses of the corresponding calibration beads.
#'   <br>Only needed when "Bead-based Normalization" is included in the argument of "normalizationM".
#' @param index_protein Character, the marker indexes for data processing and performance assessment accessed through the function "Getmarker", with manual removal of non-protein columns.
#'   <br>It is a string separated by commas, typically in the format of “channel description (channel name)”, for example: “CD126(Dy161Di), CD39(Dy162Di), CD20(Dy163Di), CD161(Dy164Di)".
#' @param DEP Character, the differentially expressed proteins used as the prior knowledge for the fourth criterion.
#' @param clusteringM Character, the method of clustering the processed data prior to plotting, including "FlowSOM" and "PhenoGraph".
#' @param ncluster Integer, the number of clusters for meta clustering in "FlowSOM".
#'   <br>Only needed when the argument of "clusteringM" is selected as "FlowSOM" while calling the function "FCprocess" or "MCprocess" for obtaining "data".
#' @param Phenograph_k Integer, the number of nearest neighbours used in PhenoGraph clustering method.
#'   <br>Only needed when the argument of "clusteringM" is selected as "PhenoGraph".
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param save_processed_res Character, the form of data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
#' @param savepath Character, the absolute path of the folder which will store files of the assessment and processed results.
#'
#' @return the data processing and performance assessment output files, which are located in the **process_res** and **assess_res** folders, respectively.
#'   ## Results of Data Processing
#'   The **process_res** folder stores the results of various data processing workflows. The form of data processing output files is decided by the parameter `save_processed_res`: "no" denotes that the results would not be saved; "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder; "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder \[*default ="one_folder"*\].
#'   ## Results of Performance Assessment
#'   The **assess_res** folder stores the results of performance assessment, which includes 3 files:
#'   <br>(**1**) **_Ranking_Table.csv** contains the overall ranking results of all data processing workflows, where the "Rank" column represents the overall ranking, and the "Value" column shows the scores for different assessment criteria.
#'   <br>(**2**) **_Ranking_Figure.pdf** visualizes the overall ranking results to help users better understand the differences among various data processing workflows. Different colors represent different performance assessment levels: dark green indicates "superior," light green indicates "good," and red indicates "poor."
#'   <br>(**3**) **_assess.RData** contains 2 lists, "table" and "table2", which provide the raw scores for different assessment criteria and performance assessment levels categorized by thresholds, respectively.
#'   ## Records
#'   In addition, **log.txt** and **info_saved.RData** files are also generated simultaneously. **log.txt** records the processing details while **info_saved.RData** records the information related to "metadata" and "index_protein".
#' @export
#'
#' @examples
#' \donttest{
#' }

MC_CSI <- function(name, datapath,
                   mergeM = "Fixed", fixedNum = 200,
                   compensationM = c("CATALYST", "CytoSpill", "None"),
                   transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                                       "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                                       "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Scale Transformation", "Truncate Transformation",
                                       "None"),
                   normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "None"),
                   signalcleanM = c("FlowAI", "FlowCut", "None"),
                   single_pos_fcs = NULL, single_pos_mass = NULL, beads_mass = c(140, 151, 153, 165, 175),
                   index_protein = NULL,

                   DEP = NULL,
                   clusteringM = "FlowSOM",
                   ncluster = 8,
                   Phenograph_k = 30,

                   cores = parallel::detectCores()/2,
                   save_processed_res = "one_folder",
                   savepath = "./"
) {
  metadata <- paste0(datapath, "/metadata.csv")
  if (save_processed_res == "one_RData") {
    MCprocess_res <- MCprocess(datapath = datapath,
                               metadata = metadata,
                               studytype = "CSI", mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass, CATALYSTM = "nnls",
                               beads_mass = beads_mass,
                               index_protein = index_protein,
                               save_processed_res = save_processed_res,
                               cores = cores)
    if (!dir.exists(paste0(savepath, "/process_res"))) {
      dir.create(paste0(savepath, "/process_res"), recursive = T)
    }
    save(MCprocess_res, file = paste0(savepath, "/process_res/", name, ".RData"))
    MCprocess_assess_res <- CSIassess(MCprocess_res, clusteringM = clusteringM, Phenograph_k = Phenograph_k, ncluster = ncluster, DEP = DEP, cores = cores)
  } else if (save_processed_res %in% c("no", "one_folder")) {
    MCprocess_assess_res <- oneStep_process_assess(
      datapath = datapath,
      metadata = metadata,
      techique = "MC",
      studytype = "CSI", mergeM = mergeM, fixedNum = fixedNum,
      compensationM = compensationM,
      transformationM = transformationM,
      normalizationM = normalizationM,
      signalcleanM = signalcleanM,
      single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass, CATALYSTM = "nnls",
      beads_mass = beads_mass,
      index_protein = index_protein,

      clusteringM = clusteringM, Phenograph_k = Phenograph_k,
      ncluster = ncluster, DEP = DEP,

      cores = cores,
      save_processed_res = save_processed_res,
      savepath = savepath
    )
  }

  if (!dir.exists(paste0(savepath, "/assess_res"))) {
    dir.create(paste0(savepath, "/assess_res"), recursive = T)
  }
  save(MCprocess_assess_res, file = paste0(savepath, "/assess_res/", name, "_assess.RData"))
  # write.csv(MCprocess_assess_res$table, file = paste0(savepath, "/assess_res/", name, "_assess_table.csv"))
  Ranking(MCprocess_assess_res, name = name, savepath = paste0(savepath, "/assess_res"))
}
