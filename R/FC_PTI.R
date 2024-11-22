
#' @title Flow Cytometry Pseudotime Trajectory Inference
#' @description FC_PTI() enables high-throughput processing of SCP data acquired from FC using up to ~960 available workflows and subsequently assesses the performance of all processing workflows based on comprehensive criteria from the perspective of PTI studies, functioning similarly to a combination of FCprocess() and PTIassess(), while utilizing less memory during execution.
#'
#' @param name Character, the filename of all the files in the "assess_res" folder which will store the assessment results.
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#' @param mergeM Character, the method of merging multiple FCS files. When multiple FCS files are selected, cells can be combined using one of the four different methods including "Fixed", "Ceil", "All" and "Min".
#'   <br>**Fixed**: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each FCS file and combined for analysis.
#'   <br>**Ceil**: up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each FCS file and combined for analysis.
#'   <br>**All**: all cells from each FCS file are combined for analysis.
#'   <br>**Min**: The minimum number of cells among all the selected FCS files are sampled from each FCS file and combined for analysis.
#' @param fixedNum Integer, the fixed number of cells to be extracted from each FCS file.
#' @param compensationM Character, the method(s) of compensation for flow cytometry data including "AutoSpill", "FlowCore", "MetaCyto" and "None". Compensation refers to the processing step of removing unwanted spillover resulting from signal crosstalk and spectral overlap across detection channels.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param transformationM Character, the method(s) of transformation for flow cytometry data including "Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value", "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation", "LnTransform", "Log Transformation", "Logicle Transformation", "QuadraticTransform", "ScaleTransform", "TruncateTransform" and "None". Transformation refers to the processing step of adjusting the data with a heavily skewed distribution to a normal distribution
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param normalizationM Character, the method(s) of normalization for flow cytometry data, including "GaussNorm", "WarpSet" "ZScore" and "None". Normalization refers to the processing step of eliminating signal decay and technical variability across all files and batches over long-term data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param signalcleanM Character, the method(s) of signal clean for flow cytometry data, including "FlowAI", "FlowClean", "FlowCut" and "None". Signal cleaning refers to the processing step of identifying and removing abrupt signal shifts and changes that derive from (i) abrupt changes in the flow rate, (ii) clogs within the capillary tubes, (iii) temporary disruptions in cytometer fluidics, and (iv) unstable data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param spillpath Character, the absolute filepath(s) of compensation beads or cells. The spillover information for a particular experiment is often obtained by running several tubes of beads or cells stained with a single color that can then be used to determine a spillover matrix for use.
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM". The filenames of the FCS files must correspond to the names of stain channels. If the original FCS files contain a pre-calculated spillover matrix as the value of the $SPILLOVER, $spillover or $SPILL keywords, this can be set as NULL.
#' @param FSC Character, the name of the forward scatter parameter.
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM".
#' @param SSC Character, the name of the side scatter parameter.
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM".
#' @param control.dir Character, the absolute path of the folder storing the FCS files of single-color controls.
#'   <br>Only needed when "AutoSpill" is included in the argument of "compensationM".
#' @param control.def.file Character, the absolute filepath of the CSV file defining the filenames and corresponding channels of the single-color controls.
#'   <br>Only needed when "AutoSpill" is included in the argument of "compensationM".
#' @param index_protein Character, the marker indexes for data processing and performance assessment accessed through the function "Getmarker", with manual removal of non-protein columns.
#'   <br>It is a string separated by commas, typically in the format of “channel description (channel name)”, for example: “CD126(Dy161Di), CD39(Dy162Di), CD20(Dy163Di), CD161(Dy164Di)".
#' @param TIM Character, the method of trajectory inference for the processed data prior to performance assessment, including "scorpius_distSpear", "scorpius_distPear", "scorpius_distManh","slingshot_tSNE", "prinCurves_tSNE", "slingshot_PCA", "slingshot_diffMaps" and "prinCurves_diffMaps".
#'   <br>For details of each method, please refer to the Methods Introduction.
#' @param pathwayhierarchy Character, the absolute filepath of the pathway hierarchy file.
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
#'   <br>(**3**) **_assess.RData** contains 2 lists, "table" and "table2", providing the raw scores for different assessment criteria and performance assessment levels categorized by thresholds, respectively.
#'   ## Records
#'   In addition, **log.txt** and **info_saved.RData** files are also generated simultaneously. **log.txt** records the processing details while **info_saved.RData** records the information related to "metadata" and "index_protein".
#' @export
#'
#' @examples
#' \donttest{
#' }

FC_PTI <- function(name, datapath,
                   mergeM = "Fixed", fixedNum = 200,
                   compensationM = c("AutoSpill", "FlowCore", "MetaCyto", "None"),
                   transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                                       "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                                       "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Scale Transformation", "Truncate Transformation",
                                       "None"),
                   normalizationM = c("GaussNorm", "WarpSet", "ZScore", "None"),
                   signalcleanM = c("FlowAI", "FlowClean", "FlowCut", "None"),
                   spillpath = NULL, FSC = "FSC-H", SSC = "SSC-H", control.dir = NULL, control.def.file = NULL,
                   index_protein = NULL,

                   TIM = NULL,
                   pathwayhierarchy = NULL,

                   cores = parallel::detectCores()/2,
                   save_processed_res = "one_folder",
                   savepath = "./"
) {

  metadata <- paste0(datapath, "/metadata.csv")
  if (save_processed_res == "one_RData") {
    FCprocess_res <- FCprocess(name = name,
                               datapath = datapath,
                               metadata = metadata,
                               studytype = "PTI", mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               spillpath = spillpath, FSC = FSC, SSC = SSC,
                               control.dir = control.dir, control.def.file = control.def.file,
                               index_protein = index_protein,
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)
    FCprocess_assess_res <- PTIassess(name = name,data = FCprocess_res, TIM = TIM,
                                      pathwayhierarchy = pathwayhierarchy,
                                      save_processed_res = save_processed_res,
                                      savepath = savepath, cores = cores)
  }  else if (save_processed_res == "one_folder") {
    FCprocess_assess_res <- oneStep_process_assess(
      name = name,
      datapath = datapath,
      metadata = metadata,
      techique = "FC",
      studytype = "PTI", mergeM = mergeM, fixedNum = fixedNum,
      compensationM = compensationM,
      transformationM = transformationM,
      normalizationM = normalizationM,
      signalcleanM = signalcleanM,
      spillpath = spillpath, FSC = FSC, SSC = SSC,
      control.dir = control.dir, control.def.file = control.def.file,
      index_protein = index_protein,

      TIM = TIM, pathwayhierarchy = pathwayhierarchy,

      cores = cores,
      save_processed_res = save_processed_res,
      savepath = savepath
    )
  }

  Ranking(FCprocess_assess_res, name = name, savepath = paste0(savepath, "/assess_res"))

}
