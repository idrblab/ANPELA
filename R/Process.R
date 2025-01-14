#' @title Process
#' @description Process() enables high-throughput processing for SCP data acquired from FC or MC by the most available workflows respectively, based on parallel computing (each workflow is distinct by combining different methods of compensation, transformation, normalization and signal clean), which facilitates the subsequent application of performance assessment, ranking and plotting.
#' @param name Character, the filename of the RData file when the "save_processed_res" parameter is set to "one_RData".
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#' @param technique Character, the technique type used in acquiring the SCP data.
#' @param studytype Character, the type of your study, including "CSI (Cell Subpopulation Identification)" and "PTI (Pseudotime Trajectory Inference)".
#' @param mergeM Character, the method of merging multiple FCS files. When multiple FCS files are selected, cells can be combined using one of the four different methods including "Fixed", "Ceil", "All" and "Min".
#'   <br>**Fixed**: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each FCS file and combined for analysis.
#'   <br>**Ceil**: up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each FCS file and combined for analysis.
#'   <br>**All**: all cells from each FCS file are combined for analysis.
#'   <br>**Min**: The minimum number of cells among all the selected FCS files are sampled from each FCS file and combined for analysis.
#' @param fixedNum Integer, the fixed number of cells to be extracted from each FCS file.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param compensationM Character, the method(s) of compensation for flow cytometry data including "AutoSpill", "CATALYST", "CytoSpill", "FlowCore", "MetaCyto" , "spillR" and "None". Compensation refers to the processing step of removing unwanted spillover resulting from signal crosstalk and spectral overlap across detection channels.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param transformationM Character, the method(s) of transformation for flow cytometry data including "Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value", "Biexponential Transformation", "Box-Cox Transformation", "Centered Log Ratio Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation", "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation" and "None". Transformation refers to the processing step of adjusting the data with a heavily skewed distribution to a normal distribution.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param normalizationM Character, the method(s) of normalization for flow cytometry data, including "Bead-based Normalization", "GaussNorm", "Mean Normalization", "Min-max Normalization", "WarpSet", "ZScore" and "None". Normalization refers to the processing step of eliminating signal decay and technical variability across all files and batches over long-term data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param signalcleanM Character, the method(s) of signal clean for flow cytometry data, including "FlowAI", "FlowClean", "FlowCut", "PeacoQC" and "None". Signal cleaning refers to the processing step of identifying and removing abrupt signal shifts and changes that derive from (i) abrupt changes in the flow rate, (ii) clogs within the capillary tubes, (iii) temporary disruptions in cytometer fluidics, and (iv) unstable data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param workflow Character, the combinations of data processing methods specified by users according to their research interests.
#'   <br>It is a vector includes one or more method combinations, typically in the format of "compensation method name_ transformation method name_ normalization method name_ signal clean method name ", for example: c("None_Biexponential Transformation_None_None","CytoSpill_FlowVS Transformation_None_FlowCut").
#' @param spillpath Character, the absolute filepath(s) of compensation beads or cells. The spillover information for a particular experiment is often obtained by running several tubes of beads or cells stained with a single color that can then be used to determine a spillover matrix for use.
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM".The filenames of the FCS files must correspond to the names of stain channels. If the original FCS files contain a pre-calculated spillover matrix as the value of the $SPILLOVER, $spillover or $SPILL keywords, this can be set as NULL.
#' @param FSC Character, the name of the forward scatter parameter.
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM".
#' @param SSC Character, the name of the side scatter parameter.
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM".
#' @param control.dir Character, the absolute path of the folder storing the FCS files of single-color controls.
#'   <br>Only needed when "AutoSpill" is included in the argument of "compensationM".
#' @param control.def.file Character, the absolute filepath of the CSV file defining the filenames and corresponding channels of the single-color controls.
#'   <br>Only needed when "AutoSpill" is included in the argument of "compensationM".
#' @param single_pos_fcs Character, the absolute filepath of the FCS file containing stained samples and control antibody-capture beads/pooled single-stained beads.
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param single_pos_mass Integer, the masses corresponding to barcode channels.
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param CATALYSTM Character, the method for solving linear system, including "flow" and "nnls".
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param beads_mass Integer, the masses of the corresponding calibration beads.
#'   <br>Only needed when "Bead-based Normalization" is included in the argument of "normalizationM".
#' @param sce_bead SingleCellExperiment, the SingleCellExperiment object for the bead experiment.
#'   <br>Only needed when "spillR" is included in the argument of "compensationM".
#' @param marker_to_barc Data frame, the table that maps the marker to the barcode in the beads experiment.
#'   <br>Only needed when "spillR" is included in the argument of "compensationM".
#' @param min_cells Integer, the minimum amount of cells (nonzero values) that should be present in one bin.
#'   <br>Only needed when "PeacoQC" is included in the argument of "signalcleanM". Lowering this parameter can affect the robustness of the peak detection.
#' @param max_bins Integer, the maximum number of bins that can be used in the cleaning process.
#'   <br>Only needed when "PeacoQC" is included in the argument of "signalcleanM". If this value is lowered, larger bins will be made.
#' @param step Integer, the step in events_per_bin to which the parameter is reduced to.
#'   <br>Only needed when "PeacoQC" is included in the argument of "signalcleanM".
#' @param excludedColumn Character, the non-protein columns names of which accessed through the function "Getmarker".
#'   <br>It is a string separated by commas, typically in the format of "channel description (channel name)", for example: "gate_source(gate_source), cell_id(cell_id), sample_id(sample_id), Time(Time), Cell_length(Cell_length), DNA-1(DNA.1.Ir191.Dd)".
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param save_processed_res Character, the format of the data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
#' @param savepath Character, the absolute path of the folder which will store the processed results.
#' @return The **process_res** folder stores the results of various data processing workflows. The form of data processing output files is decided by the parameter `save_processed_res`: "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder; "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder \[*default ="one_folder"*\].
#'   <br>In addition, the file **info_saved.RData** is also generated simultaneously, recording the information related to "metadata" and "index_protein".

#' @export
#'
#' @examples
#' \donttest{
#' }

Process <- function(
    name = "result",
    datapath,
    technique = c("MC", "FC"),
    studytype = c("CSI", "PTI"),
    mergeM = c("Fixed", "Ceil", "All", "Min"),
    fixedNum = 200,
    compensationM = c("AutoSpill", "CATALYST", "CytoSpill", "FlowCore", "MetaCyto", "spillR", "None"),
    transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                        "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                        "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation",
                        "Centered Log Ratio Transformation","None"),
    normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None"),
    signalcleanM = c("FlowAI", "FlowClean", "FlowCut", "PeacoQC", "None"),
    workflow = NULL,
    spillpath = NULL, FSC = "FSC-H", SSC = "SSC-H",
    control.dir = NULL, control.def.file = NULL,
    single_pos_fcs = NULL, single_pos_mass = NULL, CATALYSTM = c("flow", "nnls"),
    beads_mass = c(140, 151, 153, 165, 175),
    sce_bead = NULL, marker_to_barc = NULL,
    min_cells = 3, max_bins = 10, step = 10,
    excludedColumn = NULL,
    save_processed_res = "one_folder",
    savepath = "./",
    cores = floor(parallel::detectCores()/2)
){
  metadata <- paste0(datapath, "/metadata.csv")
  if (technique == "MC"){
    MCprocess_res <- MCprocess(name = name,
                               datapath = datapath,
                               metadata = metadata,
                               studytype = studytype, mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               workflow = workflow,
                               single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass, CATALYSTM = "nnls",
                               beads_mass = beads_mass,
                               sce_bead = sce_bead, marker_to_barc = marker_to_barc,
                               min_cells = min_cells, max_bins = max_bins, step =step,
                               excludedColumn = excludedColumn,
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)

    return(MCprocess_res)
  } else if (technique == "FC"){
    FCprocess_res <- FCprocess(name = name,
                               datapath = datapath,
                               metadata = metadata,
                               studytype = studytype, mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               workflow = workflow,
                               spillpath = spillpath, FSC = FSC, SSC = SSC,
                               control.dir = control.dir, control.def.file = control.def.file,
                               min_cells = min_cells, max_bins = max_bins, step =step,
                               excludedColumn = excludedColumn,
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)
    return(FCprocess_res)
  }
}
