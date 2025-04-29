
#' @title Flow Cytometry Cell Subpopulation Identification
#' @description FC_CSI() enables high-throughput processing of SCP data acquired from FC using up to ~960 available workflows and subsequently assesses the performance of all processing workflows based on comprehensive criteria from the perspective of CSI studies, functioning similarly to a combination of FCprocess() and CSIassess(), while utilizing less memory during execution.
#'
#' @param name Character, the filename of the RData file when the "save_processed_res" parameter is set to "one_RData", and the filename of all the files in the "assess_res" folder which will store the assessment results.
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#' @param technique Character, the technique type used in acquiring the SCP data.
#' @param mergeM Character, the method of merging multiple FCS files. When multiple FCS files are selected, cells can be combined using one of the four different methods including "Fixed", "Ceil", "All" and "Min".
#'   <br>**Fixed**: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each FCS file and combined for analysis.
#'   <br>**Ceil**: up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each FCS file and combined for analysis.
#'   <br>**All**: all cells from each FCS file are combined for analysis.
#'   <br>**Min**: The minimum number of cells among all the selected FCS files are sampled from each FCS file and combined for analysis.
#' @param fixedNum Integer, the fixed number of cells to be extracted from each FCS file.
#' @param compensationM Character, the method(s) of compensation for flow cytometry data including "AutoSpill", "FlowCore", "MetaCyto" and "None". Compensation refers to the processing step of removing unwanted spillover resulting from signal crosstalk and spectral overlap across detection channels.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param transformationM Character, the method(s) of transformation for flow cytometry data including "Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value", "Biexponential Transformation", "Box-Cox Transformation", "Centered Log Ratio Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation", "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation"and "None". Transformation refers to the processing step of adjusting the data with a heavily skewed distribution to a normal distribution.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param normalizationM Character, the method(s) of normalization for flow cytometry data, including "GaussNorm", "Mean Normalization", "Min-max Normalization", "WarpSet", "ZScore" and "None". Normalization refers to the processing step of eliminating signal decay and technical variability across all files and batches over long-term data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param signalcleanM Character, the method(s) of signal clean for flow cytometry data, including "FlowAI", "FlowClean", "FlowCut", "PeacoQC" and "None". Signal cleaning refers to the processing step of identifying and removing abrupt signal shifts and changes that derive from (i) abrupt changes in the flow rate, (ii) clogs within the capillary tubes, (iii) temporary disruptions in cytometer fluidics, and (iv) unstable data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param workflow Character, the combinations of data processing methods specified by users according to their research interests.
#'   <br>It is a vector includes one or more method combinations, typically in the format of "compensation method name_ transformation method name_ normalization method name_ signal clean method name ", for example: c("None_Biexponential Transformation_None_None","CytoSpill_FlowVS Transformation_None_FlowCut").
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
#' @param arcsinha Double, the argument offset coefficient ‘a’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Arcsinh Transformation’ is included in the parameter of ‘transformationM’.
#' @param arcsinhb Double, the input scaling coefficient ‘b’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Arcsinh Transformation’ is included in the parameter of ‘transformationM’.
#' @param arcsinhc Double, the output offset coefficient ‘c’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Arcsinh Transformation’ is included in the parameter of ‘transformationM’.
#' @param anna Double, the argument offset coefficient ‘a’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Asinh with Non-negative Value’ is included in the parameter of ‘transformationM’.
#' @param annb Double, the input scaling coefficient ‘b’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Asinh with Non-negative Value’ is included in the parameter of ‘transformationM’.
#' @param annc Double, the output offset coefficient ‘c’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Asinh with Non-negative Value’ is included in the parameter of ‘transformationM’.
#' @param annthreshold Double, the input cutoff value ‘threshold’ below which the input x is replaced by the threshold itself before the asinh calculation.
#'   <br>Only needed when ‘Asinh with Non-negative Value’ is included in the parameter of ‘transformationM’.
#' @param arna Double, the argument offset coefficient ‘a’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Asinh with Randomized Negative Value’ is included in the parameter of ‘transformationM’.
#' @param arnb Double, the input scaling coefficient ‘b’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Asinh with Randomized Negative Value’ is included in the parameter of ‘transformationM’.
#' @param arnc Double, the output offset coefficient ‘c’ in equation y = asinh(a + b*x) + c.
#'   <br>Only needed when ‘Asinh with Randomized Negative Value’ is included in the parameter of ‘transformationM’.
#' @param arnthreshold Double, the input cutoff value ‘threshold’ below which the input x is replaced by a small random value before the asinh calculation.
#'   <br>Only needed when ‘Asinh with Randomized Negative Value’ is included in the parameter of ‘transformationM’.
#' @param bepa Double, the positive exponential scaling coefficient ‘a’ in equation y = a*exp(b*(x-w)) - c*exp(-d*(x-w)) + f.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param bepb Double, the positive exponential rate coefficient ‘b’ in equation y = a*exp(b*(x-w)) - c*exp(-d*(x-w)) + f.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param bepc Double, the negative exponential scaling coefficient ‘c’ in equation y = a*exp(b*(x-w)) - c*exp(-d*(x-w)) + f.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param bepd Double, the negative exponential rate coefficient ‘d’ in equation y = a*exp(b*(x-w)) - c*exp(-d*(x-w)) + f.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param bepf Double, the vertical offset coefficient ‘f’ in equation y = a*exp(b*(x-w)) - c*exp(-d*(x-w)) + f.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param bepw Double, the horizontal shift coefficient ‘w’ defining the center point in equation y = a*exp(b*(x-w)) - c*exp(-d*(x-w)) + f.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param tol Double, the numerical tolerance value ‘tol’ used by the root-finding algorithm during the inversion of the function.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param maxit Integer, the maximum iterations value ‘maxit’ allowed for the root-finding algorithm during the inversion of the function.
#'   <br>Only needed when ‘Biexponential Transformation’ is included in the parameter of ‘transformationM’.
#' @param hpla Double, the scaling parameter ‘a’ determining the overall compression level and transition characteristics.
#'   <br>Only needed when ‘Hyperlog Transformation’ is included in the parameter of ‘transformationM’.
#' @param hplb Double, the linear coefficient ‘b’ controlling the width of the linear region near zero.
#'   <br>Only needed when ‘Hyperlog Transformation’ is included in the parameter of ‘transformationM’.
#' @param lntr Double, the numerator scaling coefficient ‘r’ in equation y = log(x) * (r / d).
#'   <br>Only needed when ‘Ln Transformation’ is included in the parameter of ‘transformationM’.
#' @param lntd Double, the denominator scaling coefficient ‘d’ in equation y = log(x) * (r / d).
#'   <br>Only needed when ‘Ln Transformation’ is included in the parameter of ‘transformationM’.
#' @param logbase Integer, the base of the Log Transformation.
#'   <br>Only needed when "Log Transformation" is included in the argument of "transformationM".
#' @param logr Double, the numerator scaling coefficient ‘r’ in equation y = log(x, logbase) * (r / d).
#'   <br>Only needed when ‘Log Transformation’ is included in the parameter of ‘transformationM’.
#' @param logd Double, the denominator scaling coefficient ‘d’ in equation y = log(x, logbase) * (r / d).
#'   <br> Only needed when ‘Log Transformation’ is included in the parameter of ‘transformationM’.
#' @param lgtw Double, the linear region width value ‘w’ defining the scale behavior near zero.
#'   <br> Only needed when ‘Logicle Transformation’ is included in the parameter of ‘transformationM’.
#' @param lgtt Double, the top-of-scale value ‘t’ representing the maximum expected input data value.
#'   <br> Only needed when ‘Logicle Transformation’ is included in the parameter of ‘transformationM’.
#' @param lgtm Double, the total display range ‘m’ setting the overall width of the transformed output scale.
#'   <br> Only needed when ‘Logicle Transformation’ is included in the parameter of ‘transformationM’.
#' @param lgta Double, the additional negative range ‘a’ controlling the extent of negative input values included in the display.
#'   <br> Only needed when ‘Logicle Transformation’ is included in the parameter of ‘transformationM’.
#' @param Quadratica Double, the quadratic coefficient "a" in equation y = a&#42;x^2+b&#42;x+c.
#'   <br>Only needed when "QuadraticTransform" is included in the argument of "transformationM".
#' @param Quadraticb Double, the linear coefficient "b" in equation y = a&#42;x^2+b&#42;x+c.
#'   <br>Only needed when "QuadraticTransform" is included in the argument of "transformationM".
#' @param Quadraticc Double, the intercept "c" in equation y = a&#42;x^2+b&#42;x+c.
#'   <br>Only needed when "QuadraticTransform" is included in the argument of "transformationM".
#' @param lineara Double, the multiplicative factor "a" in equation y = a&#42;x+b.
#'   <br>Only needed when "Linear Transformation" is included in the argument of "transformationM".
#' @param linearb Double, the additive factor "b" in equation y = a&#42;x+b.
#'   <br>Only needed when "Linear Transformation" is included in the argument of "transformationM".
#' @param Truncatea Double, the value at which to truncate.
#'   <br>Only needed when "TruncateTransform" is included in the argument of "transformationM".
#' @param Segment Integer, the value specifying the number of events in each segment to be analyzed.
#'   <br>Only needed when ‘FlowClean is included in the parameter of ‘signalcleanM’.
#' @param Segment2 Integer, the value representing the minimum number of cells a population must have to be included in analysis.
#'   <br>Only needed when ‘FlowClean is included in the parameter of ‘signalcleanM’.
#' @param min_cells Integer, the minimum amount of cells (nonzero values) that should be present in one bin.
#'   <br>Only needed when "PeacoQC" is included in the argument of "signalcleanM". Lowering this parameter can affect the robustness of the peak detection.
#' @param max_bins Integer, the maximum number of bins that can be used in the cleaning process.
#'   <br>Only needed when "PeacoQC" is included in the argument of "signalcleanM". If this value is lowered, larger bins will be made.
#' @param step Integer, the step in events_per_bin to which the parameter is reduced to.
#'   <br>Only needed when "PeacoQC" is included in the argument of "signalcleanM".
#' @param excludedColumn Character, the non-protein columns names of which accessed through the function "Getmarker".
#'   <br>It is a string separated by commas, typically in the format of "channel description (channel name)", for example: "gate_source(gate_source), cell_id(cell_id), sample_id(sample_id), Time(Time), Cell_length(Cell_length), DNA-1(DNA.1.Ir191.Dd)".
#' @param DEP Character, the absolute filepath of the CSV file including the differentially expressed proteins used as the prior knowledge for the fourth criterion.
#'   <br>It is a table of one column without the column name, each table cell includes one protein typically in the format of "channel description (channel name)", for example: "CD20(FITC.A)".
#' @param clusteringM Character, the method of clustering the processed data prior to plotting, including "FlowSOM" and "PhenoGraph".
#' @param ncluster Integer, the number of clusters for meta clustering in "FlowSOM".
#'   <br>Only needed when the argument of "clusteringM" is selected as "FlowSOM" while calling the function "FCprocess" or "MCprocess" for obtaining "data".
#' @param Phenograph_k Integer, the number of nearest neighbours used in PhenoGraph clustering method.
#'   <br>Only needed when the argument of "clusteringM" is selected as "PhenoGraph".
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param save_processed_res Character, the format of the data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
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

FC_CSI <- function(name = "result", datapath, technique = "FC",
                   mergeM = "Fixed", fixedNum = 200,
                   compensationM = c("AutoSpill", "FlowCore", "MetaCyto", "None"),
                   transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                                       "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                                       "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation",
                                       "Centered Log Ratio Transformation", "None"),
                   normalizationM = c("GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None"),
                   signalcleanM = c("FlowAI", "FlowClean", "FlowCut", "PeacoQC", "None"),
                   workflow = NULL,
                   spillpath = NULL, spillname = NULL, FSC = "FSC-H", SSC = "SSC-H",
                   control.dir = NULL, control.def.file = NULL,
                   arcsinha = 0, arcsinhb = 1/150, arcsinhc = 0,
                   anna = 0, annb = 1/150, annc = 0, annthreshold = 1,
                   arna = 0, arnb = 1/150, arnc = 0, arnthreshold = 1,
                   bepa = 0.5, bepb = 1, bepc = 0.5, bepd = 1, bepf = 0, bepw = 0, tol = .Machine$double.eps^0.25, maxit = as.integer(5000),
                   hpla = 1, hplb = 1,
                   lntr = 1, lntd = 1,
                   logbase = 10,logr = 1,logd = 1,
                   lgtw = 0.5, lgtt = 262144, lgtm = 4.5, lgta = 0,
                   Quadratica = 1, Quadraticb = 1, Quadraticc = 0,
                   lineara = 2, linearb = 0,
                   Truncatea = 1,
                   Segment = 200,
                   Segment2 = 200,
                   min_cells = 3, max_bins = 10, step = 10,
                   excludedColumn = NULL,

                   DEP = NULL,
                   clusteringM = "FlowSOM",
                   ncluster = 8,
                   Phenograph_k = 30,

                   cores = parallel::detectCores()/2,
                   save_processed_res = "one_folder",
                   savepath = "./ANPELA_res"
) {

  metadata <- paste0(datapath, "/metadata.csv")
  if (save_processed_res == "one_RData") {
    FCprocess_res <- FCprocess(name = name,
                               datapath = datapath,
                               metadata = metadata, technique = technique,
                               studytype = "CSI", mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               workflow = workflow,
                               spillpath = spillpath, spillname = spillname, FSC = FSC,  SSC = SSC,
                               control.dir = control.dir, control.def.file = control.def.file,
                               arcsinha, arcsinhb = arcsinhb, arcsinhc,
                               anna = anna, annb = annb, annc = annc, annthreshold = annthreshold,
                               arna = anna, arnb = arnb, arnc = arnc, arnthreshold = annthreshold,
                               bepa = bepa, bepb = bepb, bepc = bepc, bepd = bepd, bepf = bepf, bepw = bepw, tol = tol, maxit = maxit,
                               hpla = hpla, hplb = hplb,
                               lineara = lineara, linearb = linearb,
                               lntr = lntr, lntd = lntd,
                               logbase = logbase,logr = logr,logd = logd,
                               lgtw = lgtw, lgtt = lgtt, lgtm = lgtm, lgta = lgta,
                               Quadratica = Quadratica, Quadraticb = Quadraticb, Quadraticc = Quadraticc,
                               Truncatea = Truncatea,
                               Segment = Segment, Segment2 = Segment2,
                               min_cells = min_cells, max_bins = max_bins, step = step,
                               excludedColumn = excludedColumn,
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)
    FCprocess_assess_res <- CSIassess(name = name, data = FCprocess_res, clusteringM = clusteringM,
                                      Phenograph_k = Phenograph_k, ncluster = ncluster, DEP = DEP,
                                      save_processed_res = save_processed_res, savepath = savepath,
                                      cores = cores)
  } else if (save_processed_res == "one_folder") {
    FCprocess_assess_res <- oneStep_process_assess(
      name = name,
      datapath = datapath,
      metadata = metadata,
      technique = technique,
      studytype = "CSI", mergeM = mergeM, fixedNum = fixedNum,
      compensationM = compensationM,
      transformationM = transformationM,
      normalizationM = normalizationM,
      signalcleanM = signalcleanM,
      workflow = workflow,
      spillpath = spillpath, spillname = spillname, FSC = FSC,  SSC = SSC,
      control.dir = control.dir, control.def.file = control.def.file,
      arcsinha, arcsinhb = arcsinhb, arcsinhc,
      anna = anna, annb = annb, annc = annc, annthreshold = annthreshold,
      arna = anna, arnb = arnb, arnc = arnc, arnthreshold = annthreshold,
      bepa = bepa, bepb = bepb, bepc = bepc, bepd = bepd, bepf = bepf, bepw = bepw, tol = tol, maxit = maxit,
      hpla = hpla, hplb = hplb,
      lineara = lineara, linearb = linearb,
      lntr = lntr, lntd = lntd,
      logbase = logbase,logr = logr,logd = logd,
      lgtw = lgtw, lgtt = lgtt, lgtm = lgtm, lgta = lgta,
      Quadratica = Quadratica, Quadraticb = Quadraticb, Quadraticc = Quadraticc,
      Truncatea = Truncatea,
      Segment = Segment, Segment2 = Segment2,
      min_cells = min_cells, max_bins = max_bins, step = step,
      excludedColumn = excludedColumn,

      clusteringM = clusteringM, Phenograph_k = Phenograph_k,
      ncluster = ncluster, DEP = DEP,

      cores = cores,
      save_processed_res = save_processed_res,
      savepath = savepath
    )
  }

  Ranking(FCprocess_assess_res, name = name, savepath = paste0(savepath, "/assess_res"))
}
