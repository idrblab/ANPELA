#' @title Flow Cytometry Process
#' @description FCprocess() enables high-throughput processing for SCP data acquired from FC by at most ~960 available workflows based on parallel computing (each workflow is distinct by combining different methods of compensation, transformation, normalization and signal clean), which facilitates the subsequent application of performance assessment, ranking and plotting.
#'
#' @param name Character, the filename of the RData file when the "save_processed_res" parameter is set to "one_RData".
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#' @param metadata Character, the absolute filepath of the metadata file, usually located in the datapath folder. The metadata file should include the columns of "filename" and "condition" for CSI studies, and the columns of "filename" and "timepoint" for PTI studies.
#'   <br>For details on preparing the metadata file, please refer to the **sample data**.
#' @param technique Character, the technique type used in acquiring the SCP data.
#' @param studytype Character, the type of your study, including "CSI (Cell Subpopulation Identification)" and "PTI (Pseudotime Trajectory Inference)".
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
#'   <br>Only needed when "FlowCore" is included in the argument of "compensationM".The filenames of the FCS files must correspond to the names of stain channels. If the original FCS files contain a pre-calculated spillover matrix as the value of the $SPILLOVER, $spillover or $SPILL keywords, this can be set as NULL.
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
#' @param cores Integer, the number of CPU cores to be employed for performing parallel computing.
#'   <br>To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.
#' @param save_processed_res Character, the format of the data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.
#' @param savepath Character, the absolute path of the folder which will store the processed results.
#'
#' @return The **process_res** folder stores the results of various data processing workflows. The form of data processing output files is decided by the parameter `save_processed_res`: "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder; "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder \[*default ="one_folder"*\].
#'   <br>In addition, the file **info_saved.RData** is also generated simultaneously, recording the information related to "metadata" and "index_protein".
#' @export
#'
#' @examples
#' \donttest{
#' }

FCprocess <- function(name = "result",
                      datapath,
                      metadata,
                      technique = "FC",
                      studytype = c("CSI", "PTI"),
                      mergeM = c("Fixed", "Ceil", "All", "Min"),
                      fixedNum = 200,
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
                      save_processed_res = "one_folder",
                      savepath = "./ANPELA_res",
                      cores = floor(parallel::detectCores()/2), ...) {

  # dataFiles
  if (missing(datapath)) {
    stop("The parameter of 'datapath' is missing.")
  } else if (file.info(datapath)$isdir) {
    dataFiles <- list.files(datapath, pattern = ".fcs$", full.names = TRUE)
  } else if (!file.info(datapath)$isdir) {
    if (any(grepl(".fcs$", datapath))) {
      dataFiles <- datapath[grepl(".fcs$", datapath)]
    } else {
      stop("The format of parameter 'datapath' is incorrect. Please input the absolute path of the folder storing the FCS raw data files.")
    }
  }
  if (length(dataFiles) < 1) {
    stop("No FCS file found, please check the parameter of 'datapath'!")
  }


  # metadata
  if (missing(metadata)) {
    stop("The parameter of 'metadata' is missing.")
  } else if (grepl(".csv$", metadata)) {
    metadata <- read.csv(metadata)
  } else {
    stop("The format of parameter 'metadata' is incorrect. Please input the absolute filepath of the metadata file.")
  }


  # studytype
  if (missing(studytype)) {
    stop("The parameter of 'studytype' is missing.")
  } else {
    studytype <- match.arg(studytype)
  }


  # mergeM
  if (missing(mergeM)) {
    mergeM <- "Fixed"
  } else {
    mergeM <- match.arg(mergeM)
  }


  # compensationM
  if (missing(compensationM)) {
    compensationM <- c("AutoSpill", "FlowCore", "MetaCyto", "None")
  } else {
    compensationM <- match.arg(compensationM, several.ok = TRUE)
  }


  # transformationM
  if (missing(transformationM)) {
    transformationM <- c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                         "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                         "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation",
                         "Centered Log Ratio Transformation", "None")
  } else {
    transformationM <- match.arg(transformationM, several.ok = TRUE)
  }


  # normalizationM
  if (missing(normalizationM)) {
    normalizationM <- c("GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None")
  } else {
    normalizationM <- match.arg(normalizationM, several.ok = TRUE)
  }


  # signalcleanM
  if (missing(signalcleanM)) {
    signalcleanM <- c("FlowAI", "FlowClean", "FlowCut", "PeacoQC", "None")
  } else {
    signalcleanM <- match.arg(signalcleanM, several.ok = TRUE)
  }


  # control.dir
  if ("AutoSpill" %in% compensationM & is.null(workflow)|any(grepl("AutoSpill", workflow))) {
    if (is.null(control.dir)) {
      message("The parameter of 'control.dir' is missing. 'AutoSpill' compensation method can't be performed without the control files.")
      compensationM <- setdiff(compensationM, "AutoSpill")
    } else if (!file.info(control.dir)$isdir) {
      stop("The parameter of 'control.dir' is incorrect. 'AutoSpill' can't be performed without the control files. Please input the absolute path of the folder storing FCS files of single-color controls if you want to use 'AutoSpill' compensation method.")
    }
  }


  # control.def.file
  if ("AutoSpill" %in% compensationM & is.null(workflow)|any(grepl("AutoSpill", workflow))) {
    if (is.null(control.def.file)) {
      message("The parameter of 'control.def.file' is missing. 'AutoSpill' compensation method can't be performed.")
      compensationM <- setdiff(compensationM, "AutoSpill")
    } else if (!isFALSE(file.info(control.def.file)$isdir) || !grepl(".csv$", control.def.file)) {
      stop("The format of parameter 'control.def.file' is incorrect. Please input the absolute filepath of your .csv file defining the filenames and corresponding channels of the single-color controls if you want to use 'AutoSpill' compensation method..")
    }
  }


  # logbase
  if (logbase <= 0 || logbase == 1) {
    stop("The value of logbase cannot be less than or equal to 0 or equal to 1.")
  }


  # arcsinhb
  if (arcsinhb == 0) {
    stop("The value of arcsinhb cannot be 0.")
  }


  # annb
  if (annb == 0) {
    stop("The value of annb cannot be 0.")
  }


  # arnb
  if (arnb == 0) {
    stop("The value of arnb cannot be 0.")
  }


  # dataFileNames
  dataFileNames <- dataFiles


  # flag1
  flag1 <- isTRUE(all.equal(sort(limma::removeExt(basename(dataFileNames), sep=".")), sort(as.character(metadata$filename))))
  if (!flag1) {
    stop("The filenames of your FCS files are inconsistent with those of your metadata file.
         \nNote that ANPELA requires the user to input FCS files whose filename order is exactly the same as that of the metadata file.
         \nPlease input the raw FCS files & metadata in the correct format.")
  }


  # flag2
  flag2 <- switch (studytype,
                   CSI = {
                     res <- dplyr::group_by(metadata, condition) %>% dplyr::summarise(n())
                     nrow(res) == 2 && !(any((res$`n()` >= 2) == FALSE))
                   },
                   PTI = length(unique(metadata$timepoint)) >= 2
  )
  if (!flag2) {
    stop("The number of FCS file(s) is not enough for the subsequent analysis.
         \nParticularly, for two-class research, at least two samples for each class are required; for trajectory inference research, at least two time points are required.
         \nPlease refresh the page and reupload the raw FCS files & metadata in the correct format.")
  }


  # class/TI frame list
  AP2_frame_classTI <- list()
  for (i in 1:length(dataFiles)) {
    AP2_frame_classTI[[i]] <- suppressWarnings(flowCore::read.FCS(filename = dataFiles[[i]], transformation = FALSE, alter.names = TRUE))
    names(AP2_frame_classTI)[i] <- limma::removeExt(basename(dataFileNames), sep=".")[i]
  }


  # class/TI down sample expr (merged)
  exprsL <- list()
  for (i in 1:length(AP2_frame_classTI)) {
    data <- AP2_frame_classTI[[i]]@exprs
    pd <- AP2_frame_classTI[[i]]@parameters@data
    colnames(data) <- paste0(pd$desc, "(", pd$name, ")")
    rownames(data) <- paste(names(AP2_frame_classTI)[i], 1:nrow(data), sep = "_")
    exprsL[[i]] <- data
  }

  protein_list <- lapply(exprsL, function(x){
    return(colnames(x))
  })
  common_protein <- Reduce(intersect, protein_list)
  exprsL <- lapply(exprsL, function(x){
    res <- x[,common_protein]
    return(res)
  })

  set.seed(123)
  eventCountTest <- suppressWarnings(any(lapply(exprsL, function(x) if (nrow(x) < fixedNum) {
    1
  } else {
    0
  })))
  if (mergeM == "Fixed" && eventCountTest == TRUE) {
    fixedNum <- min(rapply(exprsL, nrow))
  }
  switch(mergeM, Ceil = {
    mergeFunc <- function(x) {
      if (nrow(x) < fixedNum) {
        x
      } else {
        x[sample(nrow(x), size = fixedNum, replace = FALSE), , drop = FALSE]
      }
    }
    down_exprsL <- lapply(exprsL, mergeFunc)
    merged <- do.call(rbind, down_exprsL)
  }, All = {
    down_exprsL <- exprsL
    merged <- do.call(rbind, down_exprsL)
  }, Fixed = {
    mergeFunc <- function(x) {
      x[sample(nrow(x), size = fixedNum, replace = ifelse(nrow(x) < fixedNum, TRUE, FALSE)), , drop = FALSE]
    }
    down_exprsL <- lapply(exprsL, mergeFunc)
    merged <- do.call(rbind, down_exprsL)
  }, Min = {
    minSize <- min(sapply(exprsL, nrow))
    mergeFunc <- function(x) {
      x[sample(nrow(x), size = minSize, replace = FALSE), , drop = FALSE]
    }
    down_exprsL <- lapply(exprsL, mergeFunc)
    merged <- do.call(rbind, down_exprsL)
  })
  AP2_downsample_expr_classTI <- list(down_exprsL = down_exprsL, merged = merged)
  rm(exprsL)

  common_protein1 <- gsub(".*\\(", "", common_protein)
  common_protein1 <- gsub("\\)", "", common_protein1)

  # class/TI down sample frame list
  AP2_pro0_frame <- lapply(1:length(AP2_frame_classTI), function(i, AP2_frame_classTI, AP2_downsample_expr_classTI) {
    res <- AP2_frame_classTI[[i]]
    res@exprs <- AP2_downsample_expr_classTI$down_exprsL[[i]]
    par <- flowWorkspace::pData(res@parameters)
    index <- match(common_protein1, par$name)
    flowWorkspace::pData(res@parameters) <- par[index,]
    return(res)
  }, AP2_frame_classTI = AP2_frame_classTI, AP2_downsample_expr_classTI = AP2_downsample_expr_classTI)
  names(AP2_pro0_frame) <- limma::removeExt(basename(dataFileNames), sep=".")
  rm(AP2_frame_classTI)

  # ncell_frame
  ncell_frame <- sapply(AP2_downsample_expr_classTI$down_exprsL, function(x) nrow(x))
  rm(AP2_downsample_expr_classTI)

  # spillpath
  if ("FlowCore" %in% compensationM & is.null(workflow)|any(grepl("FlowCore", workflow))) {
    if (is.null(spillpath)) {

      if (all(sapply(AP2_pro0_frame, function(x) !is.null(x@description[["SPILL"]])))) {
      } else if (all(sapply(AP2_pro0_frame, function(x) !is.null(x@description[["SPILLOVER"]])))) {
      } else if (all(sapply(AP2_pro0_frame, function(x) !is.null(x@description[["spillover"]])))) {
      } else {
        message("Some/All of your FCS files don't contain a pre-calculated spillover matrix. 'FlowCore' can't be performed without the spillover matrix.")
        compensationM <- setdiff(compensationM, "FlowCore")
      }

    } else if (file.info(spillpath)$isdir || any(!grepl(".fcs$", spillpath))) {
      stop("The format of parameter 'spillpath' is incorrect. Please input the absolute filepaths of the set of compensation controls.")
    } else {
      spillname <- gsub(".fcs", "", basename(spillpath))
    }
  }

  # MetaCyto
  if ("MetaCyto" %in% compensationM & is.null(workflow)|any(grepl("MetaCyto", workflow))) {
    if (all(sapply(AP2_pro0_frame, function(x) !is.null(x@description[["SPILL"]])))) {
    } else {
      message("Some/All of your FCS files don't contain a pre-calculated spillover matrix. 'MetaCyto' can't be performed without the spillover matrix.")
      compensationM <- setdiff(compensationM, "MetaCyto")
    }
  }


  # index_TIclass
  if (is.null(excludedColumn)) {
    cat("*************************************************************************", "\n")
    cat("The standardized column names are listed below. \n")
    cat("*************************************************************************", "\n")
    cat(paste0(common_protein, collapse = "\n"), "\n")
    cat("*************************************************************************", "\n")
    cat("Please enter the excluded column names which are separated by the comma on a single line.
          \nFor example, Time(Time), FSC-W(FSC.W), SSC-A(SSC.A)")

    excludedColumn <- readline("Now, you can select the excluded columns for data processing and performance evaluation:")
    if (excludedColumn == "") {
      stop("Note: Please select your excluded columns and then press the Enter key :(")
    }
  }
  Excluded_index_TIclass <- unlist(strsplit(excludedColumn, "\\s*,\\s*"))
  index_TIclass <- setdiff(common_protein,Excluded_index_TIclass)

  cat("The program is running. Please wait patiently.")

  # Segment, Segment2
  Segment2 <- Segment <- floor(min(ncell_frame)/3)

  # workflow
  if(is.null(workflow)){
    workflow <- expand.grid(compensation = compensationM, transformation = transformationM, normalization = normalizationM, signalclean = signalcleanM, stringsAsFactors = FALSE)
    rownames(workflow) <- paste(workflow$compensation,
                                workflow$transformation,
                                workflow$normalization,
                                workflow$signalclean, sep = "_")
  } else {
    workflow <- try(lapply(strsplit(workflow, "_"), unlist))
    if (class(workflow) == "try-error") {
      stop("The format of parameter 'workflow' is incorrect. Please input the workflow in the correct format.")
    } else {
      workflow <- as.data.frame(do.call(rbind, workflow), stringsAsFactors = FALSE)
      colnames(workflow) <- c("compensation", "transformation", "normalization", "signalclean")
      rownames(workflow) <- paste(workflow$compensation,
                                  workflow$transformation,
                                  workflow$normalization,
                                  workflow$signalclean, sep = "_")
    }
  }

  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }

  # parallel start
  opts <- list(progress = function(n) setTxtProgressBar(txtProgressBar(min = 1, max = nrow(workflow), style = 3), n))
  cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/FCprocess_log.txt"))
  doSNOW::registerDoSNOW(cl)
  time = proc.time()

  AP2_pro1_frame_classTI <- foreach::foreach(i = 1:nrow(workflow), .options.snow = opts,
                                             .packages = c("stringr", "flowStats", "flowAI", "flowCore", "flowClean", "flowCut", "flowTrans", "magrittr")) %dopar% {
                                               try(source("./processing.R"), silent = T)
                                               set.seed(123)
                                               AP2_comp_frame <- try(comp_anpela(data = AP2_pro0_frame, method = workflow[i,1], index = index_TIclass,
                                                                                 spillpath = spillpath, spillname = spillname, FSC = FSC,  SSC = SSC,
                                                                                 control.dir = control.dir, control.def.file = control.def.file), silent = T)
                                               if (class(AP2_comp_frame) == "try-error") {
                                                 # AP2_comp_frame <- AP2_pro0_frame
                                                 return(NULL)
                                               }

                                               AP2_trans_frame <- try(trans_anpela(data = AP2_comp_frame, method = workflow[i,2], index = index_TIclass,
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
                                                                                   Truncatea = Truncatea
                                               ), silent = T)
                                               if (class(AP2_trans_frame) == "try-error") {
                                                 # AP2_trans_frame <- AP2_comp_frame
                                                 return(NULL)
                                               }
                                               rm(AP2_comp_frame)
                                               proteins_excluded <- c()
                                               for (x in seq(AP2_trans_frame)) {
                                                 data_checked <- AP2_trans_frame[[x]]@exprs
                                                 data_checked[is.na(data_checked)] <- 0
                                                 proteins_excluded <- union(proteins_excluded, colnames(data_checked[, c(colMeans(apply(data_checked, 2, as.numeric) == 0, na.rm = T) > 0.99)]))
                                               }

                                               AP2_norm_frame <- try(norm_anpela(data = AP2_trans_frame, method = workflow[i,3], index = index_TIclass[!(index_TIclass %in% proteins_excluded)]), silent = T)
                                               if (class(AP2_norm_frame) == "try-error") {
                                                 # AP2_norm_frame <- AP2_trans_frame
                                                 return(NULL)
                                               }
                                               rm(AP2_trans_frame)

                                               AP2_sigcl_frame <- try(sigcl_anpela(data = AP2_norm_frame, method = workflow[i,4], index = index_TIclass,
                                                                                   Segment = Segment,
                                                                                   Segment2 = Segment2,
                                                                                   min_cells = min_cells, max_bins = max_bins, step = step, technique = technique),
                                                                      silent = T)
                                               if (class(AP2_sigcl_frame) == "try-error") {
                                                 # AP2_sigcl_frame <- AP2_norm_frame
                                                 return(NULL)
                                               }
                                               rm(AP2_norm_frame)

                                               res <- lapply(AP2_sigcl_frame, function(x) {
                                                 x@exprs[is.infinite(x@exprs)] <- 0
                                                 x@exprs[is.nan(x@exprs)] <- 0
                                                 x@exprs[is.na(x@exprs)] <- 0
                                                 return(x)
                                               })
                                               names(res) <- names(AP2_sigcl_frame)
                                               rm(AP2_sigcl_frame)
                                               gc()

                                               if (save_processed_res == "one_folder"){
                                                 if (!dir.exists(paste0(savepath, "/process_res"))) {
                                                   dir.create(paste0(savepath, "/process_res"), recursive = T)
                                                 }
                                                 if (!file.exists(paste0(savepath, "/info_saved.RData"))) {
                                                   info_saved <- list(dataFileNames = dataFileNames, metadata = metadata, index_TIclass = index_TIclass)
                                                   save(info_saved, file = paste0(savepath, "/info_saved.RData"))
                                                 }
                                                 save(res, file = paste0(savepath, "/process_res/", rownames(workflow)[i], ".RData"))
                                               }else if(save_processed_res == "one_RData"){
                                                 return(res)
                                               }
                                             }

  parallel::stopCluster(cl)
  print(proc.time()-time)
  # parallel end

  if(save_processed_res == "one_RData"){
    names(AP2_pro1_frame_classTI) <- rownames(workflow)
    FCprocess_res <- list(AP2_pro1_frame_classTI = AP2_pro1_frame_classTI, dataFileNames = dataFileNames, metadata = metadata, index_TIclass = index_TIclass)
    if (!dir.exists(paste0(savepath, "/process_res"))) {
      dir.create(paste0(savepath, "/process_res"), recursive = T)
    }
    save(FCprocess_res, file = paste0(savepath, "/process_res/", name, ".RData"))

    return(FCprocess_res)
  }
}
