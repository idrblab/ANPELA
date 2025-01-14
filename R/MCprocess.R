
#' @title Mass Cytometry Process
#' @description MCprocess() enables high-throughput processing for SCP data acquired from MC by at most ~675 available workflows based on parallel computing (each workflow is distinct by combining different methods of compensation, transformation, normalization and signal clean), which facilitates the subsequent application of performance assessment, ranking and plotting.
#'
#' @param name Character, the filename of the RData file when the "save_processed_res" parameter is set to "one_RData".
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#' @param metadata Character, the absolute filepath of the metadata file, usually located in the datapath folder. The metadata file should include the columns of "filename" and "condition" for CSI studies, and the columns of "filename" and "timepoint" for PTI studies.
#'   <br>For details on preparing the metadata file, please refer to the **sample data**.
#' @param studytype Character, the type of your study, including "CSI (Cell Subpopulation Identification)" and "PTI (Pseudotime Trajectory Inference)".
#' @param mergeM Character, the method of merging multiple FCS files. When multiple FCS files are selected, cells can be combined using one of the four different methods including "Fixed", "Ceil", "All" and "Min".
#'   <br>**Fixed**: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each FCS file and combined for analysis.
#'   <br>**Ceil**: up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each FCS file and combined for analysis.
#'   <br>**All**: all cells from each FCS file are combined for analysis.
#'   <br>**Min**: The minimum number of cells among all the selected FCS files are sampled from each FCS file and combined for analysis.
#' @param fixedNum Integer, the fixed number of cells to be extracted from each FCS file.
#' @param compensationM Character, the method(s) of compensation for flow cytometry data including "CATALYST", "CytoSpill", "spillR" and "None". Compensation refers to the processing step of removing unwanted spillover resulting from signal crosstalk and spectral overlap across detection channels.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param transformationM Character, the method(s) of transformation for flow cytometry data including "Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value", "Biexponential Transformation", "Box-Cox Transformation", "Centered Log Ratio Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation", "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation" and "None". Transformation refers to the processing step of adjusting the data with a heavily skewed distribution to a normal distribution.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param normalizationM Character, the method(s) of normalization for flow cytometry data, including "Bead-based Normalization", "GaussNorm", "Mean Normalization", "Min-max Normalization", "WarpSet", "ZScore" and "None". Normalization refers to the processing step of eliminating signal decay and technical variability across all files and batches over long-term data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param signalcleanM Character, the method(s) of signal clean for flow cytometry data, including "FlowAI", "FlowCut", "PeacoQC" and "None". Signal cleaning refers to the processing step of identifying and removing abrupt signal shifts and changes that derive from (i) abrupt changes in the flow rate, (ii) clogs within the capillary tubes, (iii) temporary disruptions in cytometer fluidics, and (iv) unstable data acquisition.
#'   <br>For details of each method, please refer to the **Methods Introduction**.
#' @param workflow Character, the combinations of data processing methods specified by users according to their research interests.
#'   <br>It is a vector includes one or more method combinations, typically in the format of "compensation method name_ transformation method name_ normalization method name_ signal clean method name ", for example: c("None_Biexponential Transformation_None_None","CytoSpill_FlowVS Transformation_None_FlowCut").
#' @param single_pos_fcs Character, the absolute filepath of the FCS file containing stained samples and control antibody-capture beads/pooled single-stained beads.
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param single_pos_mass Integer, the masses corresponding to barcode channels.
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param CATALYSTM Character, the method for solving linear system, including "flow" and "nnls".
#'   <br>Only needed when "CATALYST" is included in the argument of "compensationM".
#' @param sce_bead SingleCellExperiment, the SingleCellExperiment object for the bead experiment.
#'   <br>Only needed when "spillR" is included in the argument of "compensationM".
#' @param marker_to_barc Data frame, the table that maps the marker to the barcode in the beads experiment.
#'   <br>Only needed when "spillR" is included in the argument of "compensationM".
#' @param logbase Integer, the base of the Log Transformation.
#'   <br>Only needed when "Log Transformation" is included in the argument of "transformationM".
#' @param b1 Double, the cofactor of Arcsinh Transformation.
#'   <br>Only needed when "Arcsinh Transformation" is included in the argument of "transformationM".
#' @param b2 Double, the cofactor of Asinh with Non-negative Value.
#'   <br>Only needed when "Asinh with Non-negative Value" is included in the argument of "transformationM".
#' @param b3 Double, the cofactor of Asinh with Randomized Negative Value.
#'   <br>Only needed when "Asinh with Randomized Negative Value" is included in the argument of "transformationM".
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
#' @param beads_mass Integer, the masses of the corresponding calibration beads.
#'   <br>Only needed when "Bead-based Normalization" is included in the argument of "normalizationM".
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
#' @param savepath Character, the absolute path of the folder which will store files of the processed results.
#' @return The **process_res** folder stores the results of various data processing workflows. The form of data processing output files is decided by the parameter `save_processed_res`: "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder; "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder \[*default ="one_folder"*\].
#'   <br>In addition, the file **info_saved.RData** is also generated simultaneously, recording the information related to "metadata" and "index_protein".
#' @export
#'
#' @examples
#' \donttest{
#' }

MCprocess <- function(name = "result",
                      datapath,
                      metadata,
                      studytype = c("CSI", "PTI"),
                      mergeM = c("Fixed", "Ceil", "All", "Min"),
                      fixedNum = 200,
                      compensationM = c("CATALYST", "CytoSpill", "spillR", "None"),
                      transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                                          "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                                          "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation",
                                          "Centered Log Ratio Transformation", "None"),
                      normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None"),
                      signalcleanM = c("FlowAI", "FlowCut", "PeacoQC", "None"),
                      workflow = NULL,
                      single_pos_fcs = NULL, single_pos_mass = NULL, CATALYSTM = c("flow", "nnls"),
                      sce_bead = NULL, marker_to_barc = NULL,
                      logbase = 10,
                      b1 = 1/5,
                      b2 = 1/5,
                      b3 = 1/5,
                      Quadratica = 1, Quadraticb = 1, Quadraticc = 0,
                      lineara = 2, linearb = 0,
                      Truncatea = 1,
                      beads_mass = c(140, 151, 153, 165, 175),
                      min_cells = 3, max_bins = 10, step = 10,
                      excludedColumn = NULL,
                      save_processed_res = "one_folder",
                      savepath = "./",
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
    compensationM <- c("CATALYST", "CytoSpill", "spillR", "None")
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
    normalizationM <- c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None")
  } else {
    normalizationM <- match.arg(normalizationM, several.ok = TRUE)
  }


  # signalcleanM
  if (missing(signalcleanM)) {
    signalcleanM <- c("FlowAI", "FlowCut", "PeacoQC", "None")
  } else {
    signalcleanM <- match.arg(signalcleanM, several.ok = TRUE)
  }


  # single_pos_fcs
  if ("CATALYST" %in% compensationM & is.null(workflow)|any(grepl("CATALYST", workflow))) {
    if (is.null(single_pos_fcs)) {
      message("The parameter of 'single_pos_fcs' is missing. 'CATALYST' compensation method can't be performed.")
      compensationM <- setdiff(compensationM, "CATALYST")
    } else if (is.character(single_pos_fcs) &&  !file.exists(single_pos_fcs)) {
      stop("The parameter of 'single_pos_fcs' is incorrect. Please input the absolute filepath of the .fcs file containing stained samples and control antibody-capture beads/pooled single-stained beads.")
    }
  }


  # single_pos_mass
  if ("CATALYST" %in% compensationM & is.null(workflow)|any(grepl("CATALYST", workflow))) {
    if (!is.numeric(single_pos_mass)) {
      stop("The parameter of 'single_pos_mass' is incorrect. Please input a vector of numeric masses corresponding to barcode channels.")
    }
  }


  # CATALYSTM
  if ("CATALYST" %in% compensationM) {
    if (missing(CATALYSTM)) {
      CATALYSTM <- "nnls"
    } else {
      CATALYSTM <- match.arg(CATALYSTM)
    }
  }

  #sce_bead
  if ("spillR" %in% compensationM & is.null(workflow)|any(grepl("spillR", workflow))) {
    if (is.null(sce_bead)) {
      message("The parameter of 'sce_bead' is missing. 'spillR' compensation method can't be performed.")
      compensationM <- setdiff(compensationM, "spillR")
    } else if (!"SingleCellExperiment" %in% class(sce_bead)) {
      stop("The parameter of 'sce_bead' is incorrect. Please input the SingleCellExperiment for the bead experiment.")
    }
  }

  #marker_to_barc
  if ("spillR" %in% compensationM & is.null(workflow)|any(grepl("spillR", workflow))) {
    if (is.null(marker_to_barc)) {
      message("The parameter of 'marker_to_barc' is missing. 'spillR' compensation method can't be performed.")
      compensationM <- setdiff(compensationM, "spillR")
    } else if (!is.data.frame(marker_to_barc)) {
      stop("The parameter of 'marker_to_barc' is incorrect. Please input a dataframe that maps the marker to the barcode in the beads experiment.")
    }
  }


  # logbase
  if (logbase <= 0 || logbase == 1) {
    stop("The value of logbase cannot be less than or equal to 0 or equal to 1.")
  }


  # b1
  if (b1 == 0) {
    stop("The value of b1 cannot be 0.")
  }


  # b2
  if (b2 == 0) {
    stop("The value of b2 cannot be 0.")
  }


  # b3
  if (b3 == 0) {
    stop("The value of b3 cannot be 0.")
  }


  # beads_mass
  if ("Bead-based Normalization" %in% normalizationM & is.null(workflow)|any(grepl("Bead-based Normalization", workflow))) {
    if (missing(beads_mass)) {
      message("The parameter of 'beads_mass' is missing. 'Bead-based Normalization' normalization method can't be performed.")
      compensationM <- setdiff(normalizationM, "Bead-based Normalization")
    } else if (!is.numeric(beads_mass)) {
      stop("The parameter of 'beads_mass' is incorrect. Please input the masses of the corresponding calibration beads.")
    }
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


  # index_TIclass
  if (is.null(excludedColumn)) {
    cat("*************************************************************************", "\n")
    cat("The standardized column names are listed below. \n")
    cat("*************************************************************************", "\n")
    cat(paste0(common_protein, collapse = "\n"), "\n")
    cat("*************************************************************************", "\n")
    cat("Please enter the excluded column names which are separated by the comma on a single line.
          \nFor example, Time(Time), Cell_length(Cell_length), DNA-1(DNA.1.Ir191.Dd)")

    excludedColumn <- readline("Now, you can select the excluded columns for data processing and performance evaluation:")
    if (excludedColumn == "") {
      stop("Note: Please select your excluded columns and then press the Enter key :(")
    }
  }
  Excluded_index_TIclass <- unlist(strsplit(excludedColumn, "\\s*,\\s*"))
  index_TIclass <- setdiff(common_protein,Excluded_index_TIclass)


  cat("The program is running. Please wait patiently.")

  # Segment
  Segment <- floor(min(ncell_frame)/3)

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
  cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/MCprocess_log.txt"))
  doSNOW::registerDoSNOW(cl)
  time = proc.time()

  AP2_pro1_frame_classTI <- foreach::foreach(i = 1:nrow(workflow), .options.snow = opts,
                                             .packages = c("stringr", "flowStats", "flowAI", "flowCore", "flowCut", "CytoSpill", "tree", "flowTrans", "magrittr")) %dopar% {
                                               try(source("./processing.R"), silent = T)
                                               set.seed(123)
                                               AP2_comp_frame <- try(comp_anpela(data = AP2_pro0_frame, method = workflow[i,1], index = index_TIclass,
                                                                                 single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass, CATALYSTM = CATALYSTM,
                                                                                 sce_bead = sce_bead, marker_to_barc = marker_to_barc), silent = T)
                                               if (class(AP2_comp_frame) == "try-error") {
                                                 # AP2_comp_frame <- AP2_pro0_frame
                                                 return(NULL)
                                               }

                                               AP2_trans_frame <- try(trans_anpela(data = AP2_comp_frame, method = workflow[i,2], index = index_TIclass,
                                                                                   logbase = logbase,
                                                                                   b1 = b1,
                                                                                   b2 = b2,
                                                                                   b3 = b3,
                                                                                   Quadratica = Quadratica, Quadraticb = Quadraticb, Quadraticc = Quadraticc,
                                                                                   lineara = lineara, linearb = linearb,
                                                                                   Truncatea = Truncatea), silent = T)
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

                                               AP2_norm_frame <- try(norm_anpela(data = AP2_trans_frame, method = workflow[i,3], index = index_TIclass[!(index_TIclass %in% proteins_excluded)], beads_mass = beads_mass), silent = T)
                                               if (class(AP2_norm_frame) == "try-error") {
                                                 # AP2_norm_frame <- AP2_trans_frame
                                                 return(NULL)
                                               }
                                               rm(AP2_trans_frame)

                                               AP2_sigcl_frame <- try(sigcl_anpela(data = AP2_norm_frame, method = workflow[i,4], index = index_TIclass,
                                                                                   Segment = Segment,
                                                                                   min_cells = min_cells, max_bins = max_bins, step = step, technique = "MC"),
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
                                                 if (!file.exists(paste0(savepath, "/metadata.RData"))) {
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
    MCprocess_res <- list(AP2_pro1_frame_classTI = AP2_pro1_frame_classTI, dataFileNames = dataFileNames, metadata = metadata, index_TIclass = index_TIclass)
    if (!dir.exists(paste0(savepath, "/process_res"))) {
      dir.create(paste0(savepath, "/process_res"), recursive = T)
    }
    save(MCprocess_res, file = paste0(savepath, "/process_res/", name, ".RData"))

    return(MCprocess_res)
  }
}
