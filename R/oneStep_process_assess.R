oneStep_process_assess <- function(
    # the parameters of data processing
  datapath,
  metadata,
  technique = c("MC", "FC"),
  studytype = c("CSI", "PTI"),
  mergeM = c("Fixed", "Ceil", "All", "Min"),
  fixedNum = 200,
  compensationM = c("AutoSpill", "CATALYST", "CytoSpill", "FlowCore", "MetaCyto", "spillR", "None"),
  transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                      "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                      "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation",
                      "Centered Log Ratio Transformation", "None"),
  normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None"),
  signalcleanM = c("FlowAI", "FlowClean", "FlowCut", "PeacoQC", "None"),
  workflow = NULL,
  spillpath = NULL, FSC = "FSC-H", SSC = "SSC-H",
  control.dir = NULL, control.def.file = NULL,
  single_pos_fcs = NULL, single_pos_mass = NULL, CATALYSTM = "nnls",
  sce_bead = NULL, marker_to_barc = NULL,
  arcsinha = 0, arcsinhb = NULL, arcsinhc = 0,
  anna = 0, annb = NULL, annc = 0, annthreshold = 1,
  arna = 0, arnb = NULL, arnc = 0, arnthreshold = 1,
  bepa = 0.5, bepb = 1, bepc = 0.5, bepd = 1, bepf = 0, bepw = 0, tol = .Machine$double.eps^0.25, maxit = as.integer(5000),
  hpla = 1, hplb = 1,
  lineara = 2, linearb = 0,
  lntr = 1, lntd = 1,
  logbase = 10,logr = 1,logd = 1,
  lgtw = 0.5, lgtt = 262144, lgtm = 4.5, lgta = 0,
  Quadratica = 1, Quadraticb = 1, Quadraticc = 0,
  Truncatea = 1,
  beads_mass = NULL,
  Segment = 200,
  Segment2 = 200,
  min_cells = 3, max_bins = 10, step = 10,
  excludedColumn = NULL,

  # the parameters of workflow assessment
  name = "result",
  clusteringM = c("FlowSOM", "PhenoGraph"),
  ncluster = 8,
  Phenograph_k = 30,
  Ca_metric = "AUC",
  Cb_metric = "Silhouette coefficient (SC)",
  Cc_metric = "relative weighted consistency (CWrel)",
  ntop = NULL,
  DEP = NULL,
  marker_path = NULL, known_celltype_path = NULL,

  TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh",
          "slingshot_tSNE","slingshot_FLOWMAP", "slingshot_PCA", "slingshot_diffMaps"),
  pathwayhierarchy = NULL,
  clustering.var = NULL,

  # other parameters
  cores = floor(parallel::detectCores()/2),
  save_processed_res ="one_folder",
  savepath = "./ANPELA_res"
){

  ################################# the parameters of data processing #################################

  # dataFiles
  if (missing(datapath)) { # datapath parameter is not supplied
    stop("The parameter of 'datapath' is missing.")
  } else if (file.info(datapath)$isdir) { # provides the absolute path to the original data folder
    dataFiles <- list.files(datapath, pattern = ".fcs$", full.names = TRUE)
  } else if (!file.info(datapath)$isdir) { # provides the absolute path containing the FCS file and extracts the corresponding FCS file path
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
  if (missing(metadata)) { # No metadata parameter is provided
    stop("The parameter of 'metadata' is missing.")
  } else if (grepl(".csv$", metadata)) { # provides the absolute path to the csv file
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
    if (technique == "MC") {
      compensationM <- c("CATALYST", "CytoSpill", "spillR", "None")
    } else if (technique == "FC") {
      compensationM <- c("AutoSpill", "FlowCore", "MetaCyto", "None")
    }
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
    if (technique == "MC") {
      normalizationM <- c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None")
    } else if (technique == "FC") {
      normalizationM <- c("GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None")
    }
  } else {
    normalizationM <- match.arg(normalizationM, several.ok = TRUE)
  }


  # signalcleanM
  if (missing(signalcleanM)) {
    if (technique == "MC") {
      signalcleanM <- c("FlowAI", "FlowCut", "PeacoQC", "None")
    } else if (technique == "FC") {
      signalcleanM <- c("FlowAI", "FlowClean", "FlowCut", "PeacoQC", "None")
    }
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
  if ("CATALYST" %in% compensationM & is.null(workflow)|any(grepl("CATALYST", workflow))) {
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

  # arcsinhb
  if (is.null(arcsinhb)) {
    if (technique == "MC") {
      arcsinhb <- 1/5
    } else if (technique == "FC") {
      arcsinhb <- 1/150
    }
  } else if (arcsinhb == 0) {
    stop("The value of arcsinhb cannot be 0.")
  }


  # annb
  if (is.null(annb)) {
    if (technique == "MC") {
      annb <- 1/5
    } else if (technique == "FC") {
      annb <- 1/150
    }
  } else if (annb == 0) {
    stop("The value of annb cannot be 0.")
  }


  # arnb
  if (is.null(arnb)) {
    if (technique == "MC") {
      arnb <- 1/5
    } else if (technique == "FC") {
      arnb <- 1/150
    }
  } else if (arnb == 0) {
    stop("The value of arnb cannot be 0.")
  }


  # beads_mass
  if ("Bead-based Normalization" %in% normalizationM & is.null(workflow)|any(grepl("Bead-based Normalization", workflow))) {
    if (missing(beads_mass) | is.null(beads_mass)) {
      message("The parameter of 'beads_mass' is missing. 'Bead-based Normalization' normalization method can't be performed.")
      compensationM <- setdiff(normalizationM, "Bead-based Normalization")
    } else if (!is.numeric(beads_mass)) {
      stop("The parameter of 'beads_mass' is incorrect. Please input the masses of the corresponding calibration beads.")
    }
  }



  ################################# the parameters of workflow assessment #################################

  # clusteringM
  if (missing(clusteringM)) {
    clusteringM <- "FlowSOM"
  } else {
    clusteringM <- match.arg(clusteringM)
  }


  # ncluster
  if (clusteringM == "FlowSOM") {
    if (ncluster < 2 || ncluster %% 1 != 0) {
      stop("The number of clusters for meta clustering in FlowSOM should be positive integer number and greater than 1.")
    }
  }


  # TIM
  if (missing(TIM)) {
    TIM <- "scorpius_distSpear"
  } else {
    TIM <- match.arg(TIM)
  }


  # Ca_metric
  if (missing(Ca_metric)) {
    Ca_metric <- "AUC"
  } else {
    Ca_metric <- match.arg(Ca_metric)
  }


  # Cb_metric
  if (missing(Cb_metric)) {
    Cb_metric <- "Silhouette coefficient (SC)"
  } else {
    Cb_metric <- match.arg(Cb_metric)
  }


  # Cc_metric
  if (missing(Cc_metric)) {
    if (studytype == "CSI") {
      Cc_metric <- "relative weighted consistency (CWrel)"
    } else if (studytype == "PTI") {
      Cc_metric <- "Spearman rank correlation"
    }
  } else {
    Cc_metric <- match.arg(Cc_metric)
    if (studytype == "CSI" && Cc_metric %in% c("Spearman rank correlation", "Kendall rank correlation")){
      stop("The parameter of 'Cc_metric' is incorrect. Please select one of the following metrics:
           \n'relative weighted consistency (CWrel)', 'consistency score (CS)'")
    } else if (studytype == "PTI" && Cc_metric %in% c("relative weighted consistency (CWrel)", "consistency score (CS)")){
      stop("The parameter of 'Cc_metric' is incorrect. Please select one of the following metrics:
           \n'Spearman rank correlation', 'Kendall rank correlation'")
    }
  }



  ################################# data processing & workflow assessment #################################

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
  flag2 <- switch(studytype,
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

  # CSI/PTI down sample frame list
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
    if(technique == "MC"){
      cat("Please enter the excluded column names which are separated by the comma on a single line.
          \nFor example, Time(Time), Cell_length(Cell_length), DNA-1(DNA.1.Ir191.Dd)")
    } else if (technique == "FC"){
      cat("Please enter the excluded column names which are separated by the comma on a single line.
          \nFor example, Time(Time), FSC-W(FSC.W), SSC-A(SSC.A)")
    }
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

  # known_marker
  if (studytype == "CSI"){
    if (is.null(DEP)) {
      cat("*************************************************************************", "\n")
      cat("The standardized marker names are listed below. \n")
      cat("*************************************************************************", "\n")
      cat(paste0(colnames(AP2_pro0_frame[[1]]@exprs), collapse = "\n"), "\n")
      cat("*************************************************************************", "\n")
      cat("Please enter the marker names which are separated by comma on a single line.
        \nFor example, CD103(La139Di), CD11b(Nd144Di), CD8a(Nd146Di), CD7(Sm147Di)")
      DEP <- readline("Now, you can select the known biomarker(s) differentially expressed between two conditions (if NULL, please press the Enter key directly):")
    } else if (file.exists(DEP)) {
      DEP_data <- read.csv(DEP, header = F, stringsAsFactors = FALSE)
      DEP <-paste(as.character(unlist(DEP_data)), collapse = ",")
    } else if (!DEP == "" & !file.exists(DEP)){
      stop("The filepath of parameter 'DEP' is incorrect. Please input the absolute filepath of the csv file.")
    }
  }

  # ntop
  if (Cc_metric == "relative weighted consistency (CWrel)") {
    if (missing(ntop)) {
      ntop <- floor(length(index_TIclass)/2)
    } else if (ntop >= length(index_TIclass) || ntop %% 1 != 0 || ntop <= 0) {
      stop("The value of 'ntop' is incorrect. It should be positive whole number and less than the number of your selected markers.")
    }
  }

  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }

  # parallel start
  opts <- list(progress = function(n) setTxtProgressBar(txtProgressBar(min = 1, max = nrow(workflow), style = 3), n))
  cl <- parallel::makeCluster(cores, type = "SOCK", outfile = paste0(savepath,"/onestep_log.txt"))
  doSNOW::registerDoSNOW(cl)
  time = proc.time()

  table <- foreach::foreach(i = 1:nrow(workflow), .options.snow = opts,
                            .packages = c("dplyr", "flowCore", "foreach", "magrittr","mclust", "Rphenograph"), .combine = rbind) %dopar% {

                              ########### start data processing ###########

                              try(source("./processing.R"), silent = T)

                              AP2_comp_frame <- try(comp_anpela(data = AP2_pro0_frame, method = workflow[i,1], index = index_TIclass,
                                                                spillpath = spillpath, spillname = spillname, FSC = FSC,  SSC = SSC,
                                                                control.dir = control.dir, control.def.file = control.def.file,
                                                                single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass,
                                                                CATALYSTM = CATALYSTM,
                                                                sce_bead = sce_bead, marker_to_barc = marker_to_barc), silent = T)
                              if (class(AP2_comp_frame) == "try-error") {
                                res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                rownames(res) <- rownames(workflow)[i]
                                return(res)
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
                                                                  Truncatea = Truncatea), silent = T)
                              if (class(AP2_trans_frame) == "try-error") {
                                res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                rownames(res) <- rownames(workflow)[i]
                                return(res)
                              }
                              rm(AP2_comp_frame)

                              # Prevent too many 0 values in the data after compensation and affect normalization
                              proteins_excluded <- c()
                              for (x in seq(AP2_trans_frame)) {
                                data_checked <- AP2_trans_frame[[x]]@exprs
                                data_checked[is.na(data_checked)] <- 0
                                # Proteins that were not detected in more than 99% of the samples were recorded
                                proteins_excluded <- union(proteins_excluded, colnames(data_checked[, c(colMeans(apply(data_checked, 2, as.numeric) == 0, na.rm = T) > 0.99)]))
                              }

                              AP2_norm_frame <- try(norm_anpela(data = AP2_trans_frame, method = workflow[i,3],
                                                                index = index_TIclass[!(index_TIclass %in% proteins_excluded)],
                                                                beads_mass = beads_mass), silent = T)
                              if (class(AP2_norm_frame) == "try-error") {
                                res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                rownames(res) <- rownames(workflow)[i]
                                return(res)
                              }
                              rm(AP2_trans_frame)

                              AP2_sigcl_frame <- try(sigcl_anpela(data = AP2_norm_frame, method = workflow[i,4], index = index_TIclass,
                                                                  Segment = Segment, Segment2 = Segment2,
                                                                  min_cells = min_cells, max_bins = max_bins, step = step, technique = technique), silent = T)
                              if (class(AP2_sigcl_frame) == "try-error") {
                                res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                rownames(res) <- rownames(workflow)[i]
                                return(res)
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
                              }

                              ########### start workflow assessment ###########

                              if (studytype == "CSI") {
                                try(source("./CSI/1readfcs.R"), silent = T)
                                try(source("./CSI/2cluster.R"), silent = T)
                                try(source("./CSI/3criteria.R"), silent = T)
                                try(source("./CSI/4plot.R"), silent = T)

                                # AP2_processed_D_class
                                if (is.null(res)) {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- rownames(workflow)[i]
                                  return(res)
                                }
                                names(res) <- limma::removeExt(basename(dataFileNames), sep=".")
                                AP2_processed_D_class <- res
                                rm(res)


                                # AP2_processed_data_class0
                                data1 <- try(as.data.frame(readfcs_multi(fcsFiles = AP2_processed_D_class, mergeMethod = "all")), silent = T)
                                if (class(data1) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- rownames(workflow)[i]
                                  return(res)
                                }
                                # make.names(colnames(data1))
                                colnames(data1) <- sub("<", "_", colnames(data1))
                                colnames(data1) <- sub(">", "_", colnames(data1))
                                rm(AP2_processed_D_class)

                                condition <- metadata
                                data_with_filename <- cbind(data1, "filename" = sub("_[0-9]*$","", row.names(data1)))
                                condition_label <- as.character(condition[match(data_with_filename[,"filename"], condition$filename),2])
                                data_with_condition <- cbind(data_with_filename, "condition" = condition_label)
                                AP2_processed_data_class0 <- data_with_condition
                                rm(data1, condition, data_with_filename, condition_label, data_with_condition)


                                # AP2_processed_data_class
                                AP2_processed_data_class <- AP2_processed_data_class0[, c(index_TIclass, "filename", "condition")]
                                rm(AP2_processed_data_class0)


                                # cluster_label
                                cluster_label <- try(data_cluster(data = AP2_processed_data_class[,1:(dim(AP2_processed_data_class)[2]-2)],
                                                                  method = clusteringM, Phenograph_k = Phenograph_k,
                                                                  FlowSOM_k = ncluster, FlowSeed = 40), silent = T)
                                if (class(cluster_label) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- rownames(workflow)[i]
                                  return(res)
                                }


                                # data_with_cluster
                                data_with_cluster <- cbind(AP2_processed_data_class, "cluster" = cluster_label)


                                # test_KNN
                                subdata_cluster <- list()
                                subdata_cluster_DEG <- list()
                                train <- list()
                                test <- list()
                                KNN_res <- list()

                                for (j in 1:length(unique(cluster_label))) {

                                  subdata_cluster[[j]] <- subset(as.data.frame(data_with_cluster, stringsAsFactors = F), cluster == j)
                                  subdata_cluster[[j]]$condition <- as.factor(subdata_cluster[[j]]$condition)

                                  non_DEG <- try(feature_selection(subdata_cluster[[j]]), silent = T)
                                  if ((class(non_DEG) == "try-error") || (length(non_DEG) == length(subdata_cluster[[j]]) - 3)) next # The second case is that the number of non_DEG is equal to the number of markers in the original data, that is, subdata_cluster_DEG will have no data

                                  subdata_cluster_DEG[[j]] <- select(subdata_cluster[[j]], -one_of(non_DEG))

                                  data_condition1 <- subset(subdata_cluster_DEG[[j]], condition == unique(subdata_cluster_DEG[[j]]$condition)[1])
                                  data_condition2 <- subset(subdata_cluster_DEG[[j]], condition == unique(subdata_cluster_DEG[[j]]$condition)[2])

                                  set.seed(123)
                                  tr_ind1 <- sample(1:nrow(data_condition1), nrow(data_condition1)*0.7)
                                  tr_ind2 <- sample(1:nrow(data_condition2), nrow(data_condition2)*0.7)

                                  train[[j]] <- rbind(data_condition1[tr_ind1,], data_condition2[tr_ind2,])
                                  test[[j]] <- rbind(data_condition1[-tr_ind1,], data_condition2[-tr_ind2,])

                                  KNN_res[[j]] <- try(FNN::knn(train = train[[j]][, 1:(ncol(train[[j]]) - 3), drop = FALSE],
                                                               test = test[[j]][, 1:(ncol(test[[j]]) - 3), drop = FALSE],
                                                               cl = train[[j]]$condition, k = 5, prob = TRUE), silent = TRUE)
                                }

                                test_KNN <- list(test = test, KNN_res = KNN_res, subdata_cluster_DEG = subdata_cluster_DEG, subdata_cluster = subdata_cluster)
                                rm(subdata_cluster, subdata_cluster_DEG, train, test, KNN_res)


                                # sub_cluster_label
                                sub_cluster_label <- lapply(test_KNN$subdata_cluster_DEG, sub_cluster, FlowSeed = 40)

                                gc()

                                if (Ca_metric == "AUC") {
                                  # Criterion A Accuracy-AUC
                                  resAUC <- try(AUC(test = test_KNN$test, KNN_res = test_KNN$KNN_res), silent = T)
                                  Cauc <- try(round(sum(resAUC$auc, na.rm = TRUE)/length(test_KNN[["subdata_cluster"]]), 5), silent = T)
                                  if (class(Cauc) != "numeric") {
                                    Cauc <- NA
                                  }
                                  Ca <- Cauc
                                  rm(resAUC, Cauc)
                                } else if (Ca_metric == "F1 score") {
                                  # Criterion A Accuracy-F1 score
                                  resF1 <- try(F1_score(test = test_KNN$test, KNN_res = test_KNN$KNN_res, label = AP2_processed_data_class$condition), silent = T)
                                  CF1_score <- try(round(max(sapply(resF1, mean, na.rm = TRUE)), 5), silent = T)
                                  if (class(CF1_score) != "numeric") {
                                    CF1_score <- NA
                                  }
                                  Ca <- CF1_score

                                  rm(resF1, CF1_score)
                                }


                                if (Cb_metric == "Silhouette coefficient (SC)") {
                                  # Criterion B Silhouette coefficient (SC)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Silhouette")[[1]], silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  } else {
                                    Cb <- round((Cb + 1)/2, 5)
                                  }
                                } else if (Cb_metric == "Xie-Beni index (XB)") {
                                  # Criterion B Xie-Beni index (XB)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Xie_Beni")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "Calinski-Harabasz index (CH)") {
                                  # Criterion B Calinski-Harabasz index (CH)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Calinski_Harabasz")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "Davies-Bouldin index (DB)") {
                                  # Criterion B Davies-Bouldin index (DB)
                                  datab <- AP2_processed_data_class[, 1:(dim(AP2_processed_data_class)[2]-2)]
                                  Cb <- try(round((clusterCrit::intCriteria(as.matrix(datab), as.integer(cluster_label), "Davies_Bouldin")[[1]]), 5), silent = T)
                                  if (class(Cb) != "numeric") {
                                    Cb <- NA
                                  }
                                } else if (Cb_metric == "purity") {
                                  # Criterion B Precision-purity
                                  respurity <- try(Purity(data = test_KNN$subdata_cluster_DEG, sub_cluster_label = sub_cluster_label, FlowSeed = 40), silent = T)
                                  Cpurity <- try(round(mean(respurity, na.rm = TRUE), 5), silent = T)
                                  if (class(Cpurity) != "numeric") {
                                    Cpurity <- NA
                                  }
                                  Cb <- Cpurity
                                  rm(respurity, Cpurity)
                                } else if (Cb_metric == "Rand index (RI)") {
                                  # Criterion B Precision-RI
                                  resCRI <- try(RI(data = test_KNN$subdata_cluster_DEG, sub_cluster_label = sub_cluster_label, FlowSeed = 40), silent = T)
                                  CRI <- try(round(mean(resCRI, na.rm = TRUE), 5), silent = T)
                                  if (class(CRI) != "numeric") {
                                    CRI <- NA
                                  }
                                  Cb <- CRI
                                  rm(resCRI, CRI)
                                }
                                rm(sub_cluster_label, AP2_processed_data_class, cluster_label, datab)


                                # Criterion C Robustness-CS_pre
                                CS_preres <- try(CS_pre(test_KNN$subdata_cluster), silent = T)
                                # Robustness-CS, CW
                                if (Cc_metric == "consistency score (CS)") {
                                  resCS <- try(CSfun(CS_preres = CS_preres), silent = T)
                                } else if (Cc_metric == "relative weighted consistency (CWrel)") {
                                  resCS <- try(CWfun(CS_preres = CS_preres, top = ntop), silent = T)
                                }
                                CS <- try(round(mean(resCS$consistency, na.rm = TRUE), 5), silent = T)
                                if (class(CS) != "numeric") {
                                  CS <- NA
                                }
                                Cc <- CS
                                rm(test_KNN, CS_preres, resCS, CS)


                                # Criterion D Biological Meaning
                                if ((is.null(marker_path)||is.null(known_celltype_path))& (DEP == "" || is.null(DEP))) {
                                  Cd <- NA
                                } else if (!is.null(marker_path) && !is.null(known_celltype_path)) {
                                  CRecall <- try(round(AP2_Recall(data_with_cluster = data_with_cluster,
                                                                  marker_path = marker_path, known_celltype_path = known_celltype_path), 5), silent = T)
                                  if (class(CRecall) != "numeric") {
                                    CRecall <- NA
                                  }
                                  Cd <- CRecall
                                } else if (!is.null(DEP)) {
                                  known_marker <- unlist(strsplit(DEP, "\\s*,\\s*"))
                                  CRecall <- try(round(AP2_Recall(data_with_cluster = data_with_cluster, known_marker = known_marker), 5), silent = T)
                                  if (class(CRecall) != "numeric") {
                                    CRecall <- NA
                                  }
                                  Cd <- CRecall
                                  rm(known_marker, CRecall)
                                }
                                rm(data_with_cluster)
                                gc()

                                res <- data.frame(Ca, Cb, Cc, Cd)
                                rownames(res) <- rownames(workflow)[i]
                                return(res)

                              } else if (studytype == "PTI") {

                                try(source("./PTI/load_data2.R"), silent = T)
                                try(source("./PTI/TI_method.R"), silent = T)

                                try(source("./PTI/Bio_con_4.R"), silent = T)
                                try(source("./PTI/time_metric_3.R"), silent = T)
                                try(source("./PTI/Robustness_4.R"), silent = T)
                                try(source("./PTI/Rough_3.R"), silent = T)

                                try(source("./PTI/shift_start.R"), silent = T)
                                try(source("./PTI/calc_spline.R"), silent = T)
                                try(source("./PTI/cycle_pseudotime.R"), silent = T)
                                try(source("./PTI/reverse_pseudotime.R"), silent = T)
                                try(source("./PTI/check_pairs.R"), silent = T)
                                try(source("./PTI/ANPELA_FLOWMAP.R"))
                                try(source("./PTI/ANPELA_FLOWMAP-function.R"))
                                try(source("./PTI/plot.R"), silent = T)

                                # AP2_processed_D_TI
                                index <- stringr::str_replace_all(index_TIclass, "\\(.*", "")
                                if ("condition" %in% colnames(metadata)){
                                  res <- try(load.Data(res, index = index, measurement.condition = as.matrix(metadata$condition),
                                                       measurement.time = as.matrix(metadata$timepoint), TIM = TIM), silent = T)
                                } else {
                                  res <- try(load.Data(res, index = index, measurement.time = as.matrix(metadata$timepoint), TIM = TIM), silent = T)
                                }
                                dataset_name <-rownames(workflow)[i]
                                if (class(res) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- rownames(workflow)[i]
                                  return(res)
                                }

                                AP2_processed_D_TI <- res
                                rm(index, res)

                                # TIress
                                TIres <- try(TI(D = AP2_processed_D_TI, method = TIM,
                                                dataset_name = dataset_name,clustering.var = clustering.var))
                                print(paste0("TIres","_",i))

                                if(class(TIres) == "try-error") {
                                  res <- data.frame(Ca = NA, Cb = NA, Cc = NA, Cd = NA)
                                  rownames(res) <- rownames(workflow)[i]
                                  print(paste0("TIres","error"))
                                  return(res)
                                }



                                # Criterion A Consistency-Time
                                Ca <- try(time_metric(TIres, AP2_processed_D_TI, nruns = 100), silent = T)
                                if (class(Ca) != "numeric") {
                                  Ca <- NA
                                } else {
                                  Ca <- round(Ca, 5)
                                }


                                # Criterion B Roughness-Roughness
                                R <- try(Rough(TIres, AP2_processed_D_TI), silent = T)
                                Cb <- try(R$p.value, silent = T)
                                if (class(Cb) != "numeric") {
                                  Cb <- NA
                                } else {
                                  if (Cb > 0.05) {
                                    Cb <- 0
                                  } else {
                                    Cb <- round(1 - 20*Cb, 5)
                                  }
                                }
                                rm(R)

                                # Criterion C Robustness-Robustness
                                Rob0 <- try(Robustness(TIres, AP2_processed_D_TI, nruns = 4, cell.subset = 0.8, clustering.var = clustering.var,
                                                       dataset_name = dataset_name), silent = T)
                                if (class(Rob0) != "list") {
                                  Cc <- NA
                                } else {
                                  Cc <- switch (Cc_metric,
                                                "Spearman rank correlation" = round(Rob0$Robustness_result[1], 5),
                                                "Kendall rank correlation" = round(Rob0$Robustness_result[2], 5)
                                  )
                                }
                                rm(Rob0)


                                # Criterion D Biological Meaning-Biological consistency
                                if (!is.null(pathwayhierarchy) && file.exists(pathwayhierarchy)) {
                                  Cd <- try(suppressWarnings(Bio_con(AP2_processed_D_TI, Pathway_Hierarchy_file = pathwayhierarchy, nruns = 3, dr_method = TIres$dr_method, TIres = TIres)), silent = T)
                                  if (class(Cd) != "numeric") {
                                    Cd <- NA
                                  } else {
                                    Cd <- round(Cd, 5)
                                  }
                                } else if(is.null(pathwayhierarchy) && "condition" %in% colnames(metadata)){
                                  Cd <- try(suppressWarnings(Bio_con(AP2_processed_D_TI,  nruns = 3, dr_method = TIres$dr_method, TIres = TIres)), silent = T)
                                  if (class(Cd) != "numeric") {
                                    Cd <- NA
                                  } else {
                                    Cd <- round(Cd, 5)
                                  }
                                } else Cd <- NA

                                rm(AP2_processed_D_TI, TIres)
                                gc()

                                res <- data.frame(Ca, Cb, Cc, Cd)
                                rownames(res) <- rownames(workflow)[i]

                                return(res)
                              }
                            }


  parallel::stopCluster(cl)
  print(proc.time()-time)
  # parallel end

  if (studytype == "CSI") {

    colnames(table) <- c("Accuracy", "Tightness", "Robustness", "Correspondence")
    table2 <- table
    table2["Accuracy"][table2["Accuracy"] > 0.7] <- 10
    table2["Accuracy"][table2["Accuracy"] <= 0.7] <- 4

    table2["Tightness"][table2["Tightness"] > 0.5] <- 10
    table2["Tightness"][table2["Tightness"] <= 0.5] <- 4

    table2["Robustness"][table2["Robustness"] > 0.35] <- 10
    table2["Robustness"][table2["Robustness"] <= 0.35] <- 4

    if (!is.null(marker_path) && !is.null(known_celltype_path)) {
      table2["Correspondence"][table2["Correspondence"] > 0.7] <- 10
      table2["Correspondence"][table2["Correspondence"] <= 0.5] <- 4
    } else {
      table2["Correspondence"][table2["Correspondence"] > 0.5] <- 10
      table2["Correspondence"][table2["Correspondence"] <= 0.5] <- 4
    }

  } else if (studytype == "PTI") {

    colnames(table) <- c("Conformance", "Smoothness", "Robustness", "Correspondence")
    table2 <- table
    table2["Conformance"][table2["Conformance"] > 0.6] <- 10
    table2["Conformance"][table2["Conformance"] <= 0.6] <- 4

    table2["Smoothness"][table2["Smoothness"] > 0.8] <- 10
    table2["Smoothness"][table2["Smoothness"] <= 0.8] <- 4

    table2["Robustness"][table2["Robustness"] > 0.5] <- 10
    table2["Robustness"][table2["Robustness"] <= 0.5] <- 4

    table2["Correspondence"][table2["Correspondence"] > 0.6] <- 10
    table2["Correspondence"][table2["Correspondence"] <= 0.6] <- 4
  }

  assess_res <- list(table = table, table2 = table2)
  if (!dir.exists(paste0(savepath, "/assess_res"))) {
    dir.create(paste0(savepath, "/assess_res"), recursive = T)
  }
  save(assess_res, file = paste0(savepath, "/assess_res/", name, "_assess.RData"))
  return(assess_res)
}
