# library
library(magrittr)
library(foreach)
library(dplyr)
library(flowCore)
# source
try(source("./FCprocess.R"), silent = T)
try(source("./MCprocess.R"), silent = T)
try(source("./CSIassess.R"), silent = T)
try(source("./PTIassess.R"), silent = T)
# try(source("./PTIassess_for.R"), silent = T)
try(source("./ranking.R"), silent = T)
try(source("./processing.R"), silent = T)



Getmarker <- function(datapath) {
  files <- dir(path = datapath, pattern = "fcs$", full.names = TRUE)
  a <- flowCore::read.FCS(filename = files[[1]], transformation = FALSE, alter.names = TRUE)
  colnames(a@exprs) <- paste0(a@parameters@data$desc, "(", a@parameters@data$name, ")")
  cat(paste0("\"", paste0(as.character(colnames(a@exprs)), collapse = ", \n"), "\""))
}


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
                   lineara = 2, linearb = 0, 
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

MC_CSI <- function(name = "result", datapath, technique = "MC",
                   mergeM = "Fixed", fixedNum = 200,
                   compensationM = c("CATALYST", "CytoSpill", "spillR", "None"), 
                   transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value", 
                                       "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation", 
                                       "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation", 
                                       "Centered Log Ratio Transformation", "None"),
                   normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None"),
                   signalcleanM = c("FlowAI", "FlowCut", "PeacoQC", "None"),
                   workflow = NULL,
                   single_pos_fcs = NULL, single_pos_mass = NULL, CATALYSTM ="nnls",
                   sce_bead = NULL, marker_to_barc = NULL,
                   arcsinha = 0, arcsinhb = 1/5, arcsinhc = 0,
                   anna = 0, annb = 1/5, annc = 0, annthreshold = 1, 
                   arna = 0, arnb = 1/5, arnc = 0, arnthreshold = 1, 
                   bepa = 0.5, bepb = 1, bepc = 0.5, bepd = 1, bepf = 0, bepw = 0, tol = .Machine$double.eps^0.25, maxit = as.integer(5000),
                   hpla = 1, hplb = 1, 
                   lineara = 2, linearb = 0, 
                   lntr = 1, lntd = 1, 
                   logbase = 10,logr = 1,logd = 1,
                   lgtw = 0.5, lgtt = 262144, lgtm = 4.5, lgta = 0,
                   Quadratica = 1, Quadraticb = 1, Quadraticc = 0,
                   lineara = 2, linearb = 0,
                   Truncatea = 1,
                   beads_mass = NULL,
                   Segment = 200,
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
    MCprocess_res <- MCprocess(name = name,
                               datapath = datapath,
                               metadata = metadata, technique = technique,
                               studytype = "CSI", mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               workflow = workflow,
                               single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass,
                               CATALYSTM = CATALYSTM,
                               sce_bead = sce_bead, marker_to_barc = marker_to_barc,
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
                               beads_mass = beads_mass,
                               Segment = Segment, 
                               min_cells = min_cells, max_bins = max_bins, step = step,
                               excludedColumn = excludedColumn,
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)
    MCprocess_assess_res <- CSIassess(name = name, data = MCprocess_res, clusteringM = clusteringM, 
                                      Phenograph_k = Phenograph_k, ncluster = ncluster, DEP = DEP, 
                                      save_processed_res = save_processed_res,savepath = savepath, 
                                      cores = cores)
  } else if (save_processed_res == "one_folder") {
    MCprocess_assess_res <<- oneStep_process_assess(
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
      single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass,
      CATALYSTM = CATALYSTM,
      sce_bead = sce_bead, marker_to_barc = marker_to_barc,
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
      beads_mass = beads_mass,
      Segment = Segment, 
      min_cells = min_cells, max_bins = max_bins, step = step,
      excludedColumn = excludedColumn, 
      
      clusteringM = clusteringM, Phenograph_k = Phenograph_k, 
      ncluster = ncluster, DEP = DEP, 
      
      cores = cores,
      save_processed_res = save_processed_res,
      savepath = savepath
    )
  }
  
  Ranking(MCprocess_assess_res, name = name, savepath = paste0(savepath, "/assess_res"))
}

FC_PTI <- function(name = "result", datapath, technique = "FC",
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
                   lineara = 2, linearb = 0, 
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

                   TIM = NULL,
                   pathwayhierarchy = NULL,
                   clustering.var = NULL,
                   
                   cores = parallel::detectCores()/2,
                   save_processed_res = "one_folder",
                   savepath = "./ANPELA_res"
) {
  
  metadata <- paste0(datapath, "/metadata.csv")
  if (save_processed_res == "one_RData") {
    FCprocess_res <- FCprocess(name = name,
                               datapath = datapath, 
                               metadata = metadata, technique = technique,
                               studytype = "PTI", mergeM = mergeM, fixedNum = fixedNum,
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
    FCprocess_assess_res <- PTIassess(name = name,data = FCprocess_res, TIM = TIM, 
                                      pathwayhierarchy = pathwayhierarchy, 
                                      clustering.var = NULL,
                                      
                                      save_processed_res = save_processed_res,
                                      savepath = savepath, cores = cores)
  }  else if (save_processed_res == "one_folder") {
    FCprocess_assess_res <- oneStep_process_assess(
      name = name,
      datapath = datapath,
      metadata = metadata,
      technique = technique,
      studytype = "PTI", mergeM = mergeM, fixedNum = fixedNum,
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
      
      TIM = TIM, pathwayhierarchy = pathwayhierarchy,
      clustering.var = clustering.var, 
      cores = cores,
      save_processed_res = save_processed_res,
      savepath = savepath
    )
  }

  Ranking(FCprocess_assess_res, name = name, savepath = paste0(savepath, "/assess_res"))
  
}

MC_PTI <- function(name = "result", datapath, technique = "MC",
                   mergeM = "Fixed", fixedNum = 200,
                   compensationM = c("CATALYST", "CytoSpill", "spillR", "None"), 
                   transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value", 
                                       "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation", 
                                       "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Split Scale Transformation", "Truncate Transformation", 
                                       "Centered Log Ratio Transformation", "None"),
                   normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "Mean Normalization", "Min-max Normalization", "None"),
                   signalcleanM = c("FlowAI", "FlowCut", "PeacoQC", "None"),
                   workflow = NULL,
                   single_pos_fcs = NULL, single_pos_mass = NULL, CATALYSTM ="nnls",
                   sce_bead = NULL, marker_to_barc = NULL,
                   arcsinha = 0, arcsinhb = 1/5, arcsinhc = 0,
                   anna = 0, annb = 1/5, annc = 0, annthreshold = 1, 
                   arna = 0, arnb = 1/5, arnc = 0, arnthreshold = 1, 
                   bepa = 0.5, bepb = 1, bepc = 0.5, bepd = 1, bepf = 0, bepw = 0, tol = .Machine$double.eps^0.25, maxit = as.integer(5000),
                   hpla = 1, hplb = 1, 
                   lineara = 2, linearb = 0, 
                   lntr = 1, lntd = 1, 
                   logbase = 10,logr = 1,logd = 1,
                   lgtw = 0.5, lgtt = 262144, lgtm = 4.5, lgta = 0,
                   Quadratica = 1, Quadraticb = 1, Quadraticc = 0,
                   lineara = 2, linearb = 0,
                   Truncatea = 1,
                   beads_mass = NULL,
                   Segment = 200,
                   min_cells = 3, max_bins = 10, step = 10,
                   excludedColumn = NULL, 
                   
                   TIM = NULL,
                   pathwayhierarchy = NULL, 
                   clustering.var = NULL,
                   
                   cores = parallel::detectCores()/2,
                   save_processed_res = "one_folder",
                   savepath = "./ANPELA_res"
) {
  metadata <- paste0(datapath, "/metadata.csv")
  if (save_processed_res == "one_RData") {
    MCprocess_res <- MCprocess(name = name,
                               datapath = datapath, 
                               metadata = metadata,
                               technique = technique,
                               studytype = "PTI", mergeM = mergeM, fixedNum = fixedNum, 
                               compensationM = compensationM, 
                               transformationM = transformationM,
                               normalizationM = normalizationM, 
                               signalcleanM = signalcleanM,
                               workflow = workflow,
                               single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass,
                               CATALYSTM = CATALYSTM,
                               sce_bead = sce_bead, marker_to_barc = marker_to_barc,
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
                               beads_mass = beads_mass,
                               Segment = Segment, 
                               min_cells = min_cells, max_bins = max_bins, step = step,
                               excludedColumn = excludedColumn, 
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)
    MCprocess_assess_res <- PTIassess(name = name, data = MCprocess_res, TIM = TIM, 
                                      pathwayhierarchy = pathwayhierarchy, 
                                      clustering.var = NULL,
                                      
                                      save_processed_res = save_processed_res,
                                      savepath = savepath, cores = cores)
  }  else if (save_processed_res == "one_folder") {
    MCprocess_assess_res <- oneStep_process_assess(
      name = name,
      datapath = datapath,
      metadata = metadata,
      technique = technique,
      studytype = "PTI", mergeM = mergeM, fixedNum = fixedNum,
      compensationM = compensationM, 
      transformationM = transformationM,
      normalizationM = normalizationM, 
      signalcleanM = signalcleanM,
      workflow = workflow,
      single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass,
      CATALYSTM = CATALYSTM,
      sce_bead = sce_bead, marker_to_barc = marker_to_barc,
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
      beads_mass = beads_mass,
      Segment = Segment, 
      min_cells = min_cells, max_bins = max_bins, step = step,
      excludedColumn = excludedColumn, 
      
      TIM = TIM, pathwayhierarchy = pathwayhierarchy,
      clustering.var = clustering.var, 
      cores = cores,
      save_processed_res = save_processed_res,
      savepath = savepath
    )
  }

  Ranking(MCprocess_assess_res, name = name, savepath = paste0(savepath, "/assess_res"))
}

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
    spillpath = NULL, spillname = NULL, FSC = "FSC-H", SSC = "SSC-H",
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
    lineara = 2, linearb = 0,
    Truncatea = 1,
    beads_mass = NULL,
    Segment = 200,
    Segment2 = 200,
    min_cells = 3, max_bins = 10, step = 10,
    excludedColumn = NULL,    
    save_processed_res = "one_folder",
    savepath = "./ANPELA_res",
    cores = floor(parallel::detectCores()/2)
){
  metadata <- paste0(datapath, "/metadata.csv")
  if (technique == "MC"){
    MCprocess_res <- MCprocess(name = name,
                               datapath = datapath, 
                               metadata = metadata,
                               technique = technique,
                               studytype = studytype, mergeM = mergeM, fixedNum = fixedNum, 
                               compensationM = compensationM, 
                               transformationM = transformationM,
                               normalizationM = normalizationM, 
                               signalcleanM = signalcleanM,
                               workflow = workflow,
                               single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass,
                               CATALYSTM = CATALYSTM,
                               sce_bead = sce_bead, marker_to_barc = marker_to_barc,
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
                               beads_mass = beads_mass,
                               Segment = Segment, 
                               min_cells = min_cells, max_bins = max_bins, step = step,
                               excludedColumn = excludedColumn, 
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)

    return(MCprocess_res)
  } else if (technique == "FC"){
    FCprocess_res <- FCprocess(name = name,
                               datapath = datapath, 
                               metadata = metadata,
                               technique = technique,
                               studytype = studytype, mergeM = mergeM, fixedNum = fixedNum,
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
    return(FCprocess_res)
  }
}


Assess <- function(
    name = "result", data = NULL, respath =NULL,
    studytype = c("CSI", "PTI"),  
    clusteringM = "FlowSOM",
    ncluster = 8,
    Phenograph_k = 30,
    DEP = NULL,
    
    TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh", "slingshot_tSNE","slingshot_FLOWMAP",
            "prinCurves_tSNE", "slingshot_PCA", "slingshot_diffMaps", "prinCurves_diffMaps"),
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
  spillpath = NULL, spillname = NULL, FSC = "FSC-H", SSC = "SSC-H",
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
  lineara = 2, linearb = 0,
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
  
  TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh", "slingshot_tSNE","slingshot_FLOWMAP",
          "prinCurves_tSNE", "slingshot_PCA", "slingshot_diffMaps", "prinCurves_diffMaps"),
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
                              # # 用一个全局变量防止重复 sink
                              # if (!exists(".log_started", envir = .GlobalEnv)) {
                              #   assign(".log_started", TRUE, envir = .GlobalEnv)
                              #   
                              #   pid <- Sys.getpid()
                              #   log_file <- paste0("log_worker_", pid, ".txt")
                              #   
                              #   # 尝试打开 sink
                              #   try({
                              #     sink(log_file, split = TRUE)
                              #     sink(log_file, type = "message", append = TRUE)
                              #   }, silent = TRUE)
                              # }
                              # 
                              # # 日志输出
                              # pid <- Sys.getpid()
                              # cat(sprintf("[%s] Worker PID %d started task %d\n", Sys.time(), pid, i))
                              # 
                                            #print(i)
                                            ########### start data processing ###########
                                            
                                            try(source("./processing.R"), silent = T)
                                            # set.seed(123)
                                            # if(file.exists(paste0(savepath, "/process_res/", rownames(workflow)[i], ".RData"))) {
                                            #   load(paste0(savepath, "/process_res/", rownames(workflow)[i], ".RData"))
                                            #   load(paste0(savepath, "/info_saved.RData"))
                                            #   dataFileNames <- info_saved$dataFileNames
                                            #   metadata  <-  info_saved$metadata
                                            #   index_TIclass  <-  info_saved$index_TIclass
                                          # } else {
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
                                            # return(res)
                                            
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
                                            #}
                                            ########### start workflow assessment ###########
                                            
                                            if (studytype == "CSI") {
                                              try(source("./CSI/1readfcs.R"), silent = T)
                                              try(source("./CSI/2cluster.R"), silent = T)
                                              try(source("./CSI/3criteria.R"), silent = T)
                                              try(source("./CSI/4plot.R"), silent = T)
                                              
                                              # AP2_processed_D_class
                                              # res <- data$AP2_pro1_frame_classTI[[i]]
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
                                              # return(sub_cluster_label)
                                              
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
                                              if ((is.null(marker_path)||is.null(known_celltype_path))& DEP == "" || is.null(DEP)) {
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
                                              #save(res, file = paste0(savepath, "/each_robfp_1e80_startclus/workflow",i,"_",dataset_name, "_assess.RData"))
                                              #cat(sprintf("[%s] Worker PID %d: Finished task %d\n", Sys.time(), pid, i))
                                              return(res)
                                            }
                                            
                                            
                                          }
  
  # # 关闭 sink（每个 worker 退出时关闭）
  # parallel::clusterEvalQ(cl, {
  #   try(sink(), silent = TRUE)
  #   try(sink(type = "message"), silent = TRUE)
  #   NULL
  # })
  
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
    
    table2["Correspondence"][table2["Correspondence"] == 1] <- 10
    table2["Correspondence"][table2["Correspondence"] < 1] <- 4
  }
  
  assess_res <- list(table = table, table2 = table2)
  if (!dir.exists(paste0(savepath, "/assess_res"))) {
    dir.create(paste0(savepath, "/assess_res"), recursive = T)
  }
  save(assess_res, file = paste0(savepath, "/assess_res/", name, "_assess.RData"))
  return(assess_res)
}

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



Get_CSIres <- function(
    respath, save_processed_res ="one_folder", workflow, savepath,
    marker_path,color =c("#F39B7FFF","#3C5488FF","#7E6148FF","#B09C85FF","#8491B4FF","#4DBBD5FF","#00A087FF",
                         "#91D1C2FF","#E64B35FF","grey80",RColorBrewer::brewer.pal(12, "Set3")),
    plot = c(T,F)){
  
  #try(source("./CSI/1readfcs.R"))
  try(source("./CSI/4plot.R"))
  
  CSIres_list <- list()
  
  #load data
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
  
  
  for ( i in 1:length(workflow)){  

    dataset_name <- workflow[i]
    if (save_processed_res == "one_folder"){
      load(grep(dataset_name, datapath, value = T))
      condition_info  <- info_saved[["metadata"]][["condition"]]
      colsToUse <- info_saved[["index_TIclass"]]
    } else if (save_processed_res == "one_RData") {
      index <- stringr::str_replace_all(data$index_TIclass, "\\(.*", "")
      res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      condition_info  <- data[["metadata"]][["condition"]]
      colsToUse <- data[["index_TIclass"]]
    }
    if (class(res) == "try-error") {
      print(paste0("Can not load data for dataset: ", dataset_name))
      next
    }
    
    merge_res <- merge_multi(multi_flowFrame = res,
                             condition_info = condition_info)
    data_j <- try(as.data.frame(merge_res$data, silent = T))
    condition <- merge_res$condition
    data_rd_j <- data_j[,colsToUse] # 仅使用 MCquan_res$index_TIclass 中的列
    data_rd_j$condition <- condition
    rm(merge_res,res)
    
    set.seed(123) # t-SNE降维的随机种子
    tsne_result <- BBmisc::suppressAll(Rtsne::Rtsne(data_rd_j[,!names(data_rd_j) %in% c("condition")], 
                                                    perplexity = 40, dims = 2, verbose = TRUE,
                                                    max_iter = 1000, check_duplicates = F, num_threads= 0)) #10
    
    tsne_df <- data.frame(
      tsne_x = tsne_result$Y[,1],
      tsne_y = tsne_result$Y[,2],
      condition = as.factor(data_rd_j$condition)
    )
    
    anno_res <- cell_annotation(data = data_rd_j[,!names(data_rd_j) %in% c("condition")], 
                                marker_path = marker_path, 
                                colsToUse = colsToUse)
    tsne_df$anno <- anno_res

    CSIres_list[[dataset_name]] <- tsne_df
    save(CSIres_list, file = paste0(savepath, "/", dataset_name,"_CSIres.RData"))
    if (plot){
      p <- ggplot(tsne_df, aes(x = tsne_x, y = tsne_y, color =  .data[["anno"]])) +
        geom_point(alpha = 0.6, size = 0.8) + 
        scale_color_manual(values = color,
                           breaks = unique(tsne_df[["anno"]]))+
        theme(
          legend.position = "right",
          panel.background = element_blank(),
          plot.background = element_blank(),  # 设置背景为透明
          panel.grid.major = element_blank(),  # 删除主网格线
          panel.grid.minor = element_blank(),  # 删除次要网格线
          axis.text = element_blank(),         # 删除坐标轴的文字
          axis.ticks = element_blank(),        # 删除坐标轴的刻度线
          axis.line = element_blank(),         # 删除坐标轴的线条
          panel.border = element_blank(),       # 删除边框
          axis.title = element_blank()         # 删除坐标轴的标题
        )
      
      ggsave(paste0(savepath, "/", dataset_name,"_CSI.png"), p, width = 10, height = 10, dpi = 300, bg="transparent")
    }
  }
}


Visualize <- function(
    studytype = c("CSI","PTI"),
    respath,  save_processed_res ="one_folder", workflow, savepath,studyname,
    plot_metric = c("Ca_metric","Cb_metric","Cc_metric","Cd_metric"),
    #CSI
    clusteringM = c("FlowSOM"), ncluster = 8,
    ntop = NULL, DEP = NULL,
    #PTI
    TIM = c("scorpius_distSpear", "scorpius_distPear","scorpius_distEucl", "scorpius_distManh", 
            "slingshot_tSNE","slingshot_FLOWMAP","slingshot_PCA", "slingshot_diffMaps"),
    clustering.var = NULL, pathwayhierarchy = NULL,
    Cc_metric = c("Spearman rank correlation", "Kendall rank correlation")
    ){
  
  #load data&metadata
  datapath <- list.files(paste0(respath, "/process_res/"), pattern = "\\.RData$", full.names = T)
  if (length(datapath)==0){
    stop("The parameter of 'respath' is incorrect. The 'datafile' cannot be loaded.")
  } 
  if (save_processed_res == "one_RData") {
    assign("data",load(datapath))
    data <- get(data)
    metadata <- data[["metadata"]]
    colsToUse <- data[["index_TIclass"]]
  } else if (save_processed_res == "one_folder") {
    info_saved <- try(load(paste0(respath, "/info_saved.RData")), silent = T)
    if (class(info_saved) == "try-error") {
      stop("The parameter of 'respath' is incorrect. The 'info_saved.RData' cannot be loaded.")
    } else if (info_saved == "info_saved"){
      load(paste0(respath, "/info_saved.RData"))
    }
    metadata <- info_saved[["metadata"]]
    colsToUse <- info_saved[["index_TIclass"]]
  }
  condition_info  <- metadata[["condition"]]
  
  
  if(!dir.exists(savepath)){
    dir.create(savepath, recursive = TRUE)
  }
  
  errorplot <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", x = 0.5, y = 0.5, 
             label = paste("During the plotting process,\n", 
                           "the program had some problems and will not be able to display the output."), 
             cex = 5, color = "black") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(), 
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(), 
          panel.border = ggplot2::element_blank())
  ################CSI Visualize#####################
  if (studytype == "CSI") {
    try(source("./CSI/1readfcs.R"))
    try(source("./CSI/2cluster.R"))
    try(source("./CSI/3criteria.R"))
    try(source("./CSI/4plot.R"))
    #arg_check
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
    
    # ntop
    if (missing(ntop)) {
      ntop <- floor(length(colsToUse)/2)
    } else if (ntop >= length(colsToUse) || ntop %% 1 != 0 || ntop <= 0) {
      stop("The value of 'ntop' is incorrect. It should be positive whole number and less than the number of your selected markers.")
    }
    
    
    
    # known_marker
    if (is.null(DEP)) {
      cat("*************************************************************************", "\n")
      cat("The standardized marker names are listed below. \n")
      cat("*************************************************************************", "\n")
      if (save_processed_res == "one_RData") {
        cat(paste0(colnames(data$AP2_pro1_frame_classTI[vapply(data$AP2_pro1_frame_classTI, Negate(is.null), NA)][[1]][[1]]@exprs), collapse = "\n"), "\n")
      } else if (save_processed_res == "one_folder") {
        load(datapath[1])
        cat(paste0(colnames(res[[1]]@exprs), collapse = "\n"), "\n")
        rm(res)
      }
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
    
    for ( w in 1:length(workflow)){  
      dataset_name <- workflow[w]
      if (save_processed_res == "one_folder"){
        load(grep(dataset_name, datapath, value = T))
      } else if (save_processed_res == "one_RData") {
        res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      }
      if (class(res) == "try-error") {
        print(paste0("Can not load data for dataset: ", dataset_name))
        next
      }
      
      AP2_processed_D_class <- res
      rm(res)
      
      
      # AP2_processed_data_class0
      data1 <- try(as.data.frame(readfcs_multi(fcsFiles = AP2_processed_D_class, mergeMethod = "all")), silent = T)
      if (class(data1) == "try-error") {
        print(paste0("Can not read data for dataset: ", dataset_name))
        next
      }
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
      AP2_processed_data_class <- AP2_processed_data_class0[, c(colsToUse, "filename", "condition")]
      rm(AP2_processed_data_class0)
      
      
      # cluster_label
      cluster_label <- try(data_cluster(data = AP2_processed_data_class[,1:(dim(AP2_processed_data_class)[2]-2)], 
                                        method = clusteringM, Phenograph_k = Phenograph_k, FlowSOM_k = ncluster, 
                                        FlowSeed = 40), silent = T)
      
      if (class(cluster_label) == "try-error") {
        print(paste0("Can not cluster data for dataset: ", dataset_name))
        next
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
      # return(sub_cluster_label)
      
      gc()
      
      #CSI_Ca_Plot####
      if ("Ca_metric" %in% plot_metric) {
        resAUC <- try(AUC(test = test_KNN$test, KNN_res = test_KNN$KNN_res), silent = T)
        case_problist <- resAUC$case_problist
        
        num_rows <- round(length(which(!sapply(test_KNN$test, is.null)))/2 + 1e-10)
        if (length(which(!sapply(test_KNN$test, is.null))) == 1) num_rows <- 2
        grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Ca.png"),bg = "white",
                       width = 10, height = 5*num_rows, res=300, units ="in")
        if(length(which(!sapply(test_KNN$test, is.null))) > 1) par(mfrow=c(num_rows, 2))
        CSI_Ca_plot <- list()
        for (j in which(!sapply(test_KNN$test, is.null))) {
          CSI_Ca_plot[[j]] <- ROC_plot(i = j, test = test_KNN$test,case_problist = case_problist )
          if (any(class(CSI_Ca_plot[[j]]) == "try-error")) { 
            print(errorplot) 
          }
        }
        dev.off()
        rm(resAUC, case_problist, CSI_Ca_plot)
      }
   
    
      #CSI_Cb_Plot####
      if ("Cb_metric" %in% plot_metric){
        CSI_Cb_plot <- try(cluster_plot(data_with_cluster = data_with_cluster), silent = T)
        if (any(class(CSI_Cb_plot) == "try-error")) {
          CSI_Cb_plot <- errorplot
        }
        grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Cb.png"),bg = "white",
                       width = 10, height = 15, res=300, units ="in")
        print(CSI_Cb_plot)
        dev.off()
        rm(CSI_Cb_plot)
      }
    
      #CSI_Cc_plot####
      if ("Cc_metric" %in% plot_metric) {
        #Venn
        CS_preres <- try(CS_pre(test_KNN$subdata_cluster), silent = T)
        resCS <- try(CWfun(CS_preres = CS_preres, top = ntop), silent = T)
      
        CSI_Cc_plot_Venn <- list()
        a <- as.numeric(which(!sapply(resCS$DEGlist, is.null)))
        if (length(a) == 0) {
          message("There is no marker in the list of differentially expressed genes. 
                   The program will not be able to display the Cc_metric output.")
        } else {
          num_rows <- round(length(a)/2 + 1e-10)
          if (length(a) == 1) num_rows <- 2
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Cc_Venn.png"),bg = "white",
                         width = 10, height = 5*num_rows, res=300, units ="in")
          for (j in a) {
            CSI_Cc_plot_Venn[[j]] <- try(venn_plot(i = j, CConsistency_marker = resCS$DEGlist), silent = T)
            if (length(CSI_Cc_plot_Venn) != 0 && any(class(CSI_Cc_plot_Venn[[j]]) == "try-error")) {
              CSI_Cc_plot_Venn[[j]] <- errorplot
            } 
          }
          if (length(a) > 1){
            gridExtra::grid.arrange(grobs = CSI_Cc_plot_Venn, nrow = num_rows, ncol = 2)
          } else {
            gridExtra::grid.arrange(grobs = CSI_Cc_plot_Venn, nrow = 1, ncol = 1)
          }
          dev.off()
          rm(CS_preres, resCS, CSI_Cc_plot_Venn, a)
        }  
        #Volcano
        volcanodata <- lapply(test_KNN$subdata_cluster, volcanovalue, control = as.character(unique(metadata$condition)[1]), 
                              case = as.character(unique(metadata$condition)[2]))
        
        CSI_Cc_plot_Volcano <- list()
        non_null_indices <- as.numeric(which(!sapply(volcanodata, is.null)))
        combined_plot <- NULL
        for (j in non_null_indices) {
          CSI_Cc_plot_Volcano[[j]] <- try(volcano_plot(i = j, volcanodata = volcanodata), silent = T)
          if (any(class(CSI_Cc_plot_Volcano[[j]]) == "try-error")) {
            CSI_Cc_plot_Volcano[[j]] <- errorplot
          } 
        
          if (is.null(combined_plot)) {
            combined_plot <- CSI_Cc_plot_Volcano[[j]]
          } else {
            combined_plot <- combined_plot + CSI_Cc_plot_Volcano[[j]]
          }
        }
        
        # 根据非空元素数量设置布局
        if (length(non_null_indices) > 1) {
          num_rows <- round(length(non_null_indices) / 2+ 1e-10)
          combined_plot <- combined_plot + patchwork::plot_layout(nrow = num_rows, ncol = 2)
        } else {
          num_rows <- 2
        }
        
        # 打开 PNG 设备并保存组合后的图形
        grDevices::png(paste0(savepath, "/", studyname, "_", dataset_name, "_CSI_Cc_Volcano.png"), 
                       bg = "white",
                       width = 12, 
                       height = 5 * num_rows, 
                       res = 300, 
                       units = "in")
        print(combined_plot)
        dev.off()
        rm(volcanodata, CSI_Cc_plot_Volcano, combined_plot)
      }
      
      #CSI_Cd_plot####
      if("Cd_metric" %in% plot_metric){
        if (DEP == "" || is.null(DEP)) {
          print("The parameter 'DEP' is NULL. Please input the known biomarker(s) differentially expressed between two conditions.")
          next
        } else {
          known_marker <- unlist(strsplit(DEP, "\\s*,\\s*"))
          if (length(known_marker) != 0) {
            if (length(known_marker) > 8) {
              widthE_plot <- ceiling(length(known_marker)) * 1.4
            } else widthE_plot <- 12
            
            CSI_Cd_plot <- try(recall_plot(data_with_cluster = data_with_cluster, known_marker = known_marker), silent = T)
            if (any(class(CSI_Cd_plot) == "try-error")) {
              CSI_Cd_plot <- errorplot
            }
            grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_CSI_Cd.png"),bg = "white",
                           width = widthE_plot, height = 7, res=300, units ="in")
            print(CSI_Cd_plot)
            dev.off()
            rm(CSI_Cd_plot)
          }
        }
      }
      rm(data_with_cluster, test_KNN)
      gc()
    }
  ###########PTI Visualize#####################
  } else if (studytype == "PTI") {
    try(source("./PTI/load_data2.R"))
    try(source("./PTI/TI_method.R"))
    
    try(source("./PTI/Bio_con_4.R"))
    try(source("./PTI/time_metric_3.R"))
    try(source("./PTI/Robustness_4.R"))
    try(source("./PTI/Rough_3.R"))
    
    try(source("./PTI/shift_start.R"))
    try(source("./PTI/calc_spline.R"))
    try(source("./PTI/cycle_pseudotime.R"))
    try(source("./PTI/reverse_pseudotime.R"))
    try(source("./PTI/check_pairs.R"))
    
    try(source("./PTI/plot.R"))

    
    # TIM
    if (missing(TIM)) {
      TIM <- "scorpius_distSpear"
    } else {
      TIM <- match.arg(TIM)
    }
    
    # Cc_metric
    if (missing(Cc_metric)) {
      Cc_metric <- "Spearman rank correlation"
    } else {
      Cc_metric <- match.arg(Cc_metric)
    }
    #index
    index <- stringr::str_replace_all(colsToUse, "\\(.*", "")
    
    for (w in 1:length(workflow)){
      dataset_name <- workflow[w]
      
      #load data
      if (save_processed_res == "one_folder"){
        load(grep(dataset_name, datapath, value = T))
      } else if (save_processed_res == "one_RData") {
        res <- data[["AP2_pro1_frame_classTI"]][[dataset_name]]
      }
      
      AP2_processed_D_TI <- try(load.Data(res, index = index, measurement.time = as.matrix(metadata$timepoint), TIM = TIM), silent = T)
      if (class(res) == "try-error") {
        print(paste0("Can not load data for dataset: ", dataset_name))
        next
      }
      
      TIres <- try(TI(D = AP2_processed_D_TI, method = TIM, 
                        dataset_name = dataset_name, clustering.var = clustering.var), silent = T)

      if (class(TIres) == "try-error") {
        print(paste0( "Error in the pseudo-time trajectory inference for dataset: ", dataset_name))
        next
      }
      

      #PTI_Ca_Plot####
      #time
      if ("Ca_metric" %in% plot_metric) {
        grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Ca_time.png"),
                       bg = "white",width = 10, height = 10, res=300, units ="in")
        if (!grepl("slingshot", TIM)){
          PTI_Ca_Plot_time <- try(dr(TIres, AP2_processed_D_TI), silent = T)
          if (any(class(PTI_Ca_Plot_time) == "try-error")) {
            print(errorplot)
          } else {
            print(PTI_Ca_Plot_time) 
          }
        } else {
          PTI_Ca_Plot_time <- try(dr_multi(TIres, AP2_processed_D_TI), silent = T)
        }
        if (any(class(PTI_Ca_Plot_time) == "try-error")) {
          print(errorplot)
        }
        dev.off()
        rm(PTI_Ca_Plot_time)
      }
      #abunf_pt
      result <- TIres
      if (is.list(result$trajectory)) {
        colnames(result$pseudotime) <- paste0("Trajectory_", seq_len(ncol(result$pseudotime)))
        names(result$trajectory) <- paste0("Trajectory_", seq_len(ncol(result$pseudotime)))
        colnames(result$crv1) <- paste0("Trajectory_", seq_len(ncol(result$pseudotime)))
      }
      for(traj_idx in 1:TIres$lineages){
        if (TIres[["cr_method"]] == "run_slingshot"){
          result$pseudotime <- TIres$pseudotime[, traj_idx]
          result$trajectory <- TIres$trajectory[[paste0("Trajectory_", traj_idx)]]
          result$linInd <- as.numeric(traj_idx)
        }
        bio_meaning <- try(suppressWarnings(abund_pt_plot(TIres = result, D = AP2_processed_D_TI)), silent = T)
        
        if ("Ca_metric" %in% plot_metric) {
          heightA2_plot <- ceiling(length(unique(metadata$timepoint))/4) * 3.5
          widthA2_plot <- 12
          PTI_Ca_Plot_prot <- list()
          combined_plot <- NULL
          for (i in 1:length(index)) {
            PTI_Ca_Plot_prot[[i]] <- try(abund_pt_single(bio_meaning$to_plot, index[i], bio_meaning$dat), silent = T)
           if (any(class(PTI_Ca_Plot_prot[[i]]) == "try-error")) {
             PTI_Ca_Plot_prot[[i]] <- errorplot
           }
            if (is.null(combined_plot)) {
              combined_plot <- PTI_Ca_Plot_prot[[i]]
            } else {
              combined_plot <- combined_plot + PTI_Ca_Plot_prot[[i]]
            }
          }
          if (length(index) > 1) {
            num_rows <- round(length(index) / 2+ 1e-10)
            combined_plot <- combined_plot + patchwork::plot_layout(nrow = num_rows, ncol = 2)
          } else {
            num_rows <- 2
          }
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Ca_prot_trajectory",traj_idx,".png"),
                         bg = "white",width = widthA2_plot*2, height = heightA2_plot*num_rows, res=300, units ="in")
          print(combined_plot)
          dev.off()
          rm(PTI_Ca_Plot_prot, combined_plot)
        }
        
        # PTI_Cb_plot #####
        if ("Cb_metric" %in% plot_metric) {
          heightB_plot <- ceiling(ncol(AP2_processed_D_TI$expr)/4) * 2.75
          widthB_plot <- 12
          PTI_Cb_plot <- try(roughness_plot(result, AP2_processed_D_TI), silent = T)
          if (any(class(PTI_Cb_plot) == "try-error")) {
            PTI_Cb_plot <- errorplot
          }
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Cb_trajectory",traj_idx,".png"),
                         bg = "white",width = widthB_plot, height = heightB_plot, res=300, units ="in")
          print(PTI_Cb_plot)
          dev.off()
          rm(PTI_Cb_Plot)
        }
        
          # PTI_Cc_plot ####
        if ("Cc_metric" %in% plot_metric) {
          Rob <- try(Robustness(TIres = result, D = AP2_processed_D_TI, nruns = 4, cell.subset = 0.8,
                                clustering.var = clustering.var, dataset_name =dataset_name), silent = T)
          
          PTI_Cc_plot <- try(robustness_new_plot(input_matrix = Rob$input_matrix, 
                                                 finalMatrix = Rob$finalMatrix, method =Cc_metric), silent = T)
          if (any(class(PTI_Cc_plot) == "try-error")) {
            PTI_Cc_plot <- errorplot
          }
        
          grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Cc_trajectory",traj_idx,".png"),
                         bg = "white",width = 10, height = 10, res=300, units ="in")
          print(PTI_Cc_plot) 
          dev.off()
          rm(PTI_Cc_plot,Rob)
        }
        
        # PTI_Cd_plot ####
        if ("Cd_metric" %in% plot_metric) {
          if (!is.null(pathwayhierarchy) && file.exists(pathwayhierarchy)) {
           heightD1_plot <- ceiling(ncol(AP2_processed_D_TI$expr)/4) * 2.75
           
           PTI_Cd_plot <- try(bio_meaning$p3, silent = T)
           if (any(class(PTI_Cd_plot) == "try-error")) {
             PTI_Cd_plot <- errorplot
           }
           
           grDevices::png(paste0(savepath, "/",studyname, "_", dataset_name,"_PTI_Cd_trajectory",traj_idx,".png"),
                          bg = "white",width = 12, height = heightD1_plot, res=300, units ="in")
           print(PTI_Cd_plot) 
           dev.off()
           rm(PTI_Cd_plot)
          }
        }
        rm(bio_meaning)
      }
      
      rm(AP2_processed_D_TI, TIres)
      gc()
    }
  }
}