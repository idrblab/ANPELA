#' @title Process
#'
#' @param name
#' @param datapath
#' @param techique
#' @param studytype
#' @param mergeM
#' @param fixedNum
#' @param compensationM
#' @param transformationM
#' @param normalizationM
#' @param signalcleanM
#' @param spillpath
#' @param FSC
#' @param SSC
#' @param control.dir
#' @param control.def.file
#' @param single_pos_fcs
#' @param single_pos_mass
#' @param CATALYSTM
#' @param beads_mass
#' @param index_protein
#' @param save_processed_res
#' @param savepath
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' }

Process <- function(
    name,
    datapath,
    techique = c("MC", "FC"),
    studytype = c("CSI", "PTI"),
    mergeM = c("Fixed", "Ceil", "All", "Min"),
    fixedNum = 200,
    compensationM = c("AutoSpill", "CATALYST", "CytoSpill", "FlowCore", "MetaCyto", "None"),
    transformationM = c("Arcsinh Transformation", "Asinh with Non-negative Value", "Asinh with Randomized Negative Value",
                        "Biexponential Transformation", "Box-Cox Transformation", "FlowVS Transformation", "Hyperlog Transformation", "Linear Transformation",
                        "Ln Transformation", "Log Transformation", "Logicle Transformation", "Quadratic Transformation", "Scale Transformation", "Truncate Transformation",
                        "None"),
    normalizationM = c("Bead-based Normalization", "GaussNorm", "WarpSet", "ZScore", "None"),
    signalcleanM = c("FlowAI", "FlowClean", "FlowCut", "None"),
    spillpath = NULL, FSC = "FSC-H", SSC = "SSC-H",
    control.dir = NULL, control.def.file = NULL,
    single_pos_fcs = NULL, single_pos_mass = NULL, CATALYSTM = c("flow", "nnls"),
    beads_mass = c(140, 151, 153, 165, 175),
    index_protein = NULL,
    save_processed_res = c("no", "one_folder", "one_RData"),
    savepath = "./",
    cores = floor(parallel::detectCores()/2)
){
  metadata <- paste0(datapath, "/metadata.csv")
  if (techique == "MC"){
    MCprocess_res <- MCprocess(datapath = datapath,
                               metadata = metadata,
                               studytype = studytype, mergeM = mergeM, fixedNum = fixedNum,
                               compensationM = compensationM,
                               transformationM = transformationM,
                               normalizationM = normalizationM,
                               signalcleanM = signalcleanM,
                               single_pos_fcs = single_pos_fcs, single_pos_mass = single_pos_mass, CATALYSTM = "nnls",
                               beads_mass = beads_mass,
                               index_protein = index_protein,
                               save_processed_res = save_processed_res,
                               savepath = savepath,
                               cores = cores)
    if (save_processed_res == "one_RData"){
      if (!dir.exists(paste0(savepath, "/process_res"))) {
        dir.create(paste0(savepath, "/process_res"), recursive = T)
      }
      save(MCprocess_res, file = paste0(savepath, "/process_res/", name, ".RData"))
    }
    return(MCprocess_res)
  } else if (techique == "FC"){
    FCprocess_res <- FCprocess(datapath = datapath,
                               metadata = metadata,
                               studytype = studytype, mergeM = mergeM, fixedNum = fixedNum,
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
    if (save_processed_res == "one_RData"){
      if (!dir.exists(paste0(savepath, "/process_res"))) {
        dir.create(paste0(savepath, "/process_res"), recursive = T)
      }
      save(FCprocess_res, file = paste0(savepath, "/process_res/", name, ".RData"))
    }
    return(FCprocess_res)
  }
}

