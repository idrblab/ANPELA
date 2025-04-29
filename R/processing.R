comp_anpela <- function(data, method, index,
                        spillpath, spillname, FSC,  SSC,
                        control.dir, control.def.file,
                        single_pos_fcs, single_pos_mass, CATALYSTM,
                        sce_bead, marker_to_barc) {
  res <- switch (method,
                 "FlowCore" = flowCore_comp(frame_list = data, col_names = index, spillpath = spillpath,
                                            spillname = spillname, FSC = FSC, SSC = SSC),
                 "AutoSpill" = autospill_comp(frame_list = data, col_names = index, control.dir = control.dir, control.def.file = control.def.file),
                 "CytoSpill" = CytoSpill_comp(frame_list = data, cols = index),
                 "CATALYST" = CATALYST_comp(frame_list = data, cols = index,
                                            single_pos_fcs = single_pos_fcs,
                                            single_pos_mass = single_pos_mass,
                                            method = CATALYSTM),
                 "MetaCyto" = MetaCyto_comp(frame_list = data, col_names = index),
                 "spillR" = spillR_comp(frame_list = data, col_names = index,
                                        sce_bead = sce_bead, marker_to_barc = marker_to_barc),
                 None = data
  )
  res <- lapply(res, function(x) {
    x@exprs[is.infinite(x@exprs)] <- NA
    x@exprs[is.nan(x@exprs)] <- NA
    return(x)
  })
  names(res) <- names(data)
  return(res)
}


trans_anpela <- function(data, method, index,
                         arcsinha, arcsinhb, arcsinhc,
                         anna,annb, annc, annthreshold,
                         arna, arnb, arnc, arnthreshold,
                         bepa, bepb, bepc, bepd, bepf, bepw, tol, maxit,
                         hpla, hplb,
                         lineara, linearb,
                         lntr, lntd,
                         logbase,logr,logd,
                         lgtw, lgtt, lgtm, lgta,
                         Quadratica, Quadraticb, Quadraticc,
                         Truncatea) {
  res <- switch (method,
                 "Arcsinh Transformation" = arcsinh_trans(frame_list = data, col_names = index,a = arcsinha, b = arcsinhb ,c = arcsinhc),
                 "Asinh with Non-negative Value" = ANN_trans(frame_list = data, col_names = index,a = anna, b = annb, c = annc, threshold = annthreshold),
                 "Asinh with Randomized Negative Value" = ARN_trans(frame_list = data, col_names = index, a=arna, b = arnb, c = arnc, threshold = arnthreshold),
                 "Biexponential Transformation" = biexp_trans(frame_list = data, col_names = index,  a = bepa, b = bepb, c = bepc, d = bepd, f = bepf, w = bepw,
                                                              tol = tol, maxit = maxit),
                 "Box-Cox Transformation" = BoxCox_trans(frame_list = data, col_names = index),
                 "FlowVS Transformation" = flowVS_trans(frame_list = data, col_names = index),
                 "Hyperlog Transformation" = hyperlog_trans(frame_list = data, col_names = index,  a = hpla, b = hplb),
                 "Linear Transformation" = linear_trans(frame_list = data, col_names = index, a = lineara, b = linearb),
                 "Ln Transformation" = ln_trans(frame_list = data, col_names = index, r = lntr, d = lntd),
                 "Log Transformation" = log_trans(frame_list = data, col_names = index, logbase = as.numeric(logbase), r = logr, d = logd),
                 "Logicle Transformation" = logicle_trans(frame_list = data, col_names = index,  w = lgtw, t = lgtt, m = lgtm, a = lgta),
                 "Quadratic Transformation" = quadratic_trans(frame_list = data, col_names = index,
                                                              a = Quadratica, b = Quadraticb, c = Quadraticc),
                 "Split Scale Transformation" = splitScale_trans(frame_list = data, col_names = index),
                 "Truncate Transformation" = truncate_trans(frame_list = data, col_names = index, a = Truncatea),
                 #"Adaptive Box-Cox Transformation" = AdaptiveBoxCox_trans(frame_list = data, col_names = index),
                 "Centered Log Ratio Transformation" = CLR_trans(frame_list = data, col_names = index),
                 "None" = data
  )
  res <- lapply(res, function(x) {
    x@exprs[is.infinite(x@exprs)] <- NA
    x@exprs[is.nan(x@exprs)] <- NA
    return(x)
  })
  names(res) <- names(data)
  return(res)
}


norm_anpela <- function(data, method, index, beads_mass) {
  res <- switch (method,
                 "Bead-based Normalization" = CATALYST_norm(frame_list = data, col_names = index, beads_mass = beads_mass),
                 "GaussNorm" = gaussNorm_norm(frame_list = data, col_names = index),
                 "WarpSet" = warpSet_norm(frame_list = data, col_names = index),
                 "ZScore" = ZScore_norm(frame_list = data, col_names = index),
                 "Mean Normalization" = Mean_norm(frame_list = data, col_names = index),
                 "Min-max Normalization" = MinMax_norm(frame_list = data, col_names = index),
                 "None" = data
  )
  res <- lapply(res, function(x) {
    x@exprs[is.infinite(x@exprs)] <- NA
    x@exprs[is.nan(x@exprs)] <- NA
    return(x)
  })
  names(res) <- names(data)
  return(res)
}

sigcl_anpela <- function(data, method, index,
                         Segment,
                         Segment2,
                         min_cells, max_bins, step, technique) {
  res <- switch (method,
                 "FlowAI" = flowAI_signalC(frame_list = data, col_names = index),
                 "FlowCut" = flowCut_signalC(frame_list = data, Segment = Segment, col_names = index),
                 "FlowClean" = flowClean_signalC(frame_list = data, col_names = index, Segment = Segment2),
                 "PeacoQC" = PeacoQC_signalC(frame_list = data, col_names = index,
                                             min_cells = min_cells, max_bins = max_bins, step = step, technique = technique),
                 "None" = data
  )
  res <- lapply(res, function(x) {
    x@exprs[is.infinite(x@exprs)] <- NA
    x@exprs[is.nan(x@exprs)] <- NA
    return(x)
  })
  names(res) <- names(data)
  return(res)
}

# Compensation ------------------------------------------------------------

# Compensation: flowCore ------------------------------------------------------------
flowCore_comp <- function(frame_list, col_names, spillpath = NULL, spillname = NULL, FSC = "FSC-H", SSC = "SSC-H") {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)
  col_names1 <- make.names(col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })
  if (!is.null(spillpath) && !is.null(spillname)) {
    frames <- lapply(spillpath, flowCore::read.FCS)
    names(frames) <- spillname
    frames <- as(frames, "flowSet")
    spill_single <- flowCore::spillover(frames, unstained="unstained",
                                        fsc = FSC, ssc = SSC,
                                        stain_match = "regexpr", useNormFilt = FALSE)
    spill <- lapply(1:length(frame_list1), function(i) return(spill_single))
  } else if (all(sapply(frame_list1, function(x) !is.null(x@description[["SPILL"]])))) {
    spill <- lapply(frame_list1, function(x) return(x@description[["SPILL"]]))
  } else if (all(sapply(frame_list1, function(x) !is.null(x@description[["SPILLOVER"]])))) {
    spill <- lapply(frame_list1, function(x) return(x@description[["SPILLOVER"]]))
  } else if (all(sapply(frame_list1, function(x) !is.null(x@description[["spillover"]])))) {
    spill <- lapply(frame_list1, function(x) return(x@description[["spillover"]]))
  }

  x_comp <- lapply(1:length(frame_list1), function (i) {
    dimnames(spill[[i]])[[1]] <- make.names(dimnames(spill[[i]])[[1]])
    dimnames(spill[[i]])[[2]] <- make.names(dimnames(spill[[i]])[[2]])
    return(flowCore::compensate(x = frame_list1[[i]], spillover = spill[[i]]))
  })

  for (j in 1:length(frame_list1)) {
    frame_list1[[j]]@exprs[, col_names1] <- x_comp[[j]]@exprs[, col_names1]
  }

  frame_list1 <- lapply(frame_list1, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(frame_list1)
}

# Compensation: CytoSpill ------------------------------------------------------------
CytoSpill_comp <- function(frame_list, cols) {
  GetSpillMat <- function (data = NULL, cols, n, file = NULL, threshold = 0.1) {

    cutoffs <- CytoSpill:::.DeriveCutoffs(data, cols, n)
    model <- .EstimateSpill(data, cutoffs, cols, upperbound = threshold)
    estimates <- model[[1]]
    xcols <- model[[2]]
    spillmat <- diag(length(xcols))
    for (i in 1:length(xcols)) {
      if (!all(is.na(xcols[[i]]))) {
        for (j in 1:length(xcols[[i]])) {
          spillmat[xcols[[i]][j], i] <- ifelse(estimates[[i]][j] <
                                                 threshold, estimates[[i]][j], threshold)
        }
      }
    }
    return(spillmat)
  }

  .EstimateSpill <- function (data, cutoffs, cols, upperbound = 0.1) {
    results <- list()
    data <- data[, cols]
    xcols <- CytoSpill:::.GetFmla(data, spill_cols = .SpillColsData(data))
    for (i in 1:ncol(data)) {
      if (!is.na(xcols[[i]][1])) {
        A = as.matrix(data[which(data[, i] < cutoffs[i]),
                           xcols[[i]]])
        b = data[which(data[, i] < cutoffs[i]), i]
        x0 = runif(ncol(A), min = 0, max = upperbound)
        fn = function(x) {
          vec = A %*% x - b
          norm(vec, type = "2")
        }
        result = try(nloptr::slsqp(x0, fn, lower = rep(0,length(x0)), upper = rep(upperbound, length(x0))))
        if (isTRUE(class(result) == "try-error")) {
          result <- NULL
          xcols[[i]] <- NA
        }
        else {
          result <- result$par
        }
        results[[i]] <- result
      }
      else {
        results[[i]] <- NULL
      }
    }
    return(list(results, xcols))
  }

  .SpillColsData <- function (data, l = CATALYST::isotope_list) {
    chs <- unlist(regmatches(colnames(data), gregexpr("\\(.*\\)", colnames(data))))
    ms <- as.numeric(sapply(stringr::str_extract_all(colnames(data),"[0-9]+"), function(x) {x[length(x)]}))
    # ms <- as.numeric(regmatches(chs, gregexpr("[0-9]+", chs)))
    mets <- sapply(stringr::str_extract_all(colnames(data),"[a-zA-Z]+"), function(x) {x[length(x)-1]})
    # mets <- gsub("[[:digit:]]+Di", "", chs)
    spill_cols <- vector("list", length(ms))
    for (i in seq_along(ms)) {
      p1 <- p2 <- m1 <- m2 <- ox <- iso <- NULL
      if ((ms[i] + 1) %in% ms)
        p1 <- which(ms == (ms[i] + 1))

      if ((ms[i] + 2) %in% ms)
        p2 <- which(ms == (ms[i] + 2))

      if ((ms[i] - 1) %in% ms)
        m1 <- which(ms == (ms[i] - 1))
      if ((ms[i] - 2) %in% ms)
        m2 <- which(ms == (ms[i] - 2))
      if ((ms[i] + 16) %in% ms)
        ox <- which(ms == (ms[i] + 16))
      zy <- try({
        iso <- l[[mets[i]]]
        iso <- which(ms %in% iso[iso != ms[i]])
      }, silent = T)
      if (class(zy) == "try-error") {
        spill_cols[[i]] <- unique(c(m1, m2, p1, p2, ox))
      } else {
        spill_cols[[i]] <- unique(c(m1, m2, p1, p2, iso, ox))
      }
    }
    return(spill_cols)
  }

  res <- list()
  for (i in 1:length(frame_list)) {
    fcs_exprs <- frame_list[[i]]@exprs
    set.seed(123)
    spillmat <- suppressWarnings(GetSpillMat(fcs_exprs, cols, n = nrow(fcs_exprs)))
    data_compensated <- t(apply(fcs_exprs[, cols], 1, function(row) nnls::nnls(t(spillmat),row)$x))
    data_colnames <- colnames(fcs_exprs)
    fcs_exprs[, cols] <- data_compensated
    colnames(fcs_exprs) <- data_colnames
    res[[i]] <- frame_list[[i]]
    res[[i]]@exprs <- fcs_exprs
  }
  return(res)
}

# Compensation: CATALYST ------------------------------------------------------------
CATALYST_comp <- function(frame_list, cols, single_pos_fcs, single_pos_mass, method = c("flow", "nnls")) {
  col_names1 <- stringr::str_extract(cols, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)
  cols <- col_names1
  frame_list1 <- frame_list
  for (i in 1:length(frame_list1)) {
    colnames(frame_list1[[i]]@exprs) <- frame_list[[i]]@parameters@data$name
  }

  spill <- CATALYST::assignPrelim(x = single_pos_fcs, y = single_pos_mass) %>%
    CATALYST::estCutoffs() %>%
    CATALYST::applyCutoffs() %>%
    CATALYST::computeSpillmat()

  dimnames(spill)[[1]] <- make.names(dimnames(spill)[[1]])
  dimnames(spill)[[2]] <- make.names(dimnames(spill)[[2]])

  res <- list()
  for (i in 1:length(frame_list1)) {
    fcs_exprs <- frame_list1[[i]]
    fcs_exprs@exprs <- fcs_exprs@exprs[, cols]
    data_compensated <- switch (method,
                                flow = CATALYST::compCytof(x = fcs_exprs, y = spill, method = "flow"),
                                nnls = CATALYST::compCytof(x = fcs_exprs, y = spill, method = "nnls")
    )
    frame_list1[[i]]@exprs[, cols] <- data_compensated@exprs
    colnames(frame_list1[[i]]@exprs) <- colnames(frame_list[[i]]@exprs)
  }
  return(frame_list1)
}

# Compensation: AutoSpill ------------------------------------------------------------
autospill_comp <- function(frame_list, col_names, control.dir, control.def.file) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  # set parameters
  param <- autospill::get.autospill.param("minimal")
  # adjusted parameters
  param$worker.process.n <- 1
  # read flow controls
  flow.control <- autospill::read.flow.control(control.dir, control.def.file, param)
  # gate events before calculating spillover
  flow.gate <- autospill::gate.flow.data(flow.control, param) # case data 31.86s
  # get initial spillover matrices from untransformed data
  marker.spillover.unco.untr <- autospill::get.marker.spillover(TRUE, flow.gate, flow.control, param)
  # refine spillover matrix iteratively
  invisible(capture.output(
    refine.spillover.result <- autospill::refine.spillover(marker.spillover.unco.untr, NULL, flow.gate, flow.control, param)
  ))

  spill <- refine.spillover.result[["spillover"]]
  dimnames(spill)[[1]] <- make.names(dimnames(spill)[[1]])
  dimnames(spill)[[2]] <- make.names(dimnames(spill)[[2]])
  x_comp <- lapply(frame_list1, flowCore::compensate, spillover = spill)
  for (j in 1:length(frame_list1)) {
    frame_list1[[j]]@exprs[, col_names1] <- x_comp[[j]]@exprs[, col_names1]
  }
  frame_list1 <- lapply(frame_list1, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(frame_list1)
}

# Compensation: MetaCyto ------------------------------------------------------------
MetaCyto_comp <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  if (!is.null(flowCore::keyword(frame_list1[[1]], "SPILL")[[1]])) {
    check <- lapply(frame_list1, function(x) {
      result = isSymmetric(flowCore::keyword(x, "SPILL")[[1]])
    })
    if (unique(check)[[1]] == FALSE) {
      x_comp = lapply(frame_list1, function(x) {
        new_frame = flowCore::compensate(x, flowCore::keyword(x, "SPILL")[[1]])
        return(new_frame)
      })
    }
  }

  for (j in 1:length(x_comp)) {
    frame_list1[[j]]@exprs[, col_names1] <- x_comp[[j]]@exprs[, col_names1]
  }

  frame_list1 <- lapply(frame_list1, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(frame_list1)
}
# Compensation: spillR ------------------------------------------------------------
spillR_comp <- function(frame_list, col_names, sce_bead, marker_to_barc) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  # sce_bead <- CATALYST::prepData(ss_exp)
  # sce_bead <- CATALYST::assignPrelim(sce_bead, bc_key, verbose = FALSE)
  # sce_bead <- CATALYST::applyCutoffs(estCutoffs(sce_bead))
  # sce_bead <- CATALYST::computeSpillmat(sce_bead)
  #
  # marker_to_barc <- SummarizedExperiment::rowData(sce_bead)[, c("channel_name", "is_bc")] |>
  #   dplyr::as_tibble() |>
  #   dplyr::filter(is_bc == TRUE) |>
  #   dplyr::mutate(barcode = bc_key) |>
  #   dplyr::select(marker = channel_name, barcode)

  for (i in 1:length(frame_list1)) {
    fcs_exprs <- frame_list1[[i]]
    fcs_exprs@exprs <- fcs_exprs@exprs[, col_names1]
    sce <- CATALYST::prepData(fcs_exprs)
    data_compensated <- spillR::compCytof(sce, sce_bead, marker_to_barc, impute_value = NA)
    frame_list1[[i]]@exprs[, col_names1] <- t(data_compensated@assays@data@listData[["compcounts"]])
    colnames(frame_list1[[i]]@exprs) <- colnames(frame_list[[i]]@exprs)
    rownames(frame_list1[[i]]@exprs) <- rownames(frame_list[[j]]@exprs)
  }
  return(frame_list1)

}

# Transformation ------------------------------------------------------------

# Transformation: arcsinh transformation ------------------------------------------------------------
arcsinh_trans <- function(frame_list, col_names, a = 0, b = 1, c = 0) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })
  asinhTrans <- flowCore::arcsinhTransform(a = a, b = b, c = c)
  translist <- flowCore::transformList(col_names1, asinhTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: asinh with non-negative value ------------------------------------------------------------
ANN_trans <- function(frame_list, col_names, a = 0, b = 1, c = 0, threshold = 1) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <-lapply(frame_list, function(x, col_names1, threshold){
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))

    data <- x@exprs[, col_names1] - threshold
    data[data < 0] <- 0
    x@exprs[, col_names1] <- data
    return(x)
  }, col_names1 = col_names1, threshold = threshold)
  asinhTrans <- flowCore::arcsinhTransform(a = a, b = b, c = c)
  translist <- flowCore::transformList(col_names1, asinhTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: asinh with randomized negative value ------------------------------------------------------------
ARN_trans <- function(frame_list, col_names, a = 0, b = 1, c = 0, threshold = 1) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <-lapply(frame_list, function(x, col_names1, threshold){
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))

    data <- x@exprs[, col_names1] - threshold
    set.seed(123)
    data[data < 0] <- rnorm(length(data[data < 0]))
    x@exprs[, col_names1] <- data
    return(x)
  }, col_names1 = col_names1, threshold = threshold)

  asinhTrans <- flowCore::arcsinhTransform(a = a, b = b, c = c)
  translist <- flowCore::transformList(col_names1, asinhTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: Biexponential transformation ------------------------------------------------------------
biexp_trans <- function(frame_list, col_names, a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0,
                        tol = .Machine$double.eps^0.25, maxit = as.integer(5000)) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  biexpTrans <- flowCore::biexponentialTransform(a = a, b = b, c = c, d = d, f = f, w = w,
                                                 tol = tol, maxit = maxit)
  translist <- flowCore::transformList(col_names1, biexpTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: Box-Cox transformation ------------------------------------------------------------
# BoxCox_trans <- function(frame_list, col_names, lambda = 0.3) {
#   col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
#   col_names1 <- sub("\\(", "", col_names1)
#   col_names1 <- sub("\\)", "", col_names1)
#
#   frame_list1 <- lapply(frame_list, function(x) {
#     colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
#     colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
#     colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
#     return(x)
#   })
#
#   dataTransform <- lapply(frame_list1, function(x){
#     x@exprs[,col_names1] <- flowClust::box(x@exprs[,col_names1], lambda = lambda)
#     return(x)
#   })
#
#   dataTransform <- lapply(dataTransform, function(x) {
#     colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
#     return(x)
#   })
#   return(dataTransform)
# }
BoxCox_trans <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  dataTransform_pre <- lapply(frame_list1, flowTrans::flowTrans,
                              fun = "mclMultivBoxCox", dims = col_names1,
                              n2f = F, parameters.only = F)
  dataTransform <- frame_list1
  for (j in 1:length(dataTransform)) {
    dataTransform[[j]]@exprs <- dataTransform_pre[[j]][["result"]]@exprs
    rownames(dataTransform[[j]]@exprs) <- rownames(frame_list[[j]]@exprs)
  }

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}
# Transformation: flowVS transformation ------------------------------------------------------------
flowVS_trans <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  optimStat <- function(fs1D, cfLow=-1, cfHigh=10, MAX_BT=10^9) {
    if(cfLow>=cfHigh) {
      print("Warning: cfLow>=cfHigh, using default values")
      cfLow=-1
      cfHigh=10
    }
    cf = cfLow:cfHigh
    ncf =  length(cf)
    cfopt = rep(0,ncf-1)
    btopt = rep(0,ncf-1)

    for(i in 1:(ncf-1)) {
      tol = (exp(cf[i+1]) - exp(cf[i]))/10
      opt = suppressWarnings(optimize(f = flowVS:::flowVS1D, interval = c(exp(cf[i]),exp(cf[i+1])), fs1D, tol=tol, plot=FALSE, MAX_BT=MAX_BT))
      btopt[i] = opt$objective
      cfopt[i] = opt$minimum
    }

    minIdx = which.min(btopt)
    #now perform a local search around cfopt[minIdx] and plot
    del = cfopt[minIdx]/10
    btLocal = rep(0,11)
    btLocal[6] = btopt[minIdx]
    cfLocal = c(5:1,cfopt[minIdx],1:5)
    cfLocal[1:5] = cfopt[minIdx] - 5:1 * del
    cfLocal[7:11] = cfopt[minIdx] + 1:5 * del
    for(i in c(1:5,7:11)) {
      btLocal[i] = flowVS:::flowVS1D(cfLocal[i], fs1D)
    }

    minIdx = which.min(btLocal)
    return (cfLocal[minIdx])
  }

  estParamFlowVS <- function (fs, channels, cfLow = -1, cfHigh = 10, MAX_BT = 10^9) {
    checkmate::checkClass(fs, "flowSet")
    checkmate::checkClass(channels, "character")
    nmatch = which(channels %in% colnames(fs@frames[[names(fs@frames)[1]]]@exprs))
    if (length(nmatch) != length(channels)) {
      stop("At least one channel name is not present in the flowSet.")
    }
    cofactors = NULL
    for (col in channels) {
      fs1D = fs[, col]
      cf = optimStat(fs1D, cfLow = cfLow, cfHigh = cfHigh, MAX_BT = MAX_BT)
      cofactors = c(cofactors, cf)
    }
    return(cofactors)
  }

  data <- as(frame_list1, "flowSet")
  cofactors <- estParamFlowVS(fs = data, channels = col_names1)
  data_res <- flowVS::transFlowVS(data, channels = col_names1, cofactors)
  res <- lapply(1:length(data_res), function(i, frame_list, data_res) {
    rownames(data_res[[i]]@exprs) <- rownames(frame_list[[i]]@exprs)
    frame_list[[i]]@exprs <- data_res[[i]]@exprs
    return(frame_list[[i]])
  }, frame_list = frame_list, data_res = data_res)

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Transformation: Hyperlog transformation ------------------------------------------------------------
hyperlog_trans <- function(frame_list, col_names, a = 1, b = 1) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  hyperlogTrans <- hyperlog(col_names1, a = a, b = b)
  dataTransform <- lapply(frame_list1, function(x, hyperlogTrans) {
    data <- x@exprs
    x@exprs[, col_names1] <- matrix(as.vector(eval(hyperlogTrans)(data)), nrow = nrow(x@exprs), dimnames = list(NULL, col_names))
    return(x)
  }, hyperlogTrans = hyperlogTrans)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: linear transformation ------------------------------------------------------------
linear_trans <- function(frame_list, col_names, a = 2, b = 0) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  linearTrans <- flowCore::linearTransform(a = a, b = b)
  translist <- flowCore::transformList(col_names1, linearTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: Ln Transformation ------------------------------------------------------------
ln_trans <- function(frame_list, col_names, r = 1, d = 1) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  lnTrans <- flowCore::lnTransform(r = r, d = d)
  translist <- flowCore::transformList(col_names1, lnTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: log transformation ------------------------------------------------------------
log_trans <- function(frame_list, col_names, logbase = 10, r = 1, d = 1) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  logTrans <- flowCore::logTransform(logbase = logbase, r = r, d = d)
  translist <- flowCore::transformList(col_names1, logTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: Logicle transformation ------------------------------------------------------------
logicle_trans <- function(frame_list, col_names, w = 0.5, t = 262144, m = 4.5, a = 0) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  logicleTrans <- flowCore::logicleTransform(w = w, t = t, m = m, a = a) # 262144 for a 18 bit data range
  translist <- flowCore::transformList(col_names1, logicleTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: quadratic transformation ------------------------------------------------------------
quadratic_trans <- function(frame_list, col_names, a = 1, b = 1, c = 0) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  quadraticTrans <- flowCore::quadraticTransform(a = a, b = b, c = c)
  translist <- flowCore::transformList(col_names1, quadraticTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: Split Scale Transformation ------------------------------------------------------------
splitScale_trans <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  maxValue <- unlist(lapply(frame_list1, function(x, col_names1){max(x@exprs[,col_names1])}, col_names1 = col_names1))

  dataTransform <- lapply(1:length(frame_list1), function(i) {
    ssTransform  <- flowCore::splitScaleTransform("mySplitTransform",maxValue = maxValue[i])
    translist <- flowCore::transform(frame_list1[[i]], transformList(col_names1, ssTransform))
    return(translist)
  })
  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: truncate transformation ------------------------------------------------------------
truncate_trans <- function(frame_list, col_names, a = 1) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  truncateTrans <- flowCore::truncateTransform(a = a)
  translist <- flowCore::transformList(col_names1, truncateTrans)
  dataTransform <- lapply(frame_list1, flowCore::transform, translist = translist)

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}

# Transformation: Adaptive Box-Cox Transformation----------------------------------------------------
AdaptiveBoxCox_trans <- function(frame_list, col_names){
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  frame_list2 <- lapply(seq_along(frame_list1), function(i) {
    x <- frame_list1[[i]]
    x@exprs <- t(x@exprs)
    x@exprs <- cbind(Name = rownames(x@exprs), x@exprs)
    x@exprs <- rbind(Group = c("Group", rep(names(frame_list1)[i], ncol(x@exprs) - 1)), x@exprs)
    return(x)
  })

  dataTransform <- lapply(frame_list2, function(x) {
    exprs <- x@exprs[c("Group",col_names1),]
    temp <- ABCstats::ABCtransform(data.frame(exprs))
    x@exprs[c("Group",col_names1),] <- as.matrix(temp[,!colnames(temp) %in% "lambda"])
    return(x)
  })

  dataTransform <- lapply(seq_along(dataTransform), function(l) {
    x <- dataTransform[[l]]
    x@exprs <- t(x@exprs)
    x@exprs <- x@exprs[!rownames(x@exprs) %in% "Name",!colnames(x@exprs) %in% "Group"]
    x@exprs <- apply(x@exprs, 2, as.numeric)
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    rownames(x@exprs) <- rownames(frame_list[[l]]@exprs)
    return(x)
  })

  return(dataTransform)
}#new

# Transformation: Centered Log Ratio Transformation------------------------------------------------------------
CLR_trans <- function (frame_list, col_names){
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  clr <- function(x){
    return(log1p(x = x/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
  }

  dataTransform <- lapply(frame_list1, function(x) {
    exprs <- x@exprs[, col_names1]
    temp <- apply(exprs, 2, clr)
    x@exprs[, col_names1] <- temp
    return(x)
  })

  dataTransform <- lapply(dataTransform, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(dataTransform)
}#new

# Normalization ------------------------------------------------------------

# Normalization: Bead-based Normalization ------------------------------------------------------------
CATALYST_norm <- function(frame_list, col_names, beads_mass) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  ncell_frame <- sapply(frame_list1, function(x) nrow(x@exprs))

  frame_list1 <- CATALYST::concatFCS(frame_list1, by_time = F)

  frame_list1 <- CATALYST::normCytof(x = frame_list1, y = beads_mass, remove_beads = F, plot = F)

  res <- lapply(1:length(ncell_frame), function(i) {
    a <- sum(ncell_frame[0:(i-1)])
    b <- a + ncell_frame[i]
    frame_list[[i]]@exprs[, col_names] <- frame_list1@exprs[(a+1):b, col_names1]
    return(frame_list[[i]])
  })

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Normalization: gaussNorm ------------------------------------------------------------
gaussNorm_norm <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  data <- as(frame_list1, "flowSet")
  data_res <- flowStats::gaussNorm(flowset = data, channel.names = col_names1, max.lms = 1)$flowset
  res <- lapply(1:length(data_res), function(i) {
    frame_list[[i]]@exprs <- data_res[[i]]@exprs
    return(frame_list[[i]])
  })

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Normalization: warpSet ------------------------------------------------------------
warpSet_norm <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  data <- as(frame_list1, "flowSet")
  data_res <- flowStats::warpSet(data, col_names1)
  res <- lapply(1:length(data_res), function(i) {
    frame_list[[i]]@exprs <- data_res[[i]]@exprs
    return(frame_list[[i]])
  })

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Normalization: Z-Score ------------------------------------------------------------
ZScore_norm <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  data_res <- lapply(frame_list1, function(x) {
    exprs <- x@exprs[, col_names1]
    exprs[exprs == 0] <- NA
    temp <- scale(exprs)
    temp[is.na(exprs)] <- 0
    x@exprs[, col_names1] <- temp

    return(x)
  })


  res <- lapply(1:length(data_res), function(i) {
    frame_list[[i]]@exprs <- data_res[[i]]@exprs
    return(frame_list[[i]])
  })

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Normalization: Mean Normalization --------------------------------------------------
Mean_norm <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  data_res <- lapply(frame_list1, function(x) {
    exprs <- x@exprs[, col_names1]
    temp <-  sweep(exprs,2,apply(exprs,2,mean,na.rm=T),FUN="/")
    x@exprs[, col_names1] <- temp

    return(x)
  })

  res <- lapply(1:length(data_res), function(i) {
    frame_list[[i]]@exprs <- data_res[[i]]@exprs
    return(frame_list[[i]])
  })

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Normalization: Min-max Normalization --------------------------------------------------
MinMax_norm <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })
  mm_norm <- function(x) {(x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))}

  data_res <- lapply(frame_list1, function(x) {
    exprs <- x@exprs[, col_names1]
    temp <-  apply(exprs,2,mm_norm)
    x@exprs[, col_names1] <- temp
    return(x)
  })

  res <- lapply(1:length(data_res), function(i) {
    frame_list[[i]]@exprs <- data_res[[i]]@exprs
    return(frame_list[[i]])
  })

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}
# Signal Clean ------------------------------------------------------------

# Signal Clean: flowAI ------------------------------------------------------------
flowAI_signalC <- function(frame_list, col_names) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  flow_auto_qc <- function (fcsfiles, remove_from = "all", output = 1, timeCh = NULL,
                            timestep = NULL, second_fractionFR = 0.1, deviationFR = "MAD",
                            alphaFR = 0.01, decompFR = TRUE, ChExcludeFS = c("FSC",
                                                                             "SSC"), outlier_binsFS = FALSE, pen_valueFS = 500, max_cptFS = 3,
                            ChExcludeFM = c("FSC", "SSC"), sideFM = "both", neg_valuesFM = 1,
                            html_report = "_QC", mini_report = "QCmini", fcs_QC = "_QC",
                            fcs_highQ = FALSE, fcs_lowQ = FALSE, folder_results = "resultsQC") {
    if (is.character(fcsfiles)) {
      FileType <- toupper(strsplit(basename(fcsfiles[1]),
                                   split = "\\.")[[1]][-1])
      if (length(FileType) == 0) {
        warning("It was not possible to retrieve the file extension. The data will be processed as FCS.",
                call. = FALSE)
        FileType <- "FCS"
      } else if (FileType == "LMD") {
        set <- read.flowSet(files = fcsfiles, dataset = 2,
                            truncate_max_range = FALSE)
      } else {
        set <- read.flowSet(files = fcsfiles, truncate_max_range = FALSE)
      }
      names <- fcsfiles
    } else if (is(fcsfiles, "flowSet")) {
      FileType <- "FCS"
      set <- fcsfiles
      names <- flowCore::sampleNames(fcsfiles)
    } else if (is(fcsfiles, "flowFrame")) {
      FileType <- "FCS"
      set <- as(fcsfiles, "flowSet")
      names <- flowCore::identifier(fcsfiles)
      flowCore::sampleNames(set) <- names
    } else {
      stop("As first argument, use a flowSet or a character vector with the path of the FCS files")
    }
    N_cell_set <- flowAI:::flow_set_qc(set)
    area.color <- rep("red", length(set))
    if (missing(timeCh) || is.null(timeCh)) {
      timeCh <- flowAI:::findTimeChannel(set[[1]])
    }
    if (is.null(timeCh)) {
      warning("Impossible to retrieve the time channel automatically. The quality control can only be performed on signal acquisition and dynamic range.",
              call. = FALSE)
    }
    if (missing(timestep) || is.null(timestep)) {
      word <- which(grepl("TIMESTEP", names(keyword(set[[1]])),
                          ignore.case = TRUE))
      timestep <- as.numeric(keyword(set[[1]])[[word[1]]])
      if (!length(timestep)) {
        if (FileType == "LMD") {
          timestep <- 0.0009765625
        } else {
          warning("The TIMESTEP keyword was not found and hence it was set to 0.01. Graphs labels indicating time might not be correct",
                  call. = FALSE)
          timestep <- 0.01
        }
      }
    }
    if (second_fractionFR == "timestep") {
      second_fractionFR <- timestep
    } else if (second_fractionFR < timestep) {
      stop("The argument second_fractionFR must be greater or equal to timestep.",
           call. = FALSE)
    }
    if (folder_results != FALSE) {
      folder_results <- flowAI:::strip.sep(folder_results)
      dir.create(folder_results, showWarnings = FALSE)
      folder_results <- paste0(folder_results, .Platform$file.sep)
    } else {
      folder_results <- ""
    }
    out <- list()
    for (i in 1:length(set)) {
      filename_ext <- flowCore::identifier(set[[i]])
      filename <- sub("^([^.]*).*", "\\1", filename_ext)
      if (html_report != FALSE) {
        reportfile <- paste0(filename, html_report, ".html")
      }
      if (mini_report != FALSE) {
        minireport <- paste0(folder_results, mini_report,
                             ".txt")
        if (!file.exists(minireport)) {
          write.table(t(c("Name file", "n. of events",
                          "% anomalies", "analysis from", "% anomalies flow Rate",
                          "% anomalies Signal", "% anomalies Margins")),
                      minireport, sep = "\t", row.names = FALSE,
                      quote = FALSE, col.names = FALSE)
        }
      }
      if (fcs_QC != FALSE) {
        QC.fcs.file <- paste0(folder_results, filename, fcs_QC, ".fcs")
      }
      if (fcs_highQ != FALSE) {
        good.fcs.file <- paste0(folder_results, filename, fcs_highQ, ".fcs")
      }
      if (fcs_lowQ != FALSE) {
        bad.fcs.file <- paste0(folder_results, filename, fcs_lowQ, ".fcs")
      }
      cat(paste0("Quality control for the file: ", filename, "\n"))
      area <- area.color
      area[i] <- "blue"
      if (!is.null(timeCh)) {
        if (length(unique(exprs(set[[i]])[, timeCh])) == 1) {
          cat("The time channel contains a single value. It cannot be used to recreate the flow rate. \n")
          warning(paste0("The time channel in ", filename_ext,
                         " contains a single value. It cannot be used to recreate the flow rate. \n"),
                  call. = FALSE)
          TimeChCheck <- "single_value"
        } else {
          TimeChCheck <- NULL
        }
      } else {
        TimeChCheck <- "NoTime"
      }
      FSbinSize <- min(max(1, ceiling(nrow(set[[1]])/100)), 500)
      if (is.null(TimeChCheck)) {
        ordFCS <- flowAI:::ord_fcs_time(set[[i]], timeCh)
        rownames(ordFCS@exprs) <- rownames(exprs(set[[i]]))[order(exprs(set[[i]])[, timeCh])]
      } else {
        ordFCS <- set[[i]]
      }
      origin_cellIDs <- 1:nrow(ordFCS)
      FR_bin_arg <- list(second_fraction = second_fractionFR,
                         timeCh = timeCh, timestep = timestep)
      FR_QC_arg <- list(alpha = alphaFR, use_decomp = decompFR,
                        deviation = deviationFR)
      FS_bin_arg <- list(binSize = FSbinSize, timeCh = timeCh,
                         timestep = timestep, TimeChCheck = TimeChCheck)
      FS_QC_arg <- list(ChannelExclude = ChExcludeFS, pen_valueFS,
                        max_cptFS, outlier_binsFS)
      FM_QC_arg <- list(ChannelExclude = ChExcludeFM, side = sideFM,
                        neg_values = neg_valuesFM)
      if (is.null(TimeChCheck)) {
        FlowRateData <- try(do.call(flowAI:::flow_rate_bin, c(ordFCS, FR_bin_arg)))
        FlowRateQC <- try(do.call(flowAI:::flow_rate_check, c(ordFCS, list(FlowRateData), FR_QC_arg)))
      } else {
        FlowRateQC <- list()
        FlowRateQC$goodCellIDs <- origin_cellIDs
        FlowRateQC$res_fr_QC$badPerc <- 0
      }
      FlowSignalData <- try(do.call(flowAI:::flow_signal_bin, c(ordFCS, FS_bin_arg)))
      FlowSignalQC <- try(do.call(flowAI:::flow_signal_check, c(ordFCS, list(FlowSignalData), FS_QC_arg)))
      FlowMarginQC <- try(do.call(flowAI:::flow_margin_check, c(ordFCS, FM_QC_arg)))
      if (remove_from == "all") {
        res_list <- list(FlowRateQC, FlowSignalQC, FlowMarginQC)
        # shc
        filtered_list <- Filter(function(x) !inherits(x, "try-error"), res_list)
        if (length(filtered_list) == 0) {
          goodCellIDs <- c()
        } else if (length(filtered_list) == 1) {
          goodCellIDs <- filtered_list[[1]]$goodCellIDs
        } else {
          good_cell_ids_list <- lapply(filtered_list, function(x) x$goodCellIDs)
          goodCellIDs <- Reduce(intersect, good_cell_ids_list)
        }

        # shc
        # goodCellIDs <- intersect(FlowRateQC$goodCellIDs,
        #                          intersect(FlowSignalQC$goodCellIDs, FlowMarginQC$goodCellIDs))

        analysis <- "Flow Rate, Flow Signal and Flow Margin"
      } else if (remove_from == "FR_FS") {
        goodCellIDs <- intersect(FlowRateQC$goodCellIDs,
                                 FlowSignalQC$goodCellIDs)
        analysis <- "Flow Rate and Flow Signal"
      } else if (remove_from == "FR_FM") {
        goodCellIDs <- intersect(FlowRateQC$goodCellIDs,
                                 FlowMarginQC$goodCellIDs)
        analysis <- "Flow Rate and Flow Margin"
      } else if (remove_from == "FS_FM") {
        goodCellIDs <- intersect(FlowSignalQC$goodCellIDs,
                                 FlowMarginQC$goodCellIDs)
        analysis <- "Flow Signal and Flow Margin"
      } else if (remove_from == "FR") {
        goodCellIDs <- FlowRateQC$goodCellIDs
        analysis <- "Flow Rate"
      } else if (remove_from == "FS") {
        goodCellIDs <- FlowSignalQC$goodCellIDs
        analysis <- "Flow Signal"
      } else if (remove_from == "FM") {
        goodCellIDs <- FlowMarginQC$goodCellIDs
        analysis <- "Flow Margin"
      }
      badCellIDs <- setdiff(origin_cellIDs, goodCellIDs)
      totalBadPerc <- round(length(badCellIDs)/length(origin_cellIDs),
                            4)
      sub_exprs <- exprs(ordFCS)
      params <- parameters(ordFCS)
      keyval <- keyword(ordFCS)
      if (fcs_highQ != FALSE || output == 1) {
        goodfcs <- flowFrame(exprs = sub_exprs[goodCellIDs, ], parameters = params, description = keyval)
        rownames(goodfcs@exprs) <- rownames(sub_exprs[goodCellIDs, ])
        if (fcs_highQ != FALSE) {
          suppressWarnings(write.FCS(goodfcs, good.fcs.file))
        }
      }
      if (fcs_QC != FALSE || output == 2) {
        QCvector <- FlowSignalData$cellBinID[, "binID"]
        if (length(QCvector) > 9000)
          QCvector <- runif(length(QCvector), min = 1,
                            max = 9000)
        QCvector[badCellIDs] <- runif(length(badCellIDs),
                                      min = 10000, max = 20000)
        newFCS <- addQC(QCvector, remove_from, sub_exprs,
                        params, keyval)
        if (fcs_QC != FALSE) {
          suppressWarnings(write.FCS(newFCS, QC.fcs.file))
        }
      }
      if (length(badCellIDs) > 0 && fcs_lowQ != FALSE) {
        badfcs <- flowFrame(exprs = sub_exprs[badCellIDs, ], parameters = params, description = keyval)
        suppressWarnings(write.FCS(badfcs, bad.fcs.file))
      }
      if (mini_report != FALSE) {
        write.table(t(c(filename, as.integer(dim(set[[i]])[1]),
                        totalBadPerc * 100, analysis, FlowRateQC$res_fr_QC$badPerc *
                          100, FlowSignalQC$Perc_bad_cells$badPerc_tot *
                          100, FlowMarginQC$badPerc * 100)), minireport,
                    sep = "\t", append = TRUE, row.names = FALSE,
                    quote = FALSE, col.names = FALSE)
      }
      if (html_report != FALSE) {
        h_FS_graph <- round(0.4 * (ncol(ordFCS)), 1)
        if (!is.null(ChExcludeFS)) {
          ChannelExcludedFS <- as.character(grep(paste(ChExcludeFS,
                                                       collapse = "|"), ordFCS@parameters$name, value = TRUE))
        }
        if (!is.null(ChExcludeFM)) {
          ChannelExcludedFM <- as.character(grep(paste(ChExcludeFM,
                                                       collapse = "|"), ordFCS@parameters$name, value = TRUE))
        }
        template_path <- system.file("rmd", "autoQC_report.Rmd",
                                     package = "flowAI")
        new_template <- paste0(folder_results, filename,
                               "_template.Rmd")
        file.copy(template_path, new_template)
        if (folder_results != FALSE) {
          rmarkdown::render(new_template, html_document(),
                            output_dir = folder_results, output_file = reportfile,
                            quiet = TRUE)
        } else {
          rmarkdown::render(new_template, html_document(),
                            output_file = reportfile, quiet = TRUE)
        }
        file.remove(new_template)
      }
      if (output == 1) {
        out <- c(out, goodfcs)
      } else if (output == 2) {
        out <- c(out, newFCS)
      } else if (output == 3) {
        out[[i]] <- badCellIDs
        names(out)[i] <- filename
      }
    }
    if (output == 1 || output == 2) {
      if (length(out) == 1) {
        return(out[[1]])
      } else {
        OutSet <- as(out, "flowSet")
        flowCore::sampleNames(OutSet) <- names
        pData(OutSet) <- pData(set)
        return(OutSet)
      }
    }
    if (output == 3) {
      return(out)
    }
  }

  res <- lapply(frame_list1, flow_auto_qc,
                ChExcludeFS = setdiff(colnames(frame_list1[[1]]@exprs), col_names1),
                ChExcludeFM = setdiff(colnames(frame_list1[[1]]@exprs), col_names1),
                html_report = F, mini_report = F, fcs_QC = F,
                folder_results = FALSE
  )

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Signal Clean: flowCut ------------------------------------------------------------
flowCut_signalC <- function(frame_list, Segment = 200, col_names) {
  # col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  # col_names1 <- sub("\\(", "", col_names1)
  # col_names1 <- sub("\\)", "", col_names1)

  col_names1 <- sub("\\(.*", "", col_names)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  flowCut <- function (f, Segment = 500, Channels = NULL, Directory = NULL,
                       FileID = NULL, Plot = "Flagged Only", MaxContin = 0.1, MeanOfMeans = 0.13,
                       MaxOfMeans = 0.15, MaxValleyHgt = 0.1, MaxPercCut = 0.3,
                       LowDensityRemoval = 0.1, GateLineForce = NULL, UseOnlyWorstChannels = FALSE,
                       AmountMeanRangeKeep = 1, AmountMeanSDKeep = 2, PrintToConsole = FALSE,
                       AllowFlaggedRerun = FALSE, UseCairo = FALSE, UnifTimeCheck = 0.22,
                       RemoveMultiSD = 7, AlwaysClean = FALSE, IgnoreMonotonic = FALSE,
                       MonotonicFix = NULL, Measures = c(1:8), Verbose = FALSE) {
    start0 <- Sys.time()
    resTable <- matrix("", 17, 1)
    rownames(resTable) <- c("Is it monotonically increasing in time",
                            "Largest continuous jump", "Continuous - Pass", "Mean of % of range of means divided by range of data",
                            "Mean of % - Pass", "Max of % of range of means divided by range of data",
                            "Max of % - Pass", "Has a low density section been removed",
                            "% of low density removed", "How many segments have been removed",
                            "% of events removed from segments removed", "Worst channel",
                            "% of events removed", "FileID", "Type of Gating", "Was the file run twice",
                            "Has the file passed")
    resTable["Was the file run twice", ] <- "No"
    if (is.null(Directory)) {
      Directory <- paste0(getwd(), "/flowCut")
    }
    if (is.null(FileID)) {
      FileID <- Sys.time()
      FileID <- substr(FileID, start = 1, stop = 19)
      FileID <- gsub("-", "_", FileID)
      FileID <- gsub(":", "_", FileID)
      FileID <- gsub(" ", "__", FileID)
      if (Verbose == TRUE) {
        cat(paste0("The FileID is: ", FileID, "\n"))
      }
    }
    resTable["FileID", ] <- FileID
    if (!is(f, "flowFrame")) {
      message("f must be a flowFrame.")
      resTable["Has the file passed", ] <- "f must be a flowFrame"
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(Segment) || (Segment <= 0)) {
      message("Segment must be a number larger than 0.")
      resTable["Has the file passed", ] <- paste0("Segment",
                                                  "must be a number larger than 0.")
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (nrow(f) <= 3 * Segment) {
      message("Either your Segment size is too large or your number",
              " of cells is too small.")
      resTable["Has the file passed", ] <- paste0("Either your",
                                                  "Segment size is too large or your number of cells is too small.")
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.character(Directory)) {
      message("Directory must be a character.")
      resTable["Has the file passed", ] <- "Directory must be a character."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.character(FileID) && !is.numeric(FileID)) {
      message("FileID must be a character or a number.")
      resTable["Has the file passed", ] <- paste0("FileID",
                                                  "must be a character or a number.")
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!(Plot == "None" || Plot == "All" || Plot == "Flagged Only")) {
      message("Plot must be a character", "with one of the following",
              " options: 'None', 'All', or 'Flagged Only'.")
      resTable["Has the file passed", ] <- paste0("Plot must be a character",
                                                  " with one of the following options:", " 'None', 'All', or 'Flagged Only'.")
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(MaxContin) || (MaxContin < 0)) {
      message("MaxContin must be a number larger than 0.")
      resTable["Has the file passed", ] <- paste0("MaxContin",
                                                  " must be a number larger than 0.")
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(MeanOfMeans) || (MeanOfMeans < 0)) {
      message("MeanOfMeans must be a number larger than 0.")
      resTable["Has the file passed", ] <- paste0("MeanOfMeans",
                                                  " must be a number larger than 0.")
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(MaxOfMeans) || (MaxOfMeans < 0)) {
      message("MaxOfMeans must be a number larger than 0.")
      resTable["Has the file passed", ] <- "MaxOfMeans must be a number larger than 0."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(MaxValleyHgt) || (MaxValleyHgt < 0) || (MaxValleyHgt >
                                                            1)) {
      message("MaxValleyHgt must be a number between 0 and 1.")
      resTable["Has the file passed", ] <- "MaxValleyHgt must be a number between 0 and 1."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(MaxPercCut) || (MaxPercCut < 0) || (MaxPercCut >
                                                        1)) {
      message("MaxPercCut must be a number between 0 and 1.")
      resTable["Has the file passed", ] <- "MaxPercCut must be a number between 0 and 1."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(LowDensityRemoval) || (LowDensityRemoval <
                                           0) || (LowDensityRemoval > 1)) {
      message("LowDensityRemoval must be a number between 0 and 1.")
      resTable["Has the file passed", ] <- "LowDensityRemoval must be a number between 0 and 1."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.null(GateLineForce) && !is.numeric(GateLineForce)) {
      message("GateLineForce must be numeric.")
      resTable["Has the file passed", ] <- "GateLineForce must be numeric."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(UseOnlyWorstChannels)) {
      message("UseOnlyWorstChannels must be a logical (boolean).")
      resTable["Has the file passed", ] <- "UseOnlyWorstChannels must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (AmountMeanRangeKeep%%1 != 0) {
      message("AmountMeanRangeKeep must be an integer.")
      resTable["Has the file passed", ] <- "AmountMeanRangeKeep must be an integer."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (AmountMeanSDKeep%%1 != 0) {
      message("AmountMeanSDKeep must be an integer.")
      resTable["Has the file passed", ] <- "AmountMeanSDKeep must be an integer."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(PrintToConsole)) {
      message("PrintToConsole must be a logical (boolean).")
      resTable["Has the file passed", ] <- "PrintToConsole must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(AllowFlaggedRerun) && !is.character(Directory)) {
      message("AllowFlaggedRerun must be a logical (boolean).")
      resTable["Has the file passed", ] <- "AllowFlaggedRerun must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(UseCairo)) {
      message("UseCairo must be a logical (boolean).")
      resTable["Has the file passed", ] <- "UseCairo must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(UnifTimeCheck) || (UnifTimeCheck < 0) ||
        (UnifTimeCheck > 0.5)) {
      message("UnifTimeCheck must be numeric between 0 and 0.5.")
      resTable["Has the file passed", ] <- "UnifTimeCheck must be numeric between 0 and 0.5."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.null(RemoveMultiSD) && !is.numeric(RemoveMultiSD)) {
      message("RemoveMultiSD must be numeric.")
      resTable["Has the file passed", ] <- "RemoveMultiSD must be numeric."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(AlwaysClean)) {
      message("AlwaysClean must be a logical (boolean).")
      resTable["Has the file passed", ] <- "AlwaysClean must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(IgnoreMonotonic)) {
      message("IgnoreMonotonic must be a logical (boolean).")
      resTable["Has the file passed", ] <- "IgnoreMonotonic must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.null(MonotonicFix) && (!is.numeric(MonotonicFix) ||
                                   (MonotonicFix < 0))) {
      message("MonotonicFix must be NULL or a number greater than or equal to 0.")
      resTable["Has the file passed", ] <- "MonotonicFix must be NULL or a number greater than or equal to 0."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.numeric(Measures)) {
      message("Measures must be a numeric.")
      resTable["Has the file passed", ] <- "Measures must be numeric."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (max(Measures) > 8 || min(Measures) < 1) {
      message("Measures must be between 1 and 8.")
      resTable["Has the file passed", ] <- "Measures must be between 1 and 8."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (any(Measures%%1 != 0)) {
      message("Measures must be integers.")
      resTable["Has the file passed", ] <- "Measures must be integers."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.logical(Verbose)) {
      message("Verbose must be a logical (boolean).")
      resTable["Has the file passed", ] <- "Verbose must be a logical (boolean)."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    t.name <- flowCore::parameters(f)$name
    t.desc <- flowCore::parameters(f)$desc
    FSC.loc <- sort(unique(c(grep("fsc", tolower(t.name)), grep("fs lin",
                                                                tolower(t.name)), grep("*FS", t.name))))
    names(FSC.loc) <- NULL
    SSC.loc <- sort(unique(c(grep("ssc", tolower(t.name)), grep("ss lin",
                                                                tolower(t.name)), grep("*SS", t.name))))
    names(SSC.loc) <- NULL
    Time.loc <- sort(unique(c(grep("time", tolower(t.name)),
                              grep("time", tolower(t.desc)), grep("hdr-t", tolower(t.name)))))
    all.Time.loc <- Time.loc
    if (length(all.Time.loc) >= 2) {
      message("This file has ", length(all.Time.loc), " time channels. flowCut has selected to use ",
              t.name[Time.loc], " - ", t.desc[Time.loc], ".")
      Time.loc <- Time.loc[1]
      nonlap.loc <- all.Time.loc[which(all.Time.loc != Time.loc)]
      flowCore::parameters(f)$name[nonlap.loc] <- paste0(flowCore::parameters(f)$name[nonlap.loc],
                                                         "-Removed")
      colnames(exprs(f))[nonlap.loc] <- paste0(colnames(exprs(f))[nonlap.loc],
                                               "-Removed")
    }
    Extra.loc <- c(grep("pulse|width|length|count|sort classifier|event|phenograph|barcode",
                        tolower(t.name)))
    names(Extra.loc) <- NULL
    Extra.loc <- unique(Extra.loc)
    if (length(Extra.loc) >= 1 && Verbose == TRUE) {
      cat(paste0("Channels ", paste0(Extra.loc, collapse = ", "),
                 " are removed as they are not channels that need to be analyzed.\n"))
    }
    NoVariation <- NULL
    for (NoVar in seq_len(length(colnames(f)))) {
      if (sd(exprs(f)[, NoVar], na.rm = TRUE) == 0) {
        NoVariation <- c(NoVariation, NoVar)
      }
    }
    names(NoVariation) <- NULL
    if (length(NoVariation) >= 1 && Verbose == TRUE) {
      message("Channels ", paste0(NoVariation, collapse = ", "),
              " have no variation and have been removed from the analysis.")
    }
    if (!is.null(MonotonicFix)) {
      time.data <- exprs(f)[, Time.loc]
      if (all(time.data == cummax(time.data)) == FALSE) {
        message("Fixing file ", FileID, " for the monotonic time issue.")
        time.diff <- time.data[2:length(time.data)] - time.data[seq_len(length(time.data) -
                                                                          1)]
        diff.ind.strong <- which(time.diff < -MonotonicFix)
        diff.ind.all <- which(time.diff < 0)
        if (length(diff.ind.strong) > 0) {
          time.diffs <- abs(time.diff[diff.ind.strong])
          for (mono.ind in diff.ind.strong) {
            idx <- seq(mono.ind + 1, length(time.data))
            time.data[idx] <- time.data[idx] + time.diffs[which(mono.ind ==
                                                                  diff.ind.strong)]
          }
        }
        else {
          message("All of file ", FileID, "'s monotonic issues are larger than MonotonicFix and were not adjusted.")
        }
        if (Plot == "All" || (Plot == "Flagged Only")) {
          suppressWarnings(dir.create(paste0(Directory,
                                             "/Mono/"), recursive = TRUE))
          pngMono <- paste0(Directory, "/Mono/", gsub(".fcs",
                                                      "", FileID), "_Fix.png")
          if (PrintToConsole == FALSE) {
            if (UseCairo == TRUE) {
              CairoPNG(pngMono, width = 800, height = 600)
            }
            else {
              png(pngMono, width = 800, height = 600)
            }
          }
          par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5,
                                                            5, 3, 1))
          plot(seq_len(length(time.data)), time.data,
               pch = 19, cex = 0.2, lty = 2, col = "blue",
               main = paste0("Monotonic Time Correction: ",
                             gsub(".fcs", "", FileID)), xlab = "Cell Index",
               ylab = "Time")
          points(seq_len(length(exprs(f)[, Time.loc])),
                 exprs(f)[, Time.loc], pch = 19, cex = 0.2)
          if (length(diff.ind.all) > 0) {
            abline(v = diff.ind.all, col = "red")
          }
          if (length(diff.ind.strong) > 0) {
            abline(v = diff.ind.strong, col = "green4")
          }
          legend("bottomright", legend = c("Original",
                                           "Corrected", "Jump Fixed", "Jump Not Fixed"),
                 col = c("black", "blue", "green4", "red"),
                 lty = 1, cex = 1, lwd = c(2, 2, 1, 1))
          if (PrintToConsole == FALSE) {
            dev.off()
          } else {
            par(mfrow = c(1, 1), mar = c(5, 5, 4, 2),
                mgp = c(3, 1, 0))
          }
        }
        exprs(f)[, Time.loc] <- time.data
      }
    }
    MonotonicWithTime <- NULL
    for (MonoChan in seq_len(length(colnames(f)))) {
      if (all(exprs(f)[, MonoChan] == cummax(exprs(f)[, MonoChan])) ==
          TRUE) {
        MonotonicWithTime <- c(MonotonicWithTime, MonoChan)
      }
    }
    names(MonotonicWithTime) <- NULL
    MonotonicWithTime <- sort(unique(MonotonicWithTime))
    test <- match(all.Time.loc, MonotonicWithTime, nomatch = 0)
    if (any(test >= 1)) {
      MonotonicWithTime <- MonotonicWithTime[-test]
    }
    if (length(MonotonicWithTime) >= 1 && Verbose == TRUE) {
      message("Channels ", paste0(MonotonicWithTime, collapse = ", "),
              " are monotonically increasing in time and have been removed from the analysis.")
    }
    if (length(which(NoVariation == Time.loc)) >= 1) {
      message("Your time channel has no variation.")
      flowCore::parameters(f)$name[Time.loc] <- paste0(flowCore::parameters(f)$name[Time.loc],
                                                       "-Removed")
      colnames(exprs(f))[Time.loc] <- paste0(colnames(exprs(f))[Time.loc],
                                             "-Removed")
      Time.loc <- NULL
      if (length(which(is.na(match(all.Time.loc, NoVariation)))) >=
          1) {
        message("The first time channel will be replaced by the second time channel.")
        Time.loc <- all.Time.loc[-which(all.Time.loc ==
                                          NoVariation)][1]
        flowCore::parameters(f)$name[Time.loc] <- "Time"
        nonlap.loc <- all.Time.loc[which(all.Time.loc !=
                                           Time.loc)]
        flowCore::parameters(f)$name[nonlap.loc] <- paste0(flowCore::parameters(f)$name[nonlap.loc],
                                                           "-Removed")
        colnames(exprs(f))[nonlap.loc] <- paste0(colnames(exprs(f))[nonlap.loc],
                                                 "-Removed")
      }
    }
    if (length(Time.loc) == 0) {
      message("Your data does not have a time channel. flowCut will",
              " create one, but now flowCut will not be as fully",
              " functional as it could be. Consider recording the time",
              " for future projects.")
      exprs(f) <- cbind(exprs(f), seq_len(nrow(f)))
      colnames(exprs(f))[length(colnames(f)) + 1] <- "Time"
      pData(flowCore::parameters(f)) <- rbind(pData(flowCore::parameters(f)),
                                              c("Time", "Time", 262144, -111, 262143))
      rownames(flowCore::parameters(f))[length(colnames(f))] <- paste0("$P",
                                                                       length(colnames(f)))
      description(f)[paste0("P", length(colnames(f)), "DISPLAY")] <- "LIN"
      description(f)[paste0("flowCore_$P", length(colnames(f)),
                            "Rmax")] <- 262143
      description(f)[paste0("flowCore_$P", length(colnames(f)),
                            "Rmin")] <- 0
      description(f)[paste0("P", length(colnames(f)), "B")] <- "0"
      description(f)[paste0("P", length(colnames(f)), "R")] <- "262143"
      description(f)[paste0("P", length(colnames(f)), "N")] <- "Time"
      description(f)[paste0("P", length(colnames(f)), "G")] <- "1"
      description(f)[paste0("P", length(colnames(f)), "E")] <- "0,0"
      Time.loc <- length(colnames(f))
    }
    if (length(c(FSC.loc, SSC.loc)) == 0) {
      message("No FCS or SSC channels found.")
    }
    range_of_time <- max(exprs(f)[, Time.loc]) - min(exprs(f)[,
                                                              Time.loc])
    Time_test_passes <- TRUE
    if (range_of_time != 0) {
      uniformity_in_time_test <- abs(mean(exprs(f)[, Time.loc]) -
                                       (range_of_time/2 + min(exprs(f)[, Time.loc])))/range_of_time
      if (uniformity_in_time_test >= UnifTimeCheck) {
        message("The time channel does not appear to be distributed like an expected time channel would be.")
        Time_test_passes <- FALSE
      }
      uniformity_in_time_test_3 <- max(table(exprs(f)[, Time.loc]))
      if (uniformity_in_time_test_3 >= 0.05 * nrow(f)) {
        message("There appears to be an overwhelming amount of repeated time values.")
        Time_test_passes <- FALSE
      }
    } else {
      message("Range of the time channel is 0.")
      Time_test_passes <- FALSE
    }
    if (Time_test_passes == FALSE) {
      message("Time test(s) failed.")
      resTable["Has the file passed", ] <- "Time test(s) failed."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    CleanChan.loc <- (seq_len(ncol(f)))[-c(FSC.loc, SSC.loc,
                                           Time.loc, all.Time.loc, Extra.loc, NoVariation, MonotonicWithTime)]
    if (length(CleanChan.loc) == 0) {
      message("No marker channels to run flowCut on.")
      resTable["Has the file passed", ] <- "No marker channels to run flowCut on."
      return(list(frame = f, ind = NULL, data = resTable,
                  worstChan = NULL))
    }
    if (!is.null(Channels)) {
      if (all(is.character(Channels))) {
        Channels <- sort(unique(sapply(seq_len(length(Channels)),
                                       function(x) {
                                         # grep(tolower(Channels[x]), tolower(flowCore::parameters(f)$desc))
                                         which(tolower(flowCore::parameters(f)$desc) == tolower(Channels[x]))
                                       })))
      }
      CleanChan.loc <- intersect(CleanChan.loc, Channels)
    }
    ind.removed <- NA
    f.org <- f
    if (all(exprs(f)[, Time.loc] == cummax(exprs(f)[, Time.loc])) ==
        FALSE && !IgnoreMonotonic) {
      message("The flow frame is not monotonically increasing in time.")
      resTable["Is it monotonically increasing in time", ] <- "F"
    } else {
      resTable["Is it monotonically increasing in time", ] <- "T"
    }
    res.temp <- flowCut::removeLowDensSections(f, Time.loc = Time.loc,
                                               Segment, LowDensityRemoval, Verbose = Verbose)
    f <- res.temp$frame
    removeIndLowDens <- res.temp$rem.ind
    remove(res.temp)
    ifelse(length(removeIndLowDens) >= 1, resTable["Has a low density section been removed",
    ] <- "T", resTable["Has a low density section been removed",
    ] <- "F")
    resTable["% of low density removed", ] <- as.character(round(length(removeIndLowDens)/nrow(f.org),
                                                                 digits = 4) * 100)
    res.temp <- flowCut:::calcMeansAndSegmentsRemoved(f = f, Segment = Segment,
                                                      CleanChan.loc = CleanChan.loc, FirstOrSecond = "First",
                                                      MaxValleyHgt = MaxValleyHgt, MaxPercCut = MaxPercCut,
                                                      MaxContin = MaxContin, MeanOfMeans = MeanOfMeans, MaxOfMeans = MaxOfMeans,
                                                      GateLineForce = GateLineForce, UseOnlyWorstChannels = UseOnlyWorstChannels,
                                                      AmountMeanRangeKeep = AmountMeanRangeKeep, AmountMeanSDKeep = AmountMeanSDKeep,
                                                      RemoveMultiSD = RemoveMultiSD, AlwaysClean = AlwaysClean,
                                                      Verbose = Verbose, Measures = Measures, Time.loc = Time.loc)
    deletedSegments1 <- res.temp$deletedSegments
    quantiles1 <- res.temp$quantiles
    storeMeans1 <- res.temp$storeMeans
    meanRangePerc1 <- res.temp$meanRangePerc
    timeCentres1 <- res.temp$timeCentres
    typeOfGating <- res.temp$typeOfGating
    densityGateLine <- res.temp$densityGateLine
    cellDelete2 <- res.temp$cellDelete2
    choosenChans <- res.temp$choosenChans
    remove(res.temp)
    removed.ind <- NULL
    totalNumSeg <- floor(nrow(exprs(f))/Segment)
    if (length(deletedSegments1) == 0) {
      if (Verbose == TRUE) {
        cat("None deleted from flowCut segment removal.\n")
      }
      resTable["How many segments have been removed", ] <- as.character(0)
    } else {
      deletedSegments1 <- sort(unique(deletedSegments1), decreasing = FALSE)
      del.seg.list <- split(deletedSegments1, cumsum(c(1,
                                                       diff(deletedSegments1) != 1)))
      print.segs.rem <- sapply(seq_len(length(del.seg.list)),
                               function(x) {
                                 if (length(del.seg.list[[x]]) >= 2) {
                                   paste0(del.seg.list[[x]][1], "-", del.seg.list[[x]][length(del.seg.list[[x]])])
                                 } else {
                                   paste0(del.seg.list[[x]][1])
                                 }
                               })
      if (Verbose == TRUE) {
        cat(paste0("Removing segments ", paste0(print.segs.rem,
                                                collapse = ", "), " out of ", totalNumSeg, " segments.\n"))
      }
      resTable["How many segments have been removed", ] <- as.character(length(deletedSegments1))
      for (n in seq_len(length(deletedSegments1))) {
        if (deletedSegments1[n] == totalNumSeg)
          removed.ind <- c(removed.ind, (Segment * (deletedSegments1[n] -
                                                      1) + 1):nrow(exprs(f)))
        if (deletedSegments1[n] != totalNumSeg)
          removed.ind <- c(removed.ind, Segment * (deletedSegments1[n] -
                                                     1) + (seq_len(Segment)))
      }
      exprs(f) <- exprs(f)[-removed.ind, ]
    }
    resTable["% of events removed from segments removed", ] <- as.character(round(length(removed.ind)/nrow(f.org),
                                                                                  digits = 4) * 100)
    res.temp <- flowCut:::calcMeansAndSegmentsRemoved(f = f, Segment = Segment,
                                                      CleanChan.loc = CleanChan.loc, FirstOrSecond = "Second",
                                                      MaxValleyHgt = MaxValleyHgt, MaxPercCut = MaxPercCut,
                                                      MaxContin = MaxContin, MeanOfMeans = MeanOfMeans, MaxOfMeans = MaxOfMeans,
                                                      GateLineForce = GateLineForce, UseOnlyWorstChannels = UseOnlyWorstChannels,
                                                      AmountMeanRangeKeep = AmountMeanRangeKeep, AmountMeanSDKeep = AmountMeanSDKeep,
                                                      RemoveMultiSD = RemoveMultiSD, AlwaysClean = AlwaysClean,
                                                      Verbose = Verbose, Measures = Measures, Time.loc = Time.loc)
    quantiles2 <- res.temp$quantiles
    storeMeans2 <- res.temp$storeMeans
    meanRangePerc2 <- res.temp$meanRangePerc
    timeCentres2 <- res.temp$timeCentres
    remove(res.temp)
    maxDistJumped <- rep(0, max(CleanChan.loc))
    for (j in CleanChan.loc) {
      temp.vect <- rep(0, length(storeMeans2[[j]]) - 1)
      temp.vect <- abs(diff(storeMeans2[[j]]))/(quantiles1[[j]]["98%"] -
                                                  quantiles1[[j]]["2%"])
      maxDistJumped[j] <- max(temp.vect)
    }
    resTable["Largest continuous jump", ] <- as.character(round(max(maxDistJumped,
                                                                    na.rm = TRUE), digits = 3))
    if (resTable["Largest continuous jump", ] >= as.numeric(MaxContin)) {
      resTable["Continuous - Pass", ] <- "F"
      if (Verbose == TRUE) {
        message("The file has been flagged. The largest continuous jump ",
                "was larger than ", MaxContin * 100, "% of the range of the ",
                "2-98 percentile of the full data.")
      }
    } else {
      resTable["Continuous - Pass", ] <- "T"
    }
    resTable["Mean of % of range of means divided by range of data",
    ] <- as.character(round(mean(meanRangePerc2, na.rm = TRUE),
                            digits = 3))
    if (resTable["Mean of % of range of means divided by range of data",
    ] >= as.numeric(MeanOfMeans)) {
      if (Verbose == TRUE) {
        message("The file has been flagged. The means differ more than ",
                MeanOfMeans * 100, "% of the range of the 2-98 percentile of ",
                "the full data.")
      }
      resTable["Mean of % - Pass", ] <- "F"
    } else {
      resTable["Mean of % - Pass", ] <- "T"
    }
    resTable["Max of % of range of means divided by range of data",
    ] <- round(max(meanRangePerc2, na.rm = TRUE), digits = 3)
    worstChan <- min(which(meanRangePerc1 == max(meanRangePerc1,
                                                 na.rm = TRUE)))
    names.worschan <- flowCore::parameters(f)$name[worstChan]
    names(names.worschan) <- NULL
    resTable["Worst channel", ] <- names.worschan
    if (resTable["Max of % of range of means divided by range of data",
    ] >= MaxOfMeans) {
      if (Verbose == TRUE) {
        message("The file has been flagged. The max ranged means differ ",
                "more than ", MaxOfMeans * 100, "% of the range of the 2-98 ",
                "percentile of the full data.")
      }
      resTable["Max of % - Pass", ] <- "F"
    } else {
      resTable["Max of % - Pass", ] <- "T"
    }
    if (is.null(removed.ind) && is.null(removeIndLowDens))
      to.be.removed <- NULL
    if (is.null(removed.ind) && !is.null(removeIndLowDens))
      to.be.removed <- removeIndLowDens
    if (!is.null(removed.ind) && is.null(removeIndLowDens))
      to.be.removed <- removed.ind
    if (!is.null(removed.ind) && !is.null(removeIndLowDens)) {
      temp <- setdiff(seq_len(nrow(f.org)), removeIndLowDens)
      to.be.kept <- temp[setdiff(seq_len(length(temp)), removed.ind)]
      to.be.removed <- setdiff(seq_len(nrow(f.org)), to.be.kept)
    }
    resTable["% of events removed", ] <- as.character(round(length(to.be.removed)/nrow(f.org),
                                                            digits = 4) * 100)
    if ("TTTT" == paste0(resTable["Is it monotonically increasing in time",
    ], resTable["Continuous - Pass", ], resTable["Mean of % - Pass",
    ], resTable["Max of % - Pass", ])) {
      resTable["Has the file passed", ] <- "T"
    } else {
      resTable["Has the file passed", ] <- "F"
    }
    asterisks <- rep("", max(CleanChan.loc))
    if (UseOnlyWorstChannels == TRUE) {
      asterisks <- sapply(seq_len(max(CleanChan.loc)), function(x) {
        if (length(which(x == choosenChans)) >= 1) {
          asterisks <- " *"
        } else {
          asterisks <- ""
        }
        return(asterisks)
      })
    }
    if (Verbose == TRUE) {
      cat("Type of Gating: ", typeOfGating, ".\n", sep = "")
    }
    resTable["Type of Gating", ] <- typeOfGating
    if (resTable["Is it monotonically increasing in time", ] ==
        "T") {
      PassedMono <- "T"
    } else {
      PassedMono <- "F"
    }
    if (resTable["Continuous - Pass", ] == "T") {
      PassedCont <- "T"
    } else {
      PassedCont <- "F"
    }
    if (resTable["Mean of % - Pass", ] == "T") {
      PassedMean <- "T"
    } else {
      PassedMean <- "F"
    }
    if (resTable["Max of % - Pass", ] == "T") {
      PassedMax <- "T"
    } else {
      PassedMax <- "F"
    }
    if (resTable["Has the file passed", ] == "T") {
      FlaggedOrNot <- "Passed"
    } else {
      FlaggedOrNot <- "Flagged"
    }
    pngName <- paste0(Directory, "/", gsub(".fcs", "", FileID),
                      "_", FlaggedOrNot, "_", PassedMono, PassedCont, PassedMean,
                      PassedMax, ".png")
    if (Plot == "All" || (Plot == "Flagged Only" && FlaggedOrNot ==
                          "Flagged")) {
      z1 <- ceiling(sqrt(length(CleanChan.loc) + 2))
      if ((z1^2 - z1) >= (length(CleanChan.loc) + 2)) {
        z2 <- z1 - 1
      } else {
        z2 <- z1
      }
      suppressWarnings(dir.create(paste0(Directory), recursive = TRUE))
      if (AllowFlaggedRerun != TRUE && AllowFlaggedRerun !=
          FALSE && file.exists(AllowFlaggedRerun)) {
        pngName <- gsub(".png", "_2nd_run.png", pngName)
      }
      if (PrintToConsole == FALSE) {
        if (UseCairo == TRUE) {
          CairoPNG(filename = pngName, width = (z1) *
                     600, height = z2 * 600)
        } else {
          png(filename = pngName, width = (z1) * 600,
              height = z2 * 600)
        }
        par(mfrow = c(z2, z1), mar = c(7, 7, 4, 2), mgp = c(4,
                                                            1.5, 0), oma = c(0, 0, 5, 0))
        cex.size <- 3
      } else {
        par(mfrow = c(z2, z1), mar = c(5, 5, 4, 2), mgp = c(3,
                                                            1, 0), oma = c(0, 0, 5, 0))
        cex.size <- 1.5
      }
      for (x in CleanChan.loc) {
        plotDens(f.org, c(Time.loc, x), cex.main = cex.size,
                 cex.lab = cex.size, cex.axis = cex.size, main = paste0(round(meanRangePerc1[x],
                                                                              digits = 3), " / ", round(meanRangePerc2[x],
                                                                                                        digits = 3), " (", round(max(maxDistJumped[x]),
                                                                                                                                 digits = 3), ")", asterisks[x]))
        if (length(to.be.removed) != 0)
          points(exprs(f.org)[to.be.removed, c(Time.loc,
                                               x)], col = 1, pch = ".")
        if ((length(removeIndLowDens) != 0))
          points(exprs(f.org)[removeIndLowDens, c(Time.loc,
                                                  x)], pch = ".", cex = 1, col = "grey")
        lines(x = (timeCentres1), y = storeMeans1[[x]],
              cex.main = cex.size, cex.lab = cex.size, cex.axis = cex.size,
              lwd = 4, col = "deeppink2")
        lines(x = (timeCentres2), y = storeMeans2[[x]],
              cex.main = cex.size, cex.lab = cex.size, cex.axis = cex.size,
              lwd = 4, col = "brown")
        abline(h = c(quantiles1[[x]]["98%"], quantiles1[[x]]["2%"]),
               lwd = 4, col = "chocolate2")
        abline(h = c(quantiles2[[x]]["98%"], quantiles2[[x]]["2%"]),
               lwd = 4, col = "chocolate4")
      }
      x <- worstChan
      plotDens(f.org, c(Time.loc, x), cex.main = cex.size,
               cex.lab = cex.size, cex.axis = cex.size, main = paste0("Worst Channel without indices removed"))
      temp <- density(cellDelete2, adjust = 1)
      graphics::plot(temp, cex.main = cex.size, cex.lab = cex.size,
                     cex.axis = cex.size, main = "Density of summed measures")
      abline(v = densityGateLine, lwd = 2)
      title(main = paste0(FileID, " ", FlaggedOrNot, " ",
                          PassedMono, PassedCont, PassedMean, PassedMax),
            outer = TRUE, cex.main = cex.size + 1)
      if (PrintToConsole == FALSE) {
        dev.off()
      } else {
        par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), mgp = c(3,
                                                          1, 0))
      }
    }
    if (Verbose == TRUE) {
      if (resTable["Has the file passed", ] == "T") {
        cat("File Passed\n")
      } else {
        cat(paste0("The file has been flagged ", PassedMono,
                   PassedCont, PassedMean, PassedMax, "\n"))
      }
    }
    if (Verbose == TRUE) {
      cat("Cleaning completed in: ", TimePrint(start0), "\n",
          sep = "")
    }
    if (AllowFlaggedRerun == TRUE && resTable["Has the file passed",
    ] == "F") {
      if (Verbose == TRUE) {
        cat("Running flowCut a second time.\n")
      }
      res_flowCut <- flowCut(f = f, Segment = Segment, Channels = Channels,
                             Directory = Directory, FileID = FileID, Plot = Plot,
                             MaxContin = MaxContin, MeanOfMeans = MeanOfMeans,
                             MaxOfMeans = MaxOfMeans, MaxValleyHgt = MaxValleyHgt,
                             MaxPercCut = MaxPercCut, LowDensityRemoval = LowDensityRemoval,
                             GateLineForce = GateLineForce, UseOnlyWorstChannels = UseOnlyWorstChannels,
                             AmountMeanSDKeep = AmountMeanSDKeep, AmountMeanRangeKeep = AmountMeanRangeKeep,
                             PrintToConsole = PrintToConsole, AllowFlaggedRerun = pngName,
                             UseCairo = UseCairo, UnifTimeCheck = UnifTimeCheck,
                             RemoveMultiSD = RemoveMultiSD, AlwaysClean = AlwaysClean,
                             IgnoreMonotonic = IgnoreMonotonic, MonotonicFix = MonotonicFix,
                             Measures = Measures, Verbose = Verbose)
      if (res_flowCut$data["Has the file passed", ] == "Time test(s) failed.") {
        message("Time test(s) failed on the second run. Returning results from the first run for flowCut.")
        return(list(frame = f, ind = to.be.removed, data = resTable,
                    worstChan = worstChan))
      }
      indOfInd <- setdiff(seq_len(nrow(f.org)), to.be.removed)
      indOfInd <- sort(c(indOfInd[res_flowCut$ind], to.be.removed))
      resTableOfResTable <- res_flowCut$data
      resTableOfResTable["% of low density removed", ] <- as.character(round((nrow(f) *
                                                                                as.numeric(res_flowCut$data["% of low density removed",
                                                                                ]) + nrow(f.org) * as.numeric(resTable["% of low density removed",
                                                                                ]))/nrow(f.org), digits = 4))
      resTableOfResTable["How many segments have been removed",
      ] <- as.character(as.numeric(res_flowCut$data["How many segments have been removed",
      ]) + as.numeric(resTable["How many segments have been removed",
      ]))
      resTableOfResTable["% of events removed from segments removed",
      ] <- as.character(round((nrow(f) * as.numeric(res_flowCut$data["% of events removed from segments removed",
      ]) + nrow(f.org) * as.numeric(resTable["% of events removed from segments removed",
      ]))/nrow(f.org), digits = 4))
      resTableOfResTable["% of events removed", ] <- as.character(round((nrow(f) *
                                                                           as.numeric(res_flowCut$data["% of events removed",
                                                                           ]) + nrow(f.org) * as.numeric(resTable["% of events removed",
                                                                           ]))/nrow(f.org), digits = 4))
      resTableOfResTable["Was the file run twice", ] <- "Yes"
      return(list(frame = res_flowCut$frame, ind = indOfInd,
                  data = resTableOfResTable, worstChan = res_flowCut$worstChan))
    }
    return(list(frame = f, ind = to.be.removed, data = resTable,
                worstChan = worstChan))
  }

  res <- lapply(frame_list1, function(data, Segment = Segment, col_names1 = col_names1){
    flowCut(f = data, Segment = Segment, Channels = col_names1, Plot = "None")$frame
  }, Segment = Segment, col_names1 = col_names1)

  res <- lapply(res, function(x) {
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    return(x)
  })
  return(res)
}

# Signal Clean: flowClean ------------------------------------------------------------
flowClean_signalC <- function(frame_list, col_names, Segment = 200) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  res <- lapply(frame_list1, function(data, col_names1, Segment){
    flowClean::clean(data, vectMarkers = col_names1,
                     filePrefixWithDir = "aaaaaaaaaaaaaaaaaa", ext = "fcs", nCellCutoff = Segment)
  }, col_names1 = col_names1, Segment = Segment)

  res <- lapply(res, function(x) {
    x@exprs <- x@exprs[, -ncol(x@exprs)]
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    x@parameters@data <- x@parameters@data[-nrow(x@parameters@data),]
    return(x)
  })
  return(res)
}

# Signal Clean: PeacoQC ------------------------------------------------------------
PeacoQC_signalC <- function(frame_list, col_names, min_cells, max_bins, step, technique) {
  col_names1 <- stringr::str_extract(col_names, "\\(.*\\)")
  col_names1 <- sub("\\(", "", col_names1)
  col_names1 <- sub("\\)", "", col_names1)

  frame_list1 <- lapply(frame_list, function(x) {
    colnames(x@exprs) <- stringr::str_extract(colnames(x@exprs), "\\(.*\\)")
    colnames(x@exprs) <- sub("\\(", "", colnames(x@exprs))
    colnames(x@exprs) <- sub("\\)", "", colnames(x@exprs))
    return(x)
  })

  res <- lapply(frame_list1, function(data, col_names1, technique){
    if (technique == "FC"){
      data <- PeacoQC::RemoveMargins(ff = data, channels = col_names1, output="frame")
    }
    data <- PeacoQC::PeacoQC(ff = data, channels = col_names1, plot =F,
                             save_fcs =F, report = F, output_directory = NULL,
                             min_cells = min_cells, max_bins = max_bins, step = step)[["FinalFF"]]

    return(data)
  }, col_names1 = col_names1, technique = technique)

  res <- lapply(res, function(x) {
    x@exprs <- x@exprs[, -ncol(x@exprs)]
    colnames(x@exprs) <- colnames(frame_list[[1]]@exprs)
    x@parameters@data <- x@parameters@data[!x@parameters@data$name == "Original_ID",]
    return(x)
  })
  return(res)
}
