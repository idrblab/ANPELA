# read single fcs file
scaleData <- function(x, range = c(0, 4.5)) {
  (x - min(x))/(max(x) - min(x)) * (range[2] - range[1]) + range[1]
}
readfcs_single <- function(fcsFile) {
  if (class(fcsFile) == "list") {
    name <- names(fcsFile)
    fcs <- fcsFile[[1]]
  }
  if (class(fcsFile) == "character") {
    name <- sub(".fcs$", "", basename(fcsFile)) # abstract the file name
    fcs <- suppressWarnings(flowCore::read.FCS(fcsFile, transformation = FALSE, alter.names = TRUE)) 
  }
  
  pd <- fcs@parameters@data # abstract the protein name and its characterization
  
  exclude_channels <- grep("Time|Event", colnames(fcs@exprs), ignore.case = TRUE)
  marker_id <- setdiff(seq_along(colnames(fcs@exprs)), exclude_channels) # remove the "Time" and "Event"
  
  size_channels <- grep("FSC|SSC", colnames(fcs@exprs), ignore.case = TRUE)
  transMarker_id <- setdiff(marker_id, size_channels) # remove the "FSC" and "SSC" in transMarker_id (this variable will not be used)
  
  marker_id <- seq_along(colnames(fcs@exprs))
  
  data <- fcs@exprs
  exprs <- data[, marker_id, drop = FALSE] # data after removed the "Time" and "Event"
  
  if (length(size_channels) > 0) { # judge whether it is FC data
    if (any(size_channels %in% marker_id)) {
      used_size_channel <- size_channels[size_channels %in% marker_id]
      used_size_channel_id <- match(used_size_channel, marker_id)
      exprs[, used_size_channel_id] <- apply(exprs[, used_size_channel_id, drop = FALSE], 2,
                                             function(x) scaleData(x, range = c(0,4.5))) # scale the "FSC" and "SSC" data
    }
  }
  col_names <- paste0(pd$desc, "(", pd$name, ")")
  colnames(exprs) <- col_names[marker_id]
  # if (!nrow(exprs) == 0) {
  #   row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
  # }
  return(exprs)
}


readfcs_multi <- function(fcsFiles = files, mergeMethod = c("ceil", "all", "fixed", "min"), 
                          fixedNum = 10000, sampleSeed = 123){
  mergeMethod <- match.arg(mergeMethod)
  # exprsL <- mapply(readfcs_single, fcsFiles, SIMPLIFY = FALSE) # use "readfcs_single" function to read each file
  exprsL <- list()
  for (i in 1:length(fcsFiles)) {
    exprsL[[i]] <- readfcs_single(fcsFiles[i])
  }
  if (is.numeric(sampleSeed)) 
    set.seed(sampleSeed)
  eventCountTest <- suppressWarnings(any(lapply(exprsL, function(x) if (nrow(x) < fixedNum) { # comparation between "nrow(x)" and "fixedNum" will reflect the mergement
    1
  } else {
    0
  })))
  if (mergeMethod == "fixed" && eventCountTest == TRUE) {
    warning("One or more FCS files have less events than specified fixedNum; using lowest as fixedNum")
    fixedNum <- min(rapply(exprsL, nrow))
  }
  switch(mergeMethod, ceil = {
    mergeFunc <- function(x) {
      if (nrow(x) < fixedNum) {
        x
      } else {
        x[sample(nrow(x), size = fixedNum, replace = FALSE), , drop = FALSE]
      }
    }
    merged <- do.call(rbind, lapply(exprsL, mergeFunc))
  }, all = {
    merged <- do.call(rbind, exprsL)
  }, fixed = {
    mergeFunc <- function(x) {
      x[sample(nrow(x), size = fixedNum, replace = ifelse(nrow(x) < fixedNum, TRUE, FALSE)), , drop = FALSE]
    }
    merged <- do.call(rbind, lapply(exprsL, mergeFunc))
  }, min = {
    minSize <- min(sapply(exprsL, nrow))
    mergeFunc <- function(x) {
      x[sample(nrow(x), size = minSize, replace = FALSE), , drop = FALSE]
    }
    merged <- do.call(rbind, lapply(exprsL, mergeFunc))
  })
  return(merged)
}


