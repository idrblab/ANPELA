load.Data <- function(files.paths, index, measurement.time, measurement.condition = NULL, TIM ){
  if (class(files.paths) == "character") {
    data.raw <- flowCore::read.flowSet(files = files.paths)
  } else if (class(files.paths) == "list") {
    data.raw <- as(files.paths, "flowSet")
  }
  # data.raw <- read.flowSet(files = files.paths)
  data.info <- pData(flowCore::parameters(data.raw[[1]]))
  
  if (missing(index)) {
    index <- 1:length(data.info$desc)
  } else {
    index <- index
  }
  
  
  # make a flowSet with only the abundances of the markers of choice. (optional: arcsinh transform)
  data <- flowCore::fsApply(data.raw, function(x, cofactor=5){
    colnames(x) <- make.names(data.info$desc)
    # colnames(x) <- paste0(data.info$desc, "(", data.info$name, ")")
    index <- make.names(index)
    expr <- flowCore::exprs(x)
    expr <- expr[, index]
    flowCore::exprs(x) <- expr
    return(x)
  })
  
  # how many cells are in each dataset
  N <- flowCore::fsApply(data, function(x){ x <- dim(flowCore::exprs(x))[1]})
  
  if (TIM =="slingshot_FLOWMAP"){
    expr <- as(data, "list")
  } else {
    # make a matrix with the abundances and a vector with the time
    expr <- as.data.frame(flowCore::fsApply(data, flowCore::exprs))
  }
  timepoint <- unlist(sapply(c(1:length(N)), function(i){ rep(measurement.time[i],N[i])}, simplify = FALSE))
  
  if (is.null(measurement.condition) == FALSE) {
    condition <- unlist(sapply(c(1:length(N)), function(i){ rep(measurement.condition[i],N[i])}, simplify = FALSE))
  } else {
    condition <- NULL
  }
  
  # find the mean at each timepoint
  exp_means <- flowCore::fsApply(data,flowCore::each_col,mean)
  
  D <- list()
  D$N <- N
  D$expr <- expr
  D$exp_means <- exp_means
  D$timepoint <- timepoint
  D$ntimepoints <- length(N)
  D$info <- data.info
  if (!is.null(condition)) {D$condition <- condition}
  return (D)
}
