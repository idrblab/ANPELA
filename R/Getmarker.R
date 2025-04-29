
#' @title Get Marker
#' @description getmarker enables the preparation of marker indexes for data processing and performance assessment, with subsequent manual removal of non-protein columns as needed.
#'
#' @param datapath Character, the absolute path of the folder storing the FCS raw data files and metadata file.
#'
#' @return A character including marker indexes for data processing and performance evaluation, with subsequent manual removal of non-protein columns as needed
#' @export
#'
#' @examples
#' \donttest{
#' }
Getmarker <- function(datapath) {
  files <- dir(path = datapath, pattern = "fcs$", full.names = TRUE)
  a <- flowCore::read.FCS(filename = files[[1]], transformation = FALSE, alter.names = TRUE)
  colnames(a@exprs) <- paste0(a@parameters@data$desc, "(", a@parameters@data$name, ")")
  cat(paste0("\"", paste0(as.character(colnames(a@exprs)), collapse = ", \n"), "\""))
}
