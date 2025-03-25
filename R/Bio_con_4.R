### Calculate Biological Consistency
Bio_con <- function(D, Pathway_Hierarchy_file, nruns = 3, dr_method = TIres$dr_method, TIres = TIres){
  # Pathway_Hierarchy_file : the file with the pathway info. Put with path, extension
  Pathway_Hierarchy <- suppressMessages(as.data.frame(read.csv(file = Pathway_Hierarchy_file2, header = T, stringsAsFactors = F)))
  #Pathway_Hierarchy <- suppressMessages(as.data.frame(readr::read_delim(file = Pathway_Hierarchy_file, ";", escape_double = FALSE, trim_ws = TRUE)))
  #colnames(Pathway_Hierarchy) <- stringr::str_replace_all(colnames(Pathway_Hierarchy), "\\(.*", "")
  #Pathway_Hierarchy <- as.data.frame(apply(Pathway_Hierarchy, 2, function(x) stringr::str_replace_all(x, "\\(.*", "")))
  # 11.17 处理特殊字符
  # Pathway_Hierarchy <- stringr::str_replace_all(Pathway_Hierarchy, "[[:punct:]]", "_")
  npar <- ncol(D$expr)
  
  # decode the TI method used (cr: curve reconstruction, dr: dimensionality reduction, act: activation, transf: transformation)
  FUN <- TIres$cr_method
  #setup parallel backend to use many processors
  # cores <- detectCores()
  # cl <- makeCluster(cores[1]-1) #do not overload the computer
  # registerDoParallel(cl)
  
  finalMatrix <- foreach::foreach(j = 1:nruns, .export = FUN, .verbose = F) %do% {
    set.seed(j)
    tempMatrix <- do.call(FUN, list(D, dr_method))
    tempMatrix
  }
  
  #stop cluster
  # stopCluster(cl)
  
  temp_score <- rep(0,nruns) 
  
  for (k in 1:nruns){
    # sort the plot values in ascending pseudotime order
    ix <- order(finalMatrix[[k]]$pseudotime)
    dat <- D$expr[ix,]
    Pseudotime <- finalMatrix[[k]]$pseudotime[ix]
    if (finalMatrix[[k]]$lineages > 1){
      print(paste0("The algorithm ", TIres$cr_method, " with ", TIres$dr_method, " returned a branching trajectory. Cannot calculate the metric under Criterion Cd!"))
      return(paste0("The algorithm ", TIres$cr_method, " with ", TIres$dr_method, " returned a branching trajectory. Cannot calculate the metric under Criterion Cd!"))
    }
    
    # rescale pseudotime to [0,1] for the figures between algorithms to be comparable
    Pseudotime_zscore <- (Pseudotime - min(Pseudotime, na.rm = T))/(max(Pseudotime, na.rm = T) - min(Pseudotime, na.rm = T))
    
    # check also the reverse direction
    Pseudotime_zscore_rev <- reverse_pseudotime(Pseudotime_zscore)
    
    # 11.17
    # colnames(dat) <- stringr::str_replace_all(colnames(dat), "[[:punct:]]", "_")
    
    temp <- check_pairs(dat, Pseudotime_zscore, known_pairs = Pathway_Hierarchy)
    temp_rev <- check_pairs(dat, Pseudotime_zscore_rev, known_pairs = Pathway_Hierarchy)
    
    temp_score[k] <- max(temp, temp_rev)
  }
  return(max(temp_score)/sum(!is.na(Pathway_Hierarchy)))
}
