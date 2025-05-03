### Calculate Biological Consistency
Bio_con <- function(D, Pathway_Hierarchy_file = NULL, nruns = 3, dr_method = NULL, TIres = NULL){
  if (is.null(Pathway_Hierarchy_file)) {
    Lineage_info <- data.frame(TIres[["pseudotime"]])
    
    if ("FLOWMAP" %in% TIres$dr_method) {
      Lineage_info$condition <- as.character(TIres[["condition"]])
    } else {
      Lineage_info$condition<- D[["condition"]]
    }

    Lineage_condition <- as.table(matrix(nrow = ncol(Lineage_info)-1, ncol = length(unique(Lineage_info$condition))))
    dimnames(Lineage_condition) <- list(Lineage = colnames(Lineage_info)[1:(ncol(Lineage_info)-1)], 
                                        condition = unique(Lineage_info$condition))
    for (i in 1:(ncol(Lineage_info)-1)){
      condition_table<- table(Lineage_info[!is.na(Lineage_info[,i]),ncol(Lineage_info)])
      for(j in 1:(ncol(Lineage_condition))){
        if(grepl(colnames(Lineage_condition)[j],dimnames(condition_table))){
          Lineage_condition[i,j] <- condition_table[colnames(Lineage_condition)[j]]
        } else {
          Lineage_condition[i,j] <- 0
        }
      }
      
    }
    p.value.l <- chisq.test(Lineage_condition,rescale.p =T)[["p.value"]]
  
    if (p.value.l <= 1e-80){
      return(1)
    } else {
      return(0)
    }
  } else {
    # Pathway_Hierarchy_file : the file with the pathway info. Put with path, extension
    Pathway_Hierarchy <- suppressMessages(as.data.frame(read.csv(file = Pathway_Hierarchy_file, header = T, stringsAsFactors = F)))
    npar <- ncol(D$expr)
    FUN <- TIres$cr_method
    
    finalMatrix <- foreach::foreach(j = 1:nruns, .export = FUN, .verbose = F) %do% {
      set.seed(j)
      tempMatrix <- do.call(FUN, list(D, dr_method))
      tempMatrix
    }
    
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
      
      temp <- check_pairs(dat = dat, t = Pseudotime_zscore, known_pairs = Pathway_Hierarchy)
      temp_rev <- check_pairs(dat = dat, t = Pseudotime_zscore_rev, known_pairs = Pathway_Hierarchy)
      
      temp_score[k] <- max(temp, temp_rev)
    }
    return(max(temp_score)/sum(!is.na(Pathway_Hierarchy)))
  }
}
