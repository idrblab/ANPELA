
### Calculate Robustness
# nruns: define how many different datasets will be contrasted
# cell.subset : how many cells should be used (we used 20% )

Robustness <- function(TIres, D, nruns = 4, cell.subset = 0.9, clustering.var = NULL,dataset_name =NULL, ...){
  
  
  Robustness_result <- c()
  
  if ("FLOWMAP" %in% TIres$dr_method) {
    input_matrix <- cbind(TIres$expr, TIres$timepoint, TIres$condition, TIres$pseudotime[drop = FALSE])
    npar <- grep("timepoint", colnames(input_matrix)) - 1
    if (ncol(as.matrix(TIres$pseudotime)) == 1) {
      colnames(input_matrix)[npar+3] <- c("pseudot1")
    } else {
      colnames(input_matrix) <- gsub("Lineage", "pseudot", colnames(input_matrix))
    }
    
    colnames(input_matrix)[(npar + 1):(npar + 2)] <- c("D.timepoint", "condition")
    
    expr_matrix <- list(expr = input_matrix[,1:(npar + 2)])
    grouped_data <- list()
    for (i in 1:nrow(expr_matrix[["expr"]])) {
      condition <-  as.character(expr_matrix[["expr"]][i, "condition"])
      timepoint <-  as.character(expr_matrix[["expr"]][i, "D.timepoint"])
      name <- paste0(condition,"-",timepoint)
      
      if (!((name) %in% names(grouped_data))) {
        grouped_data[[name]] <- data.frame()
      }
      
      grouped_data[[name]] <- rbind(grouped_data[[name]], expr_matrix[["expr"]][i,1:npar])
    }
    expr_matrix <- list(expr = grouped_data)
    perfile.cells <- nrow(expr_matrix[["expr"]][[1]])
    perfile.nsub.cells <- round(cell.subset * perfile.cells)
    
  } else {
    input_matrix <- cbind(D$expr, D$timepoint, TIres$pseudotime)
    npar <- grep("timepoint", colnames(input_matrix)) - 1
    if(ncol(as.matrix(TIres$pseudotime)) == 1){
      colnames(input_matrix)[npar+2] <- c("pseudot1")  
    }else{
      colnames(input_matrix) <- gsub("Lineage","pseudot",colnames(input_matrix))
    }
    colnames(input_matrix)[npar+1] <- c("D.timepoint")
    expr_matrix <- input_matrix[,1:npar]
    timepoint.annot <- input_matrix[, npar + 1] 
  }
  
  #subset_order will have the order of each cell, after the TI method is applied on each subset.
  ncells <- nrow(input_matrix)
  nsub.cells <- round(cell.subset*ncells)
  initial_order <- subset_order <- as.data.frame(matrix(nrow =  nsub.cells, ncol = nruns))  
  
  ####
  #### ~ ~ ~ apply alg on subsampled data ~ ~ ~
  ####
  
  FUN <- TIres$cr_method
  
  finalMatrix <- foreach::foreach(j = 1:nruns, .export = FUN, .verbose = F) %do% {
    # Subsampling from the initial dataset, with stratification
    set.seed(j)
    
    if ("FLOWMAP" %in% TIres$dr_method) {
      subIndex <- sample(1:perfile.cells, perfile.nsub.cells)
      subMatrix <- list(expr = lapply(expr_matrix[["expr"]], function(df) df[subIndex,]))
      tempMatrix <-  BBmisc::suppressAll(do.call(FUN, list(subMatrix, TIres$dr_method, clustering.var = clustering.var, subsamples = FALSE, 
                                      cluster.numbers = perfile.nsub.cells,dataset_name=dataset_name)))
      input_subIndex <- unlist(lapply(subMatrix[["expr"]], function(df) as.integer(rownames(df))),use.names = FALSE)
      tempMatrix[["subIndex"]] <- input_subIndex
    } else {
      subIndex <- sample(1:ncells, nsub.cells)
      subMatrix <- list(expr = expr_matrix[subIndex, ], 
                        ntimepoints = length(unique(input_matrix$D.timepoint)),
                        timepoint = timepoint.annot[subIndex])  
      tempMatrix <- do.call(FUN, list(subMatrix, TIres$dr_method))
      tempMatrix[["subIndex"]] <- subIndex
    }
    return(tempMatrix)
  }
  
  spear_final <- kend_final <- c()
  
  for (n in 1:nruns){
    spear <- kend <- c()
    spear_max <- kend_max <- c()
    for (lineage in 1:finalMatrix[[n]]$lineages){
      subset_order[,n] <- as.matrix(finalMatrix[[n]][["pseudotime"]])[, lineage]
      for (l in 1:ncol(as.matrix(TIres$pseudotime))){ 
        initial_order[,n] <- input_matrix[,paste0("pseudot",l)][finalMatrix[[n]]$subIndex]
        spear <- c(spear, cor(initial_order[,n],subset_order[,n], use="pairwise.complete.obs", method = "spearman"))
        kend <- c(kend, cor(initial_order[,n],subset_order[,n], use="pairwise.complete.obs", method = "kendall"))
      }
      spear_max <- c(spear_max,max(spear,na.rm =T))
      kend_max <- c(kend_max,max(kend,na.rm =T))
    }
    spear_benchmark <- mean(abs(spear_max))
    kend_benchmark <-  mean(abs(kend_max))
    
    ratio <- min(finalMatrix[[n]]$lineages, TIres[["lineages"]])/max(finalMatrix[[n]]$lineages, TIres[["lineages"]])
    spear_final <- c(spear_final, ratio * spear_benchmark)
    kend_final <- c(kend_final, ratio * kend_benchmark)
    
  }
  
  Robustness_result <- c(mean(spear_final), mean(kend_final))
  names(Robustness_result) <- c("spearman", "kendall")
  
  
  return(list(Robustness_result = Robustness_result, input_matrix = input_matrix, finalMatrix = finalMatrix))  
}