#time_metric--------------------------------------------------------
time_metric <- function(TIres, D, nruns = 10) {
  max_Ps <- c()
  # load the data
  for (lineage in 1:TIres[["lineages"]]) {
    if ("FLOWMAP" %in% TIres$dr_method) {
      input_matrix <- cbind(TIres$expr, TIres$timepoint,  as.matrix(TIres$pseudotime)[, lineage, drop = FALSE])
      nruns <- min(nruns,TIres$N)
    } else {
      input_matrix <- cbind(D$expr, D$timepoint,  as.matrix(TIres$pseudotime)[, lineage, drop = FALSE])
      nruns <- min(nruns, min(D$N[D$N != 0]))
    }
    
    colnames(input_matrix)[(length(input_matrix)-1):(length(input_matrix))] <- c("D.timepoint", "pseudot")
    
    #skip NA pseudot
    input_matrix <- input_matrix[!is.na(input_matrix[,"pseudot"]),]
    
    npar <- which(colnames(input_matrix) == "D.timepoint") - 1
    exp.time <- input_matrix[, npar + 1]
    pseudo.time <- input_matrix[, npar + 2]
    
    
    P <- 0
  
    # create all possible combinations of experimental timepoints
    # measurement.times = unique(exp.time)
    measurement.times <- unique(sort(exp.time))
    exp.time.combns <- combn(length(measurement.times),2)
    
    # initialization of P( PT[i]<PT[j] | ET[i]<ET[j], for i != j)
    hits <- 0
    
    for (c in 1:dim(exp.time.combns)[2]){
      # Create all pairs from a sample of n cells where ET[i] < ET[j]
      set.seed(21)
      i <- sample(which(exp.time == measurement.times[exp.time.combns[1,c]]), nruns, replace = T)
      set.seed(22)
      j <- sample(which(exp.time == measurement.times[exp.time.combns[2,c]]), nruns, replace = T)
      pairs <- rbind(rep(i, each = nruns), rep(j, nruns))
      
      # Check for every pair whether also PT[i] < PT[j]
      for (k in 1:nruns^2){
        if (pseudo.time[pairs[1,k]] <= pseudo.time[pairs[2,k]]) 
          hits <- hits + 1
      }
    
    P <- hits / (dim(exp.time.combns)[2]*nruns^2)
  }
    max_Ps[lineage] <-  P
  }
  return(mean(max_Ps))
}