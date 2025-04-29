check_pairs <- function(dat, t, known_pairs){
  
  # add the spline fit
  spline.df = calc_spline(dat,t)
  
  max_matrix = data.frame(protein = names(dat),
                          value = apply(spline.df[,-1],2,max),
                          time =  spline.df[apply(spline.df[,-1],2,which.max),1])
  
  max_matrix <- max_matrix[order(max_matrix$time),]
  
  temp <- 0
  for (l in 1:ncol(known_pairs)){
    check_prot <- names(known_pairs)[l]
    index_examined <- which(max_matrix$protein == check_prot)
    # run through all proteins in the established Hierarchy
    for (m in 1:length(known_pairs[,l][known_pairs[,l] != ""])){
      subseq_prot <- known_pairs[m,l]
      if (is.na(subseq_prot) == FALSE){
        second_index <- which(as.character(max_matrix$protein) == as.character(subseq_prot))
        if (max_matrix$time[index_examined] < max_matrix$time[second_index]) {
          temp = temp + 1
        }
      }
    }
  }
  return(temp)
}
