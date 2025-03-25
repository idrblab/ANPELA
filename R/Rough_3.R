DeLorean.calc.roughness <- function (x) {
  N <- length(x)
  stopifnot(N > 0)
  if (0 == N) 
    return(0)
  S <- sd(x)
  if (S == 0) 
    return(0)
  sqrt(sum((x[1:(N - 1)] - x[2:N])^2)/(N - 1))/S
}

Rough <- function(TIres, D){
  TEST_list <- c()
  observedRoughnes_list <- list()
  
  # load the data
  for (lineage in 1:TIres[["lineages"]]) {
    if ("FLOWMAP" %in% TIres$dr_method) {
      input_matrix <- cbind(TIres$expr, TIres$timepoint, as.matrix(TIres$pseudotime)[, lineage, drop = FALSE])
    } else {
      input_matrix <- cbind(D$expr, D$timepoint, as.matrix(TIres$pseudotime)[, lineage, drop = FALSE])
    }
    
    colnames(input_matrix)[(length(input_matrix)-1):(length(input_matrix))] <- c("D.timepoint", "pseudot")
    
    # extract parameters
    npar <- which(colnames(input_matrix) == "D.timepoint") - 1
    data <- input_matrix[, 1:npar]
    Timepoints <- input_matrix[, npar + 1]
    Pseudotime <- input_matrix[, npar + 2]
    
    observedRoughness <- c() 
    r_perm_ofprotein <- data.frame() # 100次随机
    
    for (my_protein in 1:npar){
      a <- data[order(Pseudotime), my_protein]
      observedRoughness[my_protein] <- DeLorean.calc.roughness(a)
      names(observedRoughness)[my_protein] <- colnames(input_matrix)[my_protein]
      for (perms in 1:100) {
        set.seed(perms)
        a_perm <- sample(a)
        r_perm_ofprotein[perms, my_protein] <- DeLorean.calc.roughness(a_perm)
      }
    }
    
    PermMeanMatrix <- colMeans(r_perm_ofprotein, na.rm = FALSE, dims = 1)
    TEST <- t.test(x = observedRoughness, y = PermMeanMatrix, alternative = "less", paired = TRUE)$p.value
    
    TEST_list[lineage] <- TEST
    observedRoughnes_list[[lineage]] <- observedRoughness
  }
  return(list(p.value = mean(TEST_list), observedRoughness = observedRoughnes_list))
}