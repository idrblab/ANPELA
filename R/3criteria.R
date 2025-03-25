sub_cluster <- function(data, FlowSeed) {
  res <- try(suppressMessages(data_cluster(data = data[,1:(length(data)-3)], method = "FlowSOM", FlowSOM_k = 2, FlowSeed = FlowSeed)), 
             silent = T)
  return(res)
}

feature_selection <- function(data) {
  data3 <- dplyr::tibble(data)
  data4 <- data_ttest(data = data3)
  selected_markers <- colnames(data4)[-c(which(colnames(data4) == "filename"), which(colnames(data4) == "condition"),which(colnames(data4) == "cluster"))]
  # pvalue_res <- sapply(selected_markers, pvalue_method, data = data4, method = "t.test", paired = FALSE)
  pvalue_res <- sapply(selected_markers, pvalue_method, data = data4)
  return(names(which(pvalue_res >= 0.05)))
}

data_ttest <- function (data) {
  res <- dplyr::group_by(data, filename, condition) %>% 
    dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
  return(res)
}

# pvalue_method <- function (marker, data, 
#                            method = c("wilcox.test","t.test"), paired = FALSE) {
#   data <- data.frame(marker_expr = data[, marker], condition = data[,"condition"])
#   colnames(data) <- c("marker_expr", "condition")
#   res <- ggpubr::compare_means(marker_expr~condition, data, 
#                                method = method, paired = paired, p.adjust.method = "none")$p
# 
# }
pvalue_method <- function (marker, data) {
  data <- data.frame(marker_expr = data[, marker], condition = data[,"condition"])
  colnames(data) <- c("marker_expr", "condition")
  res <- t.test(marker_expr~condition, data)$p.value
}

subsampling <- function (data, size) {
  test.fold <- sampling::strata(c("condition"), size = size, method = "srswor", data = data)[,2]
  return(test.fold)
}

accuracy <- function(result, label) {
  total_num <- length(label)
  cluster_counter <- unique(result) # 1 2
  original_counter <- unique(label) # Biopsy PBMC
  t <- NULL
  
  for (k in cluster_counter) {
    p_k <- NULL
    for (j in original_counter) {
      count <- 0
      for (i in 1:length(result)) {
        if (result[i] == k && label[i] == j) {
          count <- count + 1
        }
      }
      p_k <- c(p_k, count)
    }
    temp_t <- max(p_k)
    t <- c(t, temp_t)
  }
  
  res <- sum(t)/total_num
  return(res)
}

CSvalue <- function (data = data, subdata = test.fold) {
  DEG <- list()
  for (j in 1:length(subdata)) { # 1:3
    data2 <- data[subdata[[j]],]
    data3 <- tibble(data2)
    data4 <- data_ttest(data = data3)
    
    selected_markers <- colnames(data4)[-c(which(colnames(data4) == "filename"), which(colnames(data4) == "condition"),which(colnames(data4) == "cluster"))]
    # pvalue_res <- sapply(selected_markers, pvalue_method, data = data4, method = "t.test", paired = FALSE)
    pvalue_res <- sapply(selected_markers, pvalue_method, data = data4)
    DEG[[letters[j]]] <- names(which(pvalue_res < 0.05))
  }
  VennList <- systemPipeR::overLapper(setlist = DEG, sep = "", type = "vennsets")@vennlist
  con.score <- 0
  for (i in 1:length(VennList)) {
    insect.n <- nchar(names(VennList[i]))
    if (insect.n < 2) next
    num.i <- 2^(insect.n - 2) * length(VennList[[i]])
    con.score <- con.score + num.i
  }
  return(list(con.score = con.score, DEG = DEG))
}

CWrel <- function(all_genelist = all_genelist, Y, n) {
  num <- 3
  partial_genelist <- vector("list", length(all_genelist))
  for (j in 1:length(all_genelist)) {
    partial_genelist[[j]] <- all_genelist[[j]][1:n]
  }
  
  partial_genelist <- as.list(partial_genelist)
  partial_genelist <- table(unlist(partial_genelist))
  
  Ff_sum <- 0
  for (k in 1:length(partial_genelist)) {
    Ff_sum <- Ff_sum + partial_genelist[k] * (partial_genelist[k] - 1)
  }
  Ff_sum <- as.numeric(Ff_sum)
  #Y <- ncol(data)
  N <- sum(partial_genelist)
  D <- N %% Y
  H <- N %% num
  
  CWrel <- (Y*(N-D+Ff_sum)-N^2+D^2)/(Y*(H^2+num*(N-H)-D)-N^2+D^2)
  return(CWrel)
}

CWvalue <- function (data = data, subdata = test.fold, top = 10) {
  DEG <- list()
  for (j in 1:length(subdata)) { # 1:3
    data2 <- data[subdata[[j]],]
    data3 <- tibble(data2)
    data4 <- data_ttest(data = data3)
    
    selected_markers <- colnames(data4)[-c(which(colnames(data4) == "filename"), which(colnames(data4) == "condition"),which(colnames(data4) == "cluster"))]
    # pvalue_res <- sapply(selected_markers, pvalue_method, data = data4, method = "t.test", paired = FALSE)
    pvalue_res <- sapply(selected_markers, pvalue_method, data = data4)
    DEG[[letters[j]]] <- names(sort(pvalue_res)[1:top])
  }
  CW_value <- try(CWrel(DEG, Y = (ncol(data)-3), n = length(DEG[[1]])))
  return(list(DEG = DEG, CW_value = CW_value))
}

# 1. Accuracy ----------------------------------------------------------------

F1_score <- function(test = test, KNN_res = KNN_res, label = condition_label) {
  F1 <- list()
  for (j in 1:2) {
    positive <- as.character(unique(label)[j])
    f1_score <- c()
    for (i in 1:length(test)) {
      if (length(unique(test[[i]]$condition)) < 2) next
      confusion_matrix <- MLmetrics::ConfusionDF(KNN_res[[i]], test[[i]]$condition)
      
      TP <- as.numeric(subset(confusion_matrix, y_true == positive & y_pred == positive)["Freq"])
      FP <- as.numeric(subset(confusion_matrix, y_true != positive & y_pred == positive)["Freq"])
      FN <- as.numeric(subset(confusion_matrix, y_true == positive & y_pred != positive)["Freq"])
      
      Precision <- TP/(TP+FP)
      Recall <- TP/(TP+FN)
      f1_score <- c(f1_score, 2 * (Precision * Recall) / (Precision + Recall))
      names(f1_score)[length(f1_score)] <- paste("cluster", i)
    }
    F1[[j]] <- f1_score
  }
  return(F1)
}

AUC <- function(test = test, KNN_res = KNN_res) {
  auc <- c()
  case_problist <- list()
  for (i in 1:length(test)) {
    if (length(unique(test[[i]]$condition)) < 2) next
    case_prob <- attributes(KNN_res[[i]])$prob
    control_label <- levels(KNN_res[[i]])[1]
    case_prob[which(KNN_res[[i]] == control_label)] <- 1- case_prob[which(KNN_res[[i]] == control_label)]
    auc <- c(auc, as.numeric(pROC::roc(test[[i]]$condition, case_prob)[["auc"]]))
    names(auc)[length(auc)] <- paste("cluster", i)
    case_problist[[i]] <- case_prob
  }
  return(list(auc = auc, case_problist = case_problist))
}
# 2. Precision ----------------------------------------------------------------
RI <- function(data = subdata_cluster_DEG, sub_cluster_label, FlowSeed = 40) {
  ri <- c()
  for (i in 1:length(data)) {
    if (class(sub_cluster_label[[i]]) == "NULL" | class(sub_cluster_label[[i]]) == "try-error") next
    x <- as.numeric(data[[i]]$condition)
    y <- as.numeric(sub_cluster_label[[i]])
    ri <- c(ri, as.numeric(fossil::rand.index(x, y)))
    names(ri)[length(ri)] <- paste("cluster", i)
  }
  return(ri)
}

Purity <- function(data = subdata_cluster_DEG, sub_cluster_label, FlowSeed = 40) {
  purity <- c()
  for (i in 1:length(data)) {
    if (class(sub_cluster_label[[i]]) == "NULL" | class(sub_cluster_label[[i]]) == "try-error") next
    purity <- c(purity, accuracy(sub_cluster_label[[i]], data[[i]]$condition))
    names(purity)[length(purity)] <- paste("cluster", i)
  }
  return(purity)
}

# 3. Robustness -----------------------------------------------------------

CS_pre <- function(subdata_cluster = subdata_cluster) {
  data <- list()
  test.fold <- list()
  for (i in 1:length(subdata_cluster)) {
    
    data[[i]] <- subdata_cluster[[i]]
    
    set.seed(123)
    test.fold[[i]] <- list()
    
    size <- as.numeric(table(data[[i]]$condition))/3
    test.fold1 <- try(subsampling(data = data[[i]], size = size), silent = T)
    if (class(test.fold1) == "try-error") next
    test.fold[[i]][[1]] <- test.fold1
    
    data.2 <- data[[i]][-test.fold1, ]
    test.fold2 <- try(subsampling(data = data.2, size = size), silent = T)
    if (class(test.fold2) == "try-error") next
    test.fold2 <- match(row.names(data.2)[test.fold2], row.names(data[[i]]))
    test.fold[[i]][[2]] <- test.fold2
    
    test.fold[[i]][[3]] <- (1:nrow(data[[i]]))[-c(test.fold1, test.fold2)]
  }
  return(list(data = data, test.fold = test.fold))
}

CSfun <- function(CS_preres = CS_preres) {
  consistency <- c()
  DEGlist <- list()
  for (i in 1:length(CS_preres$data)) {
    CS_res <- try(CSvalue(data = CS_preres$data[[i]], subdata = CS_preres$test.fold[[i]]), silent = T)
    if(class(CS_res) == "try-error") next
    
    consistency <- c(consistency, CS_res$con.score)
    names(consistency)[length(consistency)] <- paste("cluster", i)
    # DEGlist[[paste("cluster", i)]] <- CS_res$DEG
    DEGlist[[i]] <- CS_res$DEG
  }
  return(list(consistency = consistency, DEGlist = DEGlist))
}
CWfun <- function(CS_preres = CS_preres, top = 10) {
  consistency <- c()
  DEGlist <- list()
  for (i in 1:length(CS_preres$data)) {
    CW_res <- try(CWvalue(data = CS_preres$data[[i]], subdata = CS_preres$test.fold[[i]], top = top), silent = T)
    if(class(CW_res) == "try-error") next
    
    consistency <- c(consistency, CW_res$CW_value)
    names(consistency)[length(consistency)] <- paste("cluster", i)
    # DEGlist[[paste("cluster", i)]] <- CW_res$DEG
    DEGlist[[i]] <- CW_res$DEG
  }
  return(list(consistency = consistency, DEGlist = DEGlist))
}

# 4. Biological Meaning -----------------------------------------------------------
AP2_Recall <- function(data_with_cluster = data_with_cluster, known_marker = known_marker) {
  non_DEG <- feature_selection(data_with_cluster)
  DEG <- colnames(data_with_cluster)[- which(colnames(data_with_cluster) %in% c(non_DEG, "filename", "condition", "cluster"))]
  recall <- length(intersect(DEG, known_marker)) / length(known_marker)
  return(recall)
}


