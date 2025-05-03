dr_tsne <- function(data = data, tsneSeed = 42) {
  data <- as.matrix(data)
  if (is.numeric(tsneSeed)) set.seed(tsneSeed)
  mapped <- Rtsne::Rtsne(data, initial_dims = ncol(data), check_duplicates = FALSE, pca = TRUE, num_threads = parallel::detectCores()-1)$Y
  if (!is.null(mapped)) {
    colnames(mapped) <- paste("tsne", 1:2, sep = "_")
    rownames(mapped) <- row.names(data)
  }
  return(mapped)
}



# multicore_tSNE
dr_Rtsne.multicore <- function(X, dims = 2, initial_dims = 50,
                               perplexity = 350, theta = 0.5,
                               check_duplicates = TRUE,
                               pca = TRUE, max_iter = 1000, num_threads = 2, verbose = FALSE, is_distance = all(class(X)=='dist'),
                               pca_center = TRUE, pca_scale = FALSE, ...) {
  X <- model.matrix(~ . -1, data = X)
  if(dims != 2) { stop("Only 2d output is supported due to its c++ implemenation!") }
  if (!is.logical(is_distance)) { stop("is_distance should be a logical variable") }
  if (!is.numeric(theta) || (theta<0.0) || (theta>1.0) ) { stop("Incorrect theta.") }
  if (all(!is_distance & nrow(X) - 1 < 3 * perplexity)) { stop("Perplexity is too large.") }
  if(!any(class(X) %in% c("dist","matrix"))) { stop("Input X must be matrix or a dist object") }
  if (!(max_iter>0)) { stop("Incorrect number of iterations.") }

  if(is_distance) {
    # convert to matrix if dist object was supplied
    if(!is.matrix(X)) {
      X <- as.matrix(X)
    }
    if(nrow(X)!=ncol(X)) { stop("Input is not an accepted distance matrix") }
  }

  if (!(is.logical(pca_center) && is.logical(pca_scale)) ) { stop("pca_center and pca_scale should be TRUE or FALSE")}

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(initial_dims) || initial_dims<=0) { stop("Incorrect initial dimensionality.")}

  # Apply PCA
  if (pca & !is_distance) {
    pca_result <- prcomp(X,retx=TRUE,center = pca_center, scale. = pca_scale)
    X <- pca_result$x[,1:min(initial_dims,ncol(pca_result$x))]
  }
  if (check_duplicates & !is_distance){
    if (any(duplicated(X))) { stop("Remove duplicates before running TSNE.") }
  }
  # Compute Squared distance if we are using exact TSNE
  if (is_distance & theta==0.0) {
    X <- X^2
  }

  msg <- capture.output(
    res <-  Rtsne.multicore:::Rtsne_cpp(X, dims, perplexity, theta,num_threads, max_iter, is_distance)
  )
  if(verbose)
    print(msg)

  return(res)
}

# run_slingshot ----------------------------------------------------------
run_slingshot <- function(D, drMethod = c("tSNE", "PCA", "diffMaps","FLOWMAP"),
                          dataset_name = NULL,clustering.var = NULL, subsamples = 100, cluster.numbers = 100){
  # 1st step => dim reduction
  data_drMethod <- switch (drMethod,
                           tSNE = dr_Rtsne.multicore(X = D$expr, dims = 2)$Y,
                           PCA = prcomp(D$expr, scale. = FALSE, rank. =2)$x,
                           diffMaps = {
                             DifMap_output <- destiny::DiffusionMap(D$expr)
                             data.frame(cbind(DifMap_output$DC1, DifMap_output$DC2))},
                           FLOWMAP ={
                           
                             FLOWMAP_output_list <<- try(BBmisc::suppressAll(do.call(ANPELA_FLOWMAP,
                                                                                     list(files = D$expr, clustering.var = clustering.var,
                                                                                          subsamples = subsamples, cluster.numbers = cluster.numbers))),
                                                         silent = T)
                             if (class(FLOWMAP_output_list) == "try-error"){
                               data.frame()
                             } else {
                               FLOWMAP_output <<- FLOWMAP_output_list$FLOWMAP_output
                               data.frame(cbind(FLOWMAP_output[,"x"], FLOWMAP_output[,"y"]))
                             }
                           }
  )

  # 2nd step => trajectory inference
  # Clustering
  set.seed(1)
  library(mclust, quietly = TRUE)
  BBmisc::suppressAll(cl1 <- mclust::Mclust(data_drMethod)$classification)

  if (drMethod == "FLOWMAP"){
    cl1_time <- cbind(data.frame(cl1),data.frame(FLOWMAP_output$Time))
  } else {
    cl1_time <- cbind(data.frame(cl1),data.frame(D[["timepoint"]]))
  }
  colnames(cl1_time) <- c("cl1","timepoint")
  time_zero_data <- cl1_time[cl1_time$timepoint == sort(unique(cl1_time$timepoint))[1], ]
  most_common_cluster <- names(which.max(table(time_zero_data$cl1)))

  # Mapping
  lin1 <- slingshot::getLineages(data_drMethod, cl1, start.clus =most_common_cluster)
  # Construct smooth curves
  crv1 <- slingshot::getCurves(lin1)#,approx_points = 1000
  # this is your pseudotime
  pseudot <- slingshot::slingPseudotime(crv1)

  # 3rd step => Save a table with Normalized_values [1:limit], Experimental time, Pseudotime

  to_return <- list()
  to_return[["dimRed"]] <- data_drMethod
  to_return[["pseudotime"]] <- pseudot
  test <- try(crv1@metadata, silent = T)
  if (class(test) == "try-error"){
    for (i in 1:length(crv1@curves)) {
      curve <- crv1@curves[[i]]$s
      to_return[["trajectory"]] [[paste0("curve", i)]] <- curve
    }
  } else {
    for (i in 1:length(crv1@metadata[["curves"]])) {
      curve <- crv1@metadata[["curves"]][[i]]$s
      to_return[["trajectory"]] [[paste0("curve", i)]] <- curve

    }
  }

  if (class(test) == "try-error"){
    to_return[["lineages"]] <- length(crv1@lineages)
  } else {
    to_return[["lineages"]] <- length(crv1@metadata[["lineages"]])
  }

  to_return[["cr_method"]] <- "run_slingshot"
  to_return[["dr_method"]] <- drMethod
  to_return[["crv1"]] <- crv1
  if (drMethod == "FLOWMAP") {
    to_return[["timepoint"]] <- as.matrix(FLOWMAP_output$Time)
    to_return[["condition"]] <- as.matrix(FLOWMAP_output$Condition)
    to_return[["expr"]] <- FLOWMAP_output[2:(which(colnames(FLOWMAP_output) == "Time") - 1)]
    to_return[["subsamples_time"]] <- FLOWMAP_output_list$subsamples_time
  }
  return(to_return)
}


# run_scorpius -----------------------------------------------------------
run_scorpius <- function(D, drMethod = c("distPear", "distSpear", "distEucl", "distManh")) {
  # 1st step => dim reduction
  space <- switch (drMethod,
                   distPear = SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "pearson"),
                   distSpear = SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "spearman"),
                   distEucl = SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "euclidean"),
                   distManh = SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "manhattan")
  )

  # 2nd step => trajectory
  traj <- SCORPIUS::infer_trajectory(space)
  # this is your pseudot
  pseudot<- traj$time

  # 3rd step => Save a table with the results
  to_return <- list()
  to_return[["dimRed"]] <- space
  to_return[["pseudotime"]] <- pseudot
  to_return[["trajectory"]] <- traj$path
  to_return[["lineages"]] <- 1
  to_return[["cr_method"]] <- "run_scorpius"
  to_return[["dr_method"]] <- drMethod
  return(to_return)
}

# run_prinCurves ---------------------------------------------------------
run_prinCurves <- function(D, drMethod = c("tSNE", "diffMaps")){
  # 1st step => dim reduction
  data_drMethod <- switch (drMethod,
                           tSNE = dr_Rtsne.multicore(D$expr, dims = 2)$Y,
                           diffMaps = {
                             DifMap_output <- destiny::DiffusionMap(D$expr)
                             data_drMethod <- data.frame(cbind(DifMap_output$DC1, DifMap_output$DC2))}
  )

  # 2nd step => trajectory
  fit1 <- princurve::principal_curve(x = as.matrix(data_drMethod),
                                     maxit = 10,
                                     smoother = "smooth_spline")
  pseudot <- fit1$lambda

  # 3rd step => Save a table with the Abundances, Experimental time, Pseudotime values and Pseudotime trajectory
  to_return <- list()
  to_return[["dimRed"]] <- data_drMethod
  to_return[["pseudotime"]] <- pseudot
  to_return[["trajectory"]] <- fit1$s
  to_return[["lineages"]] <- 1
  to_return[["cr_method"]] <- "run_prinCurves"
  to_return[["dr_method"]] <- drMethod
  return(to_return)
}

TI <- function(D, method = c("scorpius_distSpear", "scorpius_distEucl", "scorpius_distPear", "scorpius_distManh",
                             "slingshot_tSNE", "prinCurves_tSNE", "slingshot_FLOWMAP", #flowmap
                             "slingshot_PCA", "slingshot_diffMaps", "prinCurves_diffMaps"),
               dataset_name = NULL, clustering.var = NULL, ...){
  set.seed(456)
  switch (method,
          slingshot_FLOWMAP = suppressMessages(run_slingshot(D, drMethod = "FLOWMAP",dataset_name = dataset_name, clustering.var =clustering.var)), # 0.00
          slingshot_tSNE = suppressMessages(run_slingshot(D, drMethod = "tSNE")), # 19.70=>14.81
          slingshot_PCA = suppressMessages(run_slingshot(D, drMethod = "PCA")), # 56.68=>52.74
          slingshot_diffMaps = suppressMessages(run_slingshot(D, drMethod = "diffMaps")), # 53.69

          scorpius_distPear = run_scorpius(D, drMethod = "distPear"), # 4.35=>2.62
          scorpius_distSpear = run_scorpius(D, drMethod = "distSpear"), #3.41=>2.45
          scorpius_distEucl = run_scorpius(D, drMethod = "distEucl"), # 4.03=>2.56
          scorpius_distManh = run_scorpius(D, drMethod = "distManh"), # 4.68=>2.58

          prinCurves_tSNE = run_prinCurves(D, drMethod = "tSNE"), # 24.36=>18.14
          prinCurves_diffMaps = run_prinCurves(D, drMethod = "diffMaps") # 44.81=>36.35
  )
}
