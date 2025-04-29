ANPELA_FLOWMAP <- function (mode = "multi", 
                            files, var.remove = c(), var.annotate = NULL, clustering.var= NULL, 
                            cluster.numbers = 100, transform = FALSE, cluster.mode = "hclust", 
                            distance.metric = "manhattan", 
                            minimum = 2, maximum = 10, save.folder = getwd(), subsamples = 100, 
                            name.sort = FALSE, out_folder_basename = NA, downsample = FALSE, 
                            seed.X = 1, savePDFs = FALSE, graph.out = "ForceDirected", 
                            which.palette = "bluered", dataset_name = NULL, exclude.pctile = NULL, target.pctile = NULL, 
                            target.number = NULL, target.percent = NULL, k = 10, umap.n.neighbors = 10, 
                            umap.n.components = 2,iter2 = 1000,umap_min_dist =0.01, ...) #name.sort = TRUE,subsamples = 200,maximum = 5,savePDFs = TRUE,distance.metric = "manhattan",added:transform = FALSE,dataset_name =NULL,iter2 = 1000
{ 
  
 
  
  #group fcs files by timepoints
  grouping <- sapply(names(files), function(x) gsub("[^0-9]", "", x))
  files <- split(files, grouping)#added
  cell_counts <- min(unlist(lapply(files, function(sublist) {lapply(sublist,  function(ff) nrow(ff))}),use.names = F))
  if (!subsamples == FALSE){
    if(cell_counts < subsamples){
      subsamples <- cell_counts
      cluster.numbers <- subsamples
    }
  }
  
  set.seed(seed.X)
  cat("Seed set to", seed.X, "\n")
  cat("Mode set to", mode, "\n")
  FLOWMAPR:::CheckSettings(mode, save.folder, var.remove, var.annotate,
                           clustering.var, cluster.numbers, cluster.mode, distance.metric,
                           minimum, maximum, which.palette, subsamples)
  if(is.null(clustering.var)){
    clustering.var <- files[[1]][[1]]@parameters@data[["name"]]
  }
  
  if (downsample) {
    FLOWMAPR:::CheckDownsampleSettings(exclude.pctile, target.pctile, 
                                       target.number, target.percent)
  }
  
  #setwd(save.folder)
  
  if (mode == "single") {
    check <- FLOWMAPR:::CheckModeSingle(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    runtype <- "SingleFLOWMAP"
    output.folder <- FLOWMAPR:::MakeOutFolder(runtype = runtype, k = k, 
                                              maximum = maximum, minimum = minimum, out_folder_basename = out_folder_basename)
    setwd(output.folder)
    FLOWMAPR:::PrintSummary(env = parent.frame())
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    orig.times <- FLOWMAPR:::ParseTimes(fcs.file.names, name.sort = name.sort)
    file.name <- fcs.file.names[1]
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", 
          "\n")
      fcs.files <- DownsampleFCS(fcs.file.names, clustering.var, 
                                 channel.annotate = var.annotate, channel.remove = var.remove, 
                                 exclude.pctile = exclude.pctile, target.pctile = target.pctile, 
                                 target.number = target.number, target.percent = target.percent, 
                                 transform = TRUE)
    } else {
      fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, 
                                channel.remove = var.remove, channel.annotate = var.annotate, 
                                subsamples = subsamples,transform = transform)#new added
      # fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, 
      #                           channel.remove = var.remove, channel.annotate = var.annotate, 
      #                           subsamples = subsamples)
    }
    file.clusters <- FLOWMAPR:::ClusterFCS(fcs.files = fcs.files, clustering.var = clustering.var, 
                                           numcluster = cluster.numbers, distance.metric = distance.metric, 
                                           cluster.mode = cluster.mode)
    global.file.clusters1 <<- file.clusters
    if (cluster.mode != "none") {
      cat("Upsampling all clusters to reflect Counts of entire file", 
          "\n")
      file.clusters <- FLOWMAPR:::Upsample(fcs.file.names, file.clusters, 
                                           fcs.files, var.remove, var.annotate, clustering.var)
    }
    results <- FLOWMAPR:::BuildFLOWMAPkNN(FLOWMAP.clusters = file.clusters, 
                                          k = k, min = minimum, max = maximum, distance.metric = distance.metric, 
                                          clustering.var = clustering.var)
    graph <- results$output.graph
    knn.out <- results$knn.out
    FLOWMAPR:::MakeFLOWMAPRFile(env = parent.frame())
    file.name <- paste(unlist(strsplit(basename(file.name), 
                                       "\\."))[1], "FLOW-MAP", sep = "_")
    ConvertToGraphML(output.graph = graph, file.name = file.name)
    if ("ForceDirected" %in% graph.out) {
      graph.xy <- FLOWMAPR:::RunForceDirectedLayout(graph = graph, 
                                                    mode = mode, file.name = file.name, orig.times = orig.times, 
                                                    which.palette = which.palette)
    }
    if ("UMAP" %in% graph.out) {
      FLOWMAPR:::RunUMAPlayout(knn.in = knn.out, graph = graph, file.clusters = file.clusters, 
                               umap_n_components = umap.n.components, k = k, 
                               clustering.var = clustering.var, file.name = file.name, 
                               umap_n_neighbors = umap.n.neighbors, mode = mode)
    }
  } else if (mode == "multi") {
    
    fcs.file.names <- lapply(files, function(group) {names(group)})#added
    orig.times <- ANPELA_MultiFolderParseTimes(fcs.file.names, 
                                               name.sort = name.sort)
    
    
    file.name <- fcs.file.names[[1]][1] #命名graphhml用
    
    if (downsample) {
      cat("Downsampling all files using density-based downsampling", 
          "\n")
      fcs.files <- FLOWMAPR::MultiDownsampleFCS(fcs.file.names, 
                                                clustering.var, channel.annotate = var.annotate, 
                                                channel.remove = var.remove, exclude.pctile = exclude.pctile, 
                                                target.pctile = target.pctile, target.number = target.number, 
                                                target.percent = target.percent, transform = TRUE)
    } else {
      fcs.files <- ANPELA_LoadMultiCleanFCS(files, var.remove, 
                                            var.annotate, subsamples = subsamples, transform = transform)#new added#name有点奇怪
    }
    
    fcs.files.conversion <- FLOWMAPR:::ConvertNumericLabel(fcs.files)#没搞懂干啥的
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    global.label.key <<- label.key
    file.clusters <- ANPELA_MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, 
                                            numcluster = cluster.numbers, distance.metric = distance.metric, 
                                            cluster.mode = cluster.mode)
    
    cat("Upsampling all clusters to reflect Counts of entire file", "\n")
    file.clusters <- ANPELA_MultiUpsample (files, file.clusters, 
                                           fcs.files, var.remove, var.annotate, clustering.var)
    
    
    remodel.FLOWMAP.clusters <- FLOWMAPR:::RemodelFLOWMAPClusterList(file.clusters, label.key)
    
    results <- ANPELA_BuildMultiFLOWMAPkNN(remodel.FLOWMAP.clusters, 
                                           k = k, min = minimum, max = maximum, distance.metric = distance.metric, 
                                           label.key = label.key, clustering.var = clustering.var)
    
    graph <- results$output.graph
    knn.out <- results$knn.out
    
    if ("ForceDirected" %in% graph.out) {
      graph.xy <- ANPELA_RunForceDirectedLayout(graph = graph, 
                                                mode = mode, orig.times = orig.times, 
                                                which.palette = which.palette,iter2 = iter2)
    }
    if ("UMAP" %in% graph.out) {
      graph.xy <- ANPELA_RunUMAPlayout(knn.in = knn.out, graph = graph, file.clusters = remodel.FLOWMAP.clusters, 
                                       umap_n_components = 2, k = k, clustering.var = clustering.var, 
                                       file.name = file.name, umap_n_neighbors = umap.n.neighbors, 
                                       mode = mode, umap_min_dist = umap_min_dist)
    }
  } else if (mode == "one") {
    check <- FLOWMAPR:::CheckModeOne(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    file.name <- fcs.file.names
    runtype <- "OneTimepoint"
    output.folder <- FLOWMAPR:::MakeOutFolder(runtype = runtype, k = k, 
                                              maximum = maximum, minimum = minimum)
    setwd(output.folder)
    FLOWMAPR:::PrintSummary(env = parent.frame())
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", 
          "\n")
      fcs.file <- DownsampleFCS(file.name, clustering.var, 
                                channel.annotate = var.annotate, channel.remove = var.remove, 
                                exclude.pctile = exclude.pctile, target.pctile = target.pctile, 
                                target.number = target.number, target.percent = target.percent, 
                                transform = TRUE)
    } else {
      fcs.file <- LoadCleanFCS(fcs.file.names = file.name, 
                               channel.remove = var.remove, channel.annotate = var.annotate, 
                               subsamples = subsamples,transform = transform)#new added
    }
    file.clusters <- FLOWMAPR:::ClusterFCS(fcs.files = fcs.file, clustering.var = clustering.var, 
                                           numcluster = cluster.numbers, distance.metric = distance.metric, 
                                           cluster.mode = cluster.mode)
    cat("Upsampling all clusters to reflect Counts of entire file", 
        "\n")
    file.clusters <- FLOWMAPR:::Upsample(file.name, file.clusters, 
                                         fcs.file, var.remove, var.annotate, clustering.var)
    results <- FLOWMAPR:::BuildFirstFLOWMAPkNN(FLOWMAP.clusters = file.clusters, 
                                               k = k, min = minimum, max = maximum, distance.metric = distance.metric, 
                                               clustering.var = clustering.var)
    output.graph <- results$output.graph
    output.graph <- FLOWMAPR:::AnnotateGraph(output.graph = output.graph, 
                                             FLOWMAP.clusters = file.clusters)
    graph <- output.graph
    knn.out <- list(indexes = results$indexes, distances = results$distances)
    FLOWMAPR:::MakeFLOWMAPRFile(env = parent.frame())
    file.name <- paste(unlist(strsplit(basename(file.name), 
                                       "\\."))[1], "FLOW-MAP", sep = "_")
    ConvertToGraphML(output.graph = graph, file.name = file.name)
    if ("ForceDirected" %in% graph.out) {
      graph.xy <- FLOWMAPR:::RunForceDirectedLayout(graph = graph, 
                                                    mode = mode, file.name = file.name, orig.times = orig.times, 
                                                    which.palette = which.palette)
    }
    if ("UMAP" %in% graph.out) {
      FLOWMAPR:::RunUMAPlayout(knn.in = knn.out, graph = graph, file.clusters = file.clusters, 
                               umap_n_components = 2, k = k, clustering.var = clustering.var, 
                               file.name = file.name, umap_n_neighbors = umap.n.neighbors, 
                               mode = mode)
    }
  } else if (mode == "static-multi") {
    check <- FLOWMAPR:::CheckModeSingle(files)
    cat("check", check, "\n")
    if (check[1]) {
      stop("Unknown 'files' format provided for specified mode!")
    }
    if (check[2] == "FCS") {
      fcs.file.names <- files
    }
    if (check[2] == "folder") {
      fcs.file.names <- GetFCSNames(folder = files, sort = name.sort)
    }
    file.name <- fcs.file.names[1]
    runtype <- "OneTimepoint-MultipleConditions"
    output.folder <- FLOWMAPR:::MakeOutFolder(runtype = runtype, k = k, 
                                              maximum = maximum, minimum = minimum)
    setwd(output.folder)
    FLOWMAPR:::PrintSummary(env = parent.frame())
    if (is.null(var.annotate)) {
      var.annotate <- ConstructVarAnnotate(file.name)
      assign("var.annotate", var.annotate, envir = .GlobalEnv)
    }
    if (downsample) {
      cat("Downsampling all files using SPADE downsampling", 
          "\n")
      fcs.files <- DownsampleFCS(fcs.file.names, clustering.var, 
                                 channel.annotate = var.annotate, channel.remove = var.remove, 
                                 exclude.pctile = exclude.pctile, target.pctile = target.pctile, 
                                 target.number = target.number, target.percent = target.percent, 
                                 transform = TRUE)
    } else {
      fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, 
                                channel.remove = var.remove, channel.annotate = var.annotate, 
                                subsamples = subsamples,transform = transform)#new added
      # fcs.files <- LoadCleanFCS(fcs.file.names = fcs.file.names, 
      #                           channel.remove = var.remove, channel.annotate = var.annotate, 
      #                           subsamples = subsamples)
    }
    fcs.files.old <- fcs.files
    fcs.files.list <- list()
    process.results <- FLOWMAPR:::ProcessConditions(fcs.files.old, 
                                                    fcs.file.names)
    label.key.special <- process.results$label.key.special
    fcs.files.list[[1]] <- process.results$fixed.files
    file.clusters <- FLOWMAPR:::MultiClusterFCS(list.of.files = fcs.files.list, 
                                                clustering.var = clustering.var, numcluster = cluster.numbers, 
                                                distance.metric = distance.metric, cluster.mode = cluster.mode)
    cat("Upsampling all clusters to reflect Counts of entire file", 
        "\n")
    temp.files.list <- list()
    temp.files.list[[1]] <- fcs.files
    temp.name.list <- list()
    temp.name.list[[1]] <- fcs.file.names
    file.clusters <-ANPELA_MultiUpsample(temp.name.list, file.clusters, 
                                         temp.files.list, var.remove, var.annotate, clustering.var)#new added
    # file.clusters <-FLOWMAPR::: MultiUpsample(temp.name.list, file.clusters, 
    #   temp.files.list, var.remove, var.annotate, clustering.var)
    remodel.FLOWMAP.clusters <- FLOWMAPR:::RemodelFLOWMAPClusterList(file.clusters)
    results <- FLOWMAPR:::BuildFirstMultiFLOWMAPkNN(list.of.FLOWMAP.clusters = remodel.FLOWMAP.clusters, 
                                                    k = maximum, min = minimum, max = maximum, distance.metric = distance.metric, 
                                                    clustering.var = clustering.var)
    output.graph <- results$output.graph
    output.graph <- FLOWMAPR:::AnnotateSpecialGraph(output.graph, remodel.FLOWMAP.clusters, 
                                                    label.key.special)
    graph <- output.graph
    knn.out <- list(indexes = results$indexes, distances = results$distances)
    FLOWMAPR:::MakeFLOWMAPRFile(env = parent.frame())
    file.name <- paste(unlist(strsplit(basename(file.name), 
                                       "\\."))[1], "FLOW-MAP", sep = "_")
    ConvertToGraphML(output.graph = graph, file.name = file.name)
    if ("ForceDirected" %in% graph.out) {
      graph.xy <-FLOWMAPR::: RunForceDirectedLayout(graph = graph, 
                                                    mode = mode, file.name = file.name, orig.times = orig.times, 
                                                    which.palette = which.palette)
    }
    if ("UMAP" %in% graph.out) {
      FLOWMAPR:::RunUMAPlayout(knn.in = knn.out, graph = graph, file.clusters = file.clusters, 
                               umap_n_components = 2, k = k, clustering.var = clustering.var, 
                               file.name = file.name, umap_n_neighbors = umap.n.neighbors, 
                               mode = mode)
    }
  } else {
    stop("Unknown mode!")
  }
  return(list(FLOWMAP_output = graph.xy, subsamples_time = subsamples))
}