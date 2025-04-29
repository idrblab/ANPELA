ANPELA_FLOWMAP <- function (mode = "multi", 
                            files, var.remove = c(), var.annotate = NULL, clustering.var = NULL, 
                            cluster.numbers = 100, transform = TRUE, cluster.mode = "hclust", 
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
  
  #filter out empty fcs files
  files <- lapply(files, function(sublist) {
    filtered_sublist <- lapply(sublist, function(ff) {
      if (nrow(ff) > 0) {
        return(ff)
      }
    })
    filtered_sublist <- Filter(Negate(is.null), filtered_sublist)
    return(filtered_sublist)
  })
  files <-  Filter(function(x) length(x) > 0, files)
  

  cell_counts <- min(unlist(lapply(files, function(sublist) {
    lapply(sublist, function(ff) nrow(ff))
  }), use.names = FALSE))
  
  if (!subsamples == FALSE){
    if(cell_counts < subsamples){
      subsamples <- cell_counts
      cluster.numbers <- subsamples
    }
  }
  
  set.seed(seed.X)
  cat("Seed set to", seed.X, "\n")
  cat("Mode set to", mode, "\n")
  
  if(is.null(clustering.var)){
    clustering.var <- files[[1]][[1]]@parameters@data[["name"]]
  }
  
  FLOWMAPR:::CheckSettings(mode, save.folder, var.remove, var.annotate,
                           clustering.var, cluster.numbers, cluster.mode, distance.metric,
                           minimum, maximum, which.palette, subsamples)

  if (downsample) {
    FLOWMAPR:::CheckDownsampleSettings(exclude.pctile, target.pctile, 
                                       target.number, target.percent)
  }
  #setwd(save.folder)
  if (mode == "multi") {
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
                                                target.percent = target.percent, transform = transform)
    } else {
      fcs.files <- ANPELA_LoadMultiCleanFCS(files, var.remove, 
                                            var.annotate, subsamples = subsamples, transform = transform)
    }

    fcs.files.conversion <- ANPELA_ConvertNumericLabel(fcs.files)#把condition用数字表示
    
    fixed.fcs.files <- fcs.files.conversion$fixed.list.FCS.files
    label.key <- fcs.files.conversion$label.key
    global.label.key <<- label.key
    file.clusters <- ANPELA_MultiClusterFCS(fixed.fcs.files, clustering.var = clustering.var, 
                                            numcluster = cluster.numbers, distance.metric = distance.metric, 
                                            cluster.mode = cluster.mode)

    cat("Upsampling all clusters to reflect Counts of entire file", "\n")
    file.clusters <- ANPELA_MultiUpsample (files, file.clusters, 
                                           fcs.files, var.remove, var.annotate, clustering.var)
    
    remodel.FLOWMAP.clusters <- ANPELA_RemodelFLOWMAPClusterList(file.clusters, label.key)
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
  } else {
    stop("Unknown mode!")
  }
  return(list(FLOWMAP_output = graph.xy, subsamples_time = subsamples))
}