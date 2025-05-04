#ANPELA_MultiUpsample---------------------------------------------------
ANPELA_MultiUpsample <- function (file.names, FLOWMAP.clusters, fcs.files, var.remove,
                                  var.annotate, clustering.var, transform = TRUE)
{
  fixed.FLOWMAP.clusters <- FLOWMAP.clusters
  for (t in 1:length(FLOWMAP.clusters)) {
    for (f in 1:length(FLOWMAP.clusters[[t]]$cluster.counts)) {
      original.fcs.file <- ANPELA_LoadCleanFCS(file.names[[t]][f],
                                               channel.remove = var.remove, channel.annotate = var.annotate,
                                               subsamples = FALSE,transform = transform)
      original.fcs.file <- original.fcs.file[[1]]
      original.fcs.file <- original.fcs.file[, clustering.var]
      downsample.fcs.file <- fcs.files[[t]][[f]]
      downsample.fcs.file <- downsample.fcs.file[, clustering.var]
      cluster.assign <- as.integer(as.matrix(FLOWMAP.clusters[[t]]$cell.assgn[[f]]))
      downsample.fcs.file$cluster <- cluster.assign
      cluster.medians <-  downsample.fcs.file %>% dplyr::group_by(cluster) %>%
        dplyr::summarise_all(median)
      all.cells.assign <- RANN::nn2(cluster.medians[, 2:ncol(cluster.medians)],
                                    original.fcs.file,
                                    k = 1)
      fixed.counts <- table(all.cells.assign[["nn.idx"]])
      fixed.counts <- as.data.frame(as.matrix(fixed.counts))

      if (length(fixed.counts) < length(cluster.medians$cluster)) {
        missing.cluster <- setdiff(cluster.medians$cluster,
                                   all.cells.assign$nn.idx)
        new.cluster <- data.frame(V1 = rep(0,length(missing.cluster)), row.names = missing.cluster)
        fixed.counts <- rbind(fixed.counts, new.cluster)
      }

      colnames(fixed.counts) <- c("Counts")
      fixed.FLOWMAP.clusters[[t]]$cluster.counts[[f]]$Counts <- fixed.counts
      rm(original.fcs.file, downsample.fcs.file, cluster.assign,
         fixed.counts)
    }
  }
  return(fixed.FLOWMAP.clusters)
}

#ANPELA_MakeOutFolder---------------------------------------------------
ANPELA_MakeOutFolder <- function (dataset_name,runtype, maximum, minimum, k = "", out_folder_basename = NA, graph.out, iter2)
{
  output <- graph.out#new added
  if (is.na(out_folder_basename)) {
    name <- gsub(" ", "_", Sys.time(), fixed = TRUE)
    name <- gsub(":", ".", name, fixed = TRUE)
    #new added
    if(is.null(dataset_name)){
    	output.folder <- paste("iter2",iter2,output,"max", maximum, "_k", k, name,
                                runtype, "run", sep = "_")
    } else {output.folder <- paste(dataset_name,"iter2",iter2, "max", maximum, "_k", k,
                                   sep = "_")
    }

    print(output.folder)
    dir.create(output.folder)
    cat("output.folder is", output.folder, "\n")
    return(output.folder)
  } else {
    name <- gsub(" ", "_", Sys.time(), fixed = TRUE)
    name <- gsub(":", ".", name, fixed = TRUE)
    output.folder <- paste(out_folder_basename, name, runtype,
                           "run", sep = "_")
    print(output.folder)
    dir.create(output.folder)
    cat("output.folder is", output.folder, "\n")
    return(output.folder)
  }
}

# ANPELA_RunForceDirectedLayout ----------------------------------------------------

ANPELA_RunForceDirectedLayout <- function (mode, graph, orig.times = NULL, which.palette = NULL,iter2 = 1000) ##new added iter2 = 1000;delete:file.name
{
  graph.xy <- ANPELA_ForceDirectedXY(graph = graph,iter2 = iter2)#new added

  if (mode != "one" && mode != "one-special") {
    fixed.graph <- FLOWMAPR:::ConvertOrigTime(graph.xy, orig.times)
  } else {
    fixed.graph <- graph.xy
  }

  all.attr <- igraph::get.vertex.attribute(fixed.graph)
  this.df <- c()
  for (x in 1:length(all.attr)) {
    this.df <- cbind(this.df, all.attr[[x]])
    this.df <- as.data.frame(this.df)
  }
  colnames(this.df) <- names(all.attr)

  return(this.df)
}

#ANPELA_ForceDirectedXY---------------------------------------------------
ANPELA_ForceDirectedXY <- function (graph,iter2 = 1000)
{

  igraph::V(graph)$size <- rep(20, igraph::vcount(graph))
  force.graph1 <- ANPELA_layout_forceatlas2(graph, iter = 10000,
                                           stopping.tolerance = 0.001, prevent.overlap = FALSE)
  graph.with.xy <- graph
  igraph::V(graph.with.xy)$x <- force.graph1$lay[, 1]
  igraph::V(graph.with.xy)$y <- force.graph1$lay[, 2]

  force.graph2 <- ANPELA_layout_forceatlas2(graph.with.xy,
                                           iter = iter2, stopping.tolerance = 0.001, prevent.overlap = TRUE)#new added

  igraph::V(graph.with.xy)$x <- force.graph2$lay[, 1]
  igraph::V(graph.with.xy)$y <- force.graph2$lay[, 2]
  return(graph.with.xy)
}


#ANPELA_MultiFolderParseTimes---------------------------------------------------
ANPELA_MultiFolderParseTimes <- function (fcs.file.names, name.sort)
{
  times <- c()
  for (i in 1:length(fcs.file.names)) {
    if (length(fcs.file.names[[i]]) > 1) {
      warning("Multiple FCS files per timepoint, only using the first for time label!")
      times <- c(times, FLOWMAPR:::ParseTimes(fcs.file.names[[i]][1],
                                              name.sort))
    } else {
      times <- c(times, FLOWMAPR:::ParseTimes(fcs.file.names[[i]],
                                              name.sort))
    }
  }
  return(times)
}

#ANPELA_LoadMultiCleanFCS----------------------------------------------------
ANPELA_LoadMultiCleanFCS <- function (list.of.file.names, channel.remove, channel.annotate,
                                      subsamples = 1000, transform = TRUE)
{
  list.of.FCS.files <- list()
  subsamp.orig <- subsamples
  for (t in 1:length(list.of.file.names)) {
    fcs.file.names <- list.of.file.names[[t]]
    if (length(subsamp.orig) == 1) {
      if (subsamp.orig != FALSE) {
        cat("Subsampling all files to:", subsamp.orig,
            "\n")
        subsample.new <- rep(subsamp.orig, times = length(fcs.file.names))
        subsamples <- subsample.new
      }
    } else if (length(subsamp.orig) > 1) {
      subsamples <- subsamp.orig[[t]]
    }
    list.of.FCS.files[[t]] <- ANPELA_LoadCleanFCS(fcs.file.names,
                                                  channel.remove, channel.annotate, subsamples, transform)#added

    f.names <- c()
    for (i in 1:length(list.of.FCS.files[[t]])) {
      Time <- rep(as.numeric(t), times = dim(list.of.FCS.files[[t]][[i]])[1])
      list.of.FCS.files[[t]][[i]] <- cbind.data.frame(list.of.FCS.files[[t]][[i]],
                                                      Time, stringsAsFactors = FALSE)
      this.name <- names(fcs.file.names)[i]
      this.name <- gsub(".fcs", "", this.name)
      if (grepl("-", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "-"))[1]
      }
      if (grepl("\\.", this.name)) {
        this.name <- unlist(strsplit(this.name, split = "\\."))[1]
      }
      Condition <- rep(this.name, times = dim(list.of.FCS.files[[t]][[i]])[1])
      list.of.FCS.files[[t]][[i]] <- cbind.data.frame(list.of.FCS.files[[t]][[i]],
                                                      Condition, stringsAsFactors = FALSE)
      f.names <- c(f.names, this.name)
      rm(Time, Condition)
    }
    names(list.of.FCS.files[[t]]) <- f.names
  }
  return(list.of.FCS.files)
}

#ANPELA_LoadCleanFCS---------------------------------------------------
ANPELA_LoadCleanFCS <- function (fcs.file.names, channel.remove, channel.annotate,
                                 subsamples = 1000, transform = TRUE)
{
  clean.fcs.files <- list()
  print("Subamples variable: ")
  print(subsamples)
  if (length(subsamples) == 1) {
    if (subsamples != FALSE) {
      cat("Subsampling all files to:", subsamples, "\n")
      subsample.new <- rep(subsamples, times = length(fcs.file.names))
      subsamples <- subsample.new
    }
  } else {
    subsamples <- subsamples
  }
  for (i in 1:length(fcs.file.names)) {
    current.file <- names(fcs.file.names)[i]

    cat("Reading FCS file data from:", current.file, "\n")
    if (length(subsamples) == 1) {
      if (subsamples == FALSE) {
        fcs.file <-fcs.file.names[[i]]

        if(class(fcs.file) == "data.frame"){
          fcs.file <- as.data.frame(fcs.file)
        } else {
          fcs.file <- as.data.frame(exprs(fcs.file))
        }
        print(nrow(fcs.file))
      } else {
        cat("Subsampling", current.file, "to", subsamples[i],
            "cells\n")
        fcs.file <-fcs.file.names[[i]]

        fcs.file <- as.data.frame(exprs(fcs.file))
        subsample.ids <- runif(subsamples[i], min = 1,
                               max = nrow(fcs.file))
        fcs.file <- fcs.file[subsample.ids, ]
        rownames(fcs.file) <- c(1:nrow(fcs.file))
      }
    } else {
      cat("Subsampling", current.file, "to", subsamples[i],
          "cells\n")
      fcs.file <-fcs.file.names[[i]]

      fcs.file <- as.data.frame(exprs(fcs.file))
      subsample.ids <- runif(subsamples[i], min = 1, max = nrow(fcs.file))
      fcs.file <- fcs.file[subsample.ids, ]
      rownames(fcs.file) <- c(1:nrow(fcs.file))
    }
    global.colnames.pre.fix <<- colnames(fcs.file)

    cat("Removing unnecessary channel names from:", current.file,
        "\n")
    fcs.file <- subset(fcs.file, select = colnames(fcs.file)[!colnames(fcs.file) %in%
                                                               channel.remove])
    if (transform) {
	  cat("Transforming data from:", current.file, "\n")
      fcs.file <- apply(fcs.file, 2, FLOWMAPR:::Asinh)
    }

    fcs.file <- as.data.frame(fcs.file)
    fcs.file <- FLOWMAPR:::RemoveExistingTimeVar(fcs.file)
    clean.fcs.files[[i]] <- fcs.file
    rm(fcs.file)
  }
  return(clean.fcs.files)
}

#ANPELA_ConvertNumericLabel----------------------------------------------------
ANPELA_ConvertNumericLabel <- function (list.of.clean.FCS.files.with.labels)
{
  label.key <- list()
  fixed.list.FCS.files <- list()
  for (i in 1:length(list.of.clean.FCS.files.with.labels)) {
    label.key[[i]] <- names(list.of.clean.FCS.files.with.labels[[i]])
    tmp.label.fcs <- list.of.clean.FCS.files.with.labels[[i]]
    for (n in 1:length(tmp.label.fcs)) {
      inds <- matrix(integer(0))
      for (l in 1:length(label.key[[i]])){
        if (length(inds) == 0){
          inds <- which(tmp.label.fcs[[n]] == label.key[[i]][l],
                        arr.ind = TRUE)
        } else {
          next
        }
      }
      col.ind <- unname(inds[1, 2])
      tmp.label.fcs[[n]][, col.ind] <- as.numeric(n)
    }
    fixed.list.FCS.files[[i]] <- tmp.label.fcs
  }
  names(fixed.list.FCS.files) <- names(list.of.clean.FCS.files.with.labels)
  return(list(fixed.list.FCS.files = fixed.list.FCS.files,
              label.key = label.key))
}

#ANPELA_MultiClusterFCS---------------------------------------------------
ANPELA_MultiClusterFCS <- function (list.of.files, clustering.var, numcluster, distance.metric,
                                    cluster.mode)
{
  list.of.FLOWMAP.clusters <- list()
  numcluster.orig <- numcluster
  if (length(numcluster.orig) == 1) {
    cat("Clustering all files to:", numcluster, "\n")
  }
  for (t in 1:length(list.of.files)) {
    cat("Clustering all files from time", t, "\n")
    fcs.files <- list.of.files[[t]]
    if (length(numcluster.orig) == 1) {
      numcluster.new <- rep(numcluster.orig, times = length(fcs.files))
      numcluster <- numcluster.new
    } else {
      numcluster <- numcluster.orig[[t]]
    }
    if (length(numcluster) != length(fcs.files)) {
      stop("Cluster number not specified for all FCS files!")
    }
    file.clusters <- ANPELA_ClusterFCS(fcs.files, clustering.var,
                                       numcluster, distance.metric, cluster.mode)
    list.of.FLOWMAP.clusters[[t]] <- file.clusters
  }
  return(list.of.FLOWMAP.clusters)
}
#ANPELA_ClusterFCS---------------------------------------------------
ANPELA_ClusterFCS <- function (fcs.files, clustering.var, numcluster, distance.metric,
                               cluster.mode)
{
  full.clusters <- data.frame()
  table.breaks <- c()
  table.lengths <- c()
  cluster.medians <- list()
  cluster.counts <- list()
  cell.assgn <- list()
  if (length(numcluster) == 1) {
    cat("Clustering all files to:", numcluster, "\n")
    numcluster.new <- rep(numcluster, times = length(fcs.files))
    numcluster <- numcluster.new
  }
  if (length(numcluster) != length(fcs.files)) {
    stop("Cluster number not specified for all FCS files!")
  }
  for (i in 1:length(fcs.files)) {
    current.file <- fcs.files[[i]]
    global.current.file <<- current.file
    global.clustering.var <<- clustering.var
    cat("Clustering data from file", i, "\n")
    cluster.counts[[i]] <- data.frame()
    cluster.medians[[i]] <- data.frame()
    tmp.FCS.for.cluster <- subset(x = current.file, select = clustering.var)
    cat("Subsetting for clustering channels only", "\n")
    cat("Using clustering method:", cluster.mode, "\n")
    if (cluster.mode == "hclust") {
      cluster.results <- FLOWMAPR:::HclustClustering(current.file = current.file,
                                                     tmp.FCS.for.cluster = tmp.FCS.for.cluster, distance.metric = distance.metric,
                                                     numcluster = numcluster[i])
    }
    else if (cluster.mode == "kmeans") {
      cluster.results <- FLOWMAPR:::KMeansClustering(current.file = current.file,
                                                     tmp.FCS.for.cluster = tmp.FCS.for.cluster, distance.metric = distance.metric,
                                                     numcluster = numcluster[i])
    }
    else if (cluster.mode == "none") {
      new.counts <- data.frame()
      cluster.results <- list(tmp.cell.assgn = data.frame(1:dim(tmp.FCS.for.cluster)[1]),
                              new.medians = tmp.FCS.for.cluster, new.counts = data.frame(rep.int(1,
                                                                                                 dim(tmp.FCS.for.cluster)[1])))
    } else {
      stop("Unrecognized clustering method!")
    }
    cell.assgn[[i]] <- cluster.results$tmp.cell.assgn
    cluster.medians[[i]] <- cluster.results$new.medians
    cluster.counts[[i]] <- cluster.results$new.counts
    colnames(cluster.counts[[i]]) <- c("Counts")
    table.lengths <- append(table.lengths, nrow(cluster.medians[[i]]))
    full.clusters <- rbind(full.clusters, cluster.medians[[i]])
    table.breaks <- append(table.breaks, nrow(full.clusters))
  }
  for (i in 1:length(cluster.medians)) {
    rownames(cluster.medians[[i]]) <- seq(1, table.lengths[i])
  }
  rownames(full.clusters) <- seq(1, dim(full.clusters)[1])
  FLOWMAP.clusters <- FLOWMAPR:::FLOWMAPcluster(full.clusters = full.clusters,
                                                table.breaks = table.breaks, table.lengths = table.lengths,
                                                cluster.medians = cluster.medians, cluster.counts = cluster.counts,
                                                cell.assgn = cell.assgn)
  return(FLOWMAP.clusters)
}

#ANPELA_BuildMultiFLOWMAPkNN---------------------------------------------------
ANPELA_BuildMultiFLOWMAPkNN <- function (remodel.FLOWMAP.clusters, k, min, max, distance.metric,
                                         label.key, clustering.var)
{
  global.remodel.FLOWMAP.clusters <<- remodel.FLOWMAP.clusters
  first.results <- FLOWMAPR:::BuildFirstMultiFLOWMAPkNN(remodel.FLOWMAP.clusters,
                                                        k, min, max, distance.metric = distance.metric, clustering.var)
  output.graph <- first.results$output.graph
  knn.indexes <- list()
  knn.distances <- list()
  knn.indexes[[1]] <- first.results$indexes
  knn.distances[[1]] <- first.results$distances
  table.lengths <- remodel.FLOWMAP.clusters$table.lengths
  table.breaks <- c(0, remodel.FLOWMAP.clusters$table.breaks)
  for (n in 1:(length(remodel.FLOWMAP.clusters$cluster.medians) -
               1)) {
    offset <- table.breaks[n] + 1
    cat(paste0("Starting next round, offset = ", as.character(offset),
               "\n"))
    cat(paste0("Building FLOWMAP from", n, "to", n + 1,
               "\n"))
    clusters <- rbind(remodel.FLOWMAP.clusters$cluster.medians[[n]],
                      remodel.FLOWMAP.clusters$cluster.medians[[n + 1]])
    clusters <- subset(clusters, select = clustering.var)
    results <- ANPELA_BaseBuildKNN(clusters = clusters, table.breaks = table.breaks,
                                   offset = offset, output.graph = output.graph, n = n,
                                   k = k, min = min, max = max, distance.metric = distance.metric)
    output.graph <- results$output.graph
    knn.indexes[[n]] <- data.frame(cbind(knn.indexes[[n]],
                                         results$indexes[1:table.lengths[n], ]))
    knn.distances[[n]] <- data.frame(cbind(knn.distances[[n]],
                                           results$distances[1:table.lengths[n], ]))
    knn.indexes[[n + 1]] <- results$indexes[(table.lengths[n] +
                                               1):nrow(results$indexes), ]
    knn.distances[[n + 1]] <- results$distances[(table.lengths[n] +
                                                   1):nrow(results$distances), ]
  }
  knn.indexes[[n + 1]] <- data.frame(knn.indexes[[n + 1]])
  knn.distances[[n + 1]] <- data.frame(knn.distances[[n +
                                                        1]])
  global.pre.simplify <<- knn.indexes
  keep_ind = lapply(knn.indexes, function(y) apply(y, 1, function(x) {
    to.keep <- c()
    for (i in unique(x)) {
      to.keep <- rbind(to.keep, head(which(x == i), 1))
    }
    to.keep
  }))
  global.keep.ind <<- keep_ind
  knn.ind.sim <- lapply(1:(length(keep_ind) - 1), function(x) {
    knn.indexes.keep <- c()
    for (i in 1:length(keep_ind[[x]])) {
      knn.indexes.keep <- dplyr::bind_rows(knn.indexes.keep,
                                           knn.indexes[[x]][i, keep_ind[[x]][[i]]])
    }
    print("kept")
    knn.indexes.keep
  })
  knn.ind.sim <- dplyr::bind_rows(knn.ind.sim, knn.indexes[[length(knn.indexes)]])
  knn.dist.sim <- lapply(1:(length(keep_ind) - 1), function(x) {
    knn.distances.keep <- c()
    for (i in 1:length(keep_ind[[x]])) {
      knn.distances.keep <- dplyr::bind_rows(knn.distances.keep,
                                             knn.distances[[x]][i, keep_ind[[x]][[i]]])
    }
    knn.distances.keep
  })
  knn.dist.sim <- dplyr::bind_rows(knn.dist.sim, knn.distances[[length(knn.distances)]])
  print("post-simplify")
  knn.out <- list(indexes = knn.ind.sim, distances = knn.dist.sim)
  global.knn.out.pre <<- knn.out
  print("processing knn output")
  for (i in 1:ncol(knn.out$indexes)) knn.out$indexes[sapply(knn.out$indexes[, i], is.null), i] <- 0
  for (i in 1:ncol(knn.out$distances)) knn.out$distances[sapply(knn.out$distances[, i], is.null), i] <- 10000
  knn.out$indexes <- as.data.frame(matrix(unlist(knn.out$indexes),
                                          nrow = length(unlist(knn.out$indexes[1]))))
  knn.out$distances <- as.data.frame(matrix(unlist(knn.out$distances),
                                            nrow = length(unlist(knn.out$distances[1]))))

  dist_order = t(apply(knn.out$distances, 1, order))
  for(i in 1:nrow(knn.out$indexes)) knn.out$indexes[i,] <- knn.out$indexes[i,dist_order[i,]] #order indexes by increasing distance (necessary to incorporate second timepoint connections into overall order)
  for(i in 1:nrow(knn.out$distances)) knn.out$distances[i,] <- knn.out$distances[i,][dist_order[i,]] #increasing order distances


  knn.out$indexes <- as.data.frame(matrix(unlist(knn.out$indexes),
                                          nrow = length(unlist(knn.out$indexes[1]))))
  knn.out$distances <- as.data.frame(matrix(unlist(knn.out$distances),
                                            nrow = length(unlist(knn.out$distances[1]))))
  print("finished processing knn output")
  global.knn.out <<- knn.out
  global.pre.output.graph <<- output.graph
  output.graph.final <- igraph::graph.empty()
  vertices.edges <- igraph::as_edgelist(output.graph)
  vertices.edges <- as.matrix(data.frame(lapply(data.frame(vertices.edges),
                                                function(x) as.numeric(as.character(x)))))
  output.graph.final <- igraph::graph_from_edgelist(vertices.edges,
                                                    directed = FALSE)
  igraph::E(output.graph.final)$weight <- unlist(lapply(igraph::E(output.graph)$weight,
                                                        function(x) 1/as.numeric(x)))
  output.graph.final <- ANPELA_AnnotateMultiGraph(output.graph, remodel.FLOWMAP.clusters,
                                                      label.key)
  output.graph.final <- igraph::simplify(output.graph.final)
  return(list(output.graph = output.graph.final, knn.out = knn.out))
}
#ANPELA_BaseBuildKNN---------------------------------------------------
ANPELA_BaseBuildKNN <- function (clusters, table.breaks, offset, n, k, min, max, output.graph,
                                 edgelist.save, distance.metric)
{
  cat(paste0("Clusters length:", as.character(dim(clusters)[1]),
             "\n"))
  if (k < max) {
    cat(paste0("Input k (", k, ") less than max edge number (",
               max, ") -- changing k to max edge number"))
    k <- max
  }
  if (distance.metric == "manhattan") {
    cat("Manhattan distance no longer supported. Using euclinean distance.\n")
    nns <- RANN::nn2(data = clusters, k = k + 1, searchtype = "priority",
                     eps = 0.1)
  }
  else if (distance.metric == "euclidean") {
    nns <- RANN::nn2(data = clusters, k = k + 1, searchtype = "priority",
                     eps = 0.1)
  }
  temp_nnids.df <- as.data.frame(nns$nn.idx)
  temp_nndists.df <- as.data.frame(nns$nn.dists)
  nn.ids.df <- temp_nnids.df[, 2:length(temp_nnids.df)]
  nn.dists.df <- temp_nndists.df[, 2:length(temp_nndists.df)]
  if (offset > 0) {
    cat(paste0("In BaseBuildKNN, offset = ", as.character(offset),
               "\n"))
    cat(paste0("In BaseBuildKNN, table.breaks[n + 2] = ",
               as.character(table.breaks[n + 2]), "\n"))
    row.names(nn.ids.df) <- offset:table.breaks[n + 2]
    nn.ids.df[] <- lapply(nn.ids.df, function(x) x + table.breaks[n])
    row.names(nn.dists.df) <- offset:table.breaks[n + 2]
  }
  numcluster <- nrow(clusters)
  normalized.densities <- FLOWMAPR:::KnnDensity(k = k, min, max, n = n,
                                                nn.ids.df = nn.ids.df, nn.dists.df = nn.dists.df, numcluster = numcluster,
                                                table.breaks = table.breaks, offset = offset)
  results <- FLOWMAPR:::DrawNormalizedEdgesKnn(output.graph = output.graph,
                                               nn.ids.df = nn.ids.df, nn.dists.df = nn.dists.df, normalized.densities = normalized.densities,
                                               n = n, table.breaks = table.breaks, offset = offset)
  output.graph <- results$output.graph
  edgelist.save <- results$edgelist.with.distances
  knn.indexes <- results$indexes
  knn.distances <- results$distances
  if (!igraph::is_connected(output.graph)) {
    print("Graph has disconnected components!")
    if (offset == 0) {
      results <- FLOWMAPR:::FirstConnectSubgraphs(output.graph = output.graph,
                                                  edge.list = edgelist.save, offset = offset,
                                                  table.breaks = table.breaks, n = n, distance.metric = distance.metric,
                                                  clusters = clusters)
    }
    else if (offset > 0) {
      results <- ANPELA_ConnectSubgraphs(output.graph = output.graph,
                                         edge.list = edgelist.save, offset = offset,
                                         table.breaks = table.breaks, n = n, distance.metric = distance.metric,
                                         clusters = clusters)
    }
    output.graph <- results$output.graph
    edgelist.save <- results$edgelist.with.distances
  }
  return(list(output.graph = output.graph, indexes = knn.indexes,
              distances = knn.distances))
}
#ANPELA_ConnectSubgraphs---------------------------------------------------
ANPELA_ConnectSubgraphs <- function (output.graph, edge.list, offset, table.breaks, n,
                                     distance.metric, clusters)
{
  print("in connected")
  orig.nrow.edge.list <- nrow(edge.list)
  time.prox.graph <- igraph::graph.empty()
  edge.list$id_num <- edge.list$id - offset + 1
  edge.list$index_num <- edge.list$index - offset + 1
  time.prox.graph <- igraph::graph_from_edgelist(as.matrix(cbind(edge.list$id_num,
                                                         edge.list$index_num)), directed = FALSE)
  igraph::E(time.prox.graph)$weight <- edge.list$dist
  time.prox.graph <- igraph::set_vertex_attr(time.prox.graph,
                                     "name", index = igraph::V(time.prox.graph), as.character(offset:table.breaks[n +
                                                                                                            2]))
  subgraphs.ls <- igraph::decompose(time.prox.graph)
  y <- 0
  to.add.df <- NA
  while (length(subgraphs.ls) >= 2) {
    y <- y + 1
    print(y)
    if (length(igraph::as_edgelist(subgraphs.ls[[2]])) != 0) {
      print("in loop")
      rownames(clusters) <- c(offset:table.breaks[n +
                                                    2])
      i <- 1
      j <- 1
      to.add.df <- data.frame()
      for (p in 1:length(subgraphs.ls)) {
        si.graph <- subgraphs.ls[[p]]
        print(p)
        inter.sub.nns.df <- data.frame()
        for (m in 1:length(subgraphs.ls)) {
          if (m > length(subgraphs.ls)) {
            break
          }
          if (m == p) {
            next
          }
          sj.graph <- subgraphs.ls[[m]]
          if (distance.metric == "manhattan") {
            print("Manhattan distance no longer supported. Using euclinean distance.")
            inter.sub.nns <- RANN::nn2(data = clusters[igraph::V(sj.graph)$name,
            ], query = clusters[igraph::V(si.graph)$name,
            ], k = 1, searchtype = "priority", eps = 0.1)
          }
          else if (distance.metric == "euclidean") {
            inter.sub.nns <- RANN::nn2(data = clusters[igraph::V(sj.graph)$name,
            ], query = clusters[igraph::V(si.graph)$name,
            ], k = 1, searchtype = "priority", eps = 0.1)
          }
          temp.inter.sub.nns.df <- data.frame(id = igraph::V(si.graph)$name)
          temp.inter.sub.nns.df <- cbind(temp.inter.sub.nns.df,
                                         as.data.frame(igraph::V(sj.graph)$name[as.numeric(inter.sub.nns$nn.idx)]),
                                         as.data.frame(inter.sub.nns$nn.dists))
          colnames(temp.inter.sub.nns.df) <- c("id",
                                               "index", "dist")
          inter.sub.nns.df <- rbind(inter.sub.nns.df,
                                    temp.inter.sub.nns.df)
        }
        inter.sub.nns.df <- inter.sub.nns.df[order(inter.sub.nns.df$dist),
        ]
        to.add.df <- rbind(to.add.df, inter.sub.nns.df[1,
        ])
      }
      to.add.df <- data.frame(lapply(data.frame(to.add.df),
                                     function(x) as.numeric(as.character(x))))
      to.add.df$id_num <- to.add.df$id - offset + 1
      to.add.df$index_num <- to.add.df$index - offset +
        1
      edge.list <- rbind(edge.list, to.add.df)
      output.graph.update <- igraph::graph.empty()
      output.graph.update <- igraph::graph_from_edgelist(as.matrix(edge.list[,
                                                                     c("id_num", "index_num")]), directed = FALSE)
      igraph::E(output.graph.update)$weight <- edge.list$dist
      output.graph.update <- igraph::set.vertex.attribute(output.graph.update,
                                                  "name", index = igraph::V(output.graph.update), as.character(offset:table.breaks[n +
                                                                                                                             2]))
      subgraphs.ls <- igraph::decompose.graph(output.graph.update)
      subgraphs.el.ls <- list()
      for (ind in 1:length(subgraphs.ls)) {
        subgraphs.el.ls[[ind]] <- igraph::as_edgelist(subgraphs.ls[[ind]])
      }
    }
    else {
      break
    }
  }
  if (any(!is.na(to.add.df))) {
    output.graph.update <- igraph::graph.empty()
    og.el <- igraph::as_edgelist(output.graph)
    colnames(og.el) <- c("id", "index")
    og.el <- rbind(og.el[, 1:2], edge.list[(orig.nrow.edge.list +
                                              1):nrow(edge.list), c("id", "index")])
    output.graph.update <- igraph::graph_from_edgelist(as.matrix(og.el),
                                               directed = FALSE)
    ogu.weights <- c(igraph::E(output.graph)$weight, edge.list[(orig.nrow.edge.list +
                                                          1):nrow(edge.list), "dist"])
    igraph::E(output.graph.update)$weight <- ogu.weights
    igraph::V(output.graph.update)$name <- 1:length(igraph::V(output.graph.update))
    output.graph <- output.graph.update
  }
  edgelist.with.distances <- edge.list[, c("id", "index",
                                           "dist")]
  colnames(edgelist.with.distances) <- c("row.inds", "col.inds",
                                         "values")
  results <- list(output.graph = output.graph, edgelist.with.distances = edge.list)
  return(results)
}

#ANPELA_RunUMAPlayout---------------------------------------------------
ANPELA_RunUMAPlayout <- function (graph, knn.in, file.clusters, clustering.var, file.name = file.name,
                                  umap_n_neighbors, k, umap_n_components, mode, umap_min_dist) #umap_min_dist
{
  cat("running UMAP\n")
  self_col_ind <- c(1:nrow(knn.in$indexes))
  self_col_dist <- rep(0, nrow(knn.in$indexes))
  knn.in$indexes <- data.frame(cbind(self_col_ind, data.frame(knn.in$indexes)))
  knn.in$indexes <- as.matrix(knn.in$indexes)
  knn.in$distances <- data.frame(cbind(self_col_dist, knn.in$distances))
  knn.in$distances <- as.matrix(knn.in$distances)
  knn <- list()
  if (umap_n_neighbors > k) {
    umap_n_neighbors <- k
    print("umap_n_neighboors must be <= k, has been changed to k")
  }
  umap_n_neighbors <- umap_n_neighbors + 1
  knn[["idx"]] <- as.matrix(knn.in$indexes[, 1:umap_n_neighbors])
  knn[["dist"]] <- as.matrix(knn.in$distances[, 1:umap_n_neighbors])
  umap.out <- uwot::umap(file.clusters$full.clusters, ret_nn = TRUE,
                         nn_method = knn, verbose = TRUE, n_components = umap_n_components,
                         min_dist=umap_min_dist, spread=1)

  all.attr <- igraph::get.vertex.attribute(graph)
  this.df <- c()
  for (x in 1:length(all.attr)) {
    this.df <- cbind(this.df, all.attr[[x]])
    this.df <- as.data.frame(this.df)
  }
  colnames(this.df) <- names(all.attr)
  colnames(umap.out$embedding) <- c("x", "y")
  this.df <- cbind(this.df, umap.out$embedding)
  return(this.df)

}

#ANPELA_layout_forceatlas2---------------------------------------------------
Rcpp::sourceCpp("./src/forceatlas2.cpp") # load ANPELA_layout_forceatlas2Cpp
ANPELA_layout_forceatlas2 <- function (G, ew.influence = 1, kgrav = 1, iter = 1000, prevent.overlap = FALSE,
                                       fixed = NULL, stopping.tolerance = 0.001, barnes.hut = FALSE)
{
  v.count <- igraph::vcount(G)
  if (v.count >= 2000) {
    barnes.hut <- TRUE
  }
  if (v.count > 2000) {
    stopping.tolerance <- 0.01
  } else if (v.count > 800) {
    stopping.tolerance <- 0.005
  } else stopping.tolerance <- 0.001

  if (is.null(fixed)) {
    fixed <- rep(FALSE, v.count)
  }

  lay <- NULL

  if (is.null(igraph::get.vertex.attribute(G, "x"))) {
    lay <- matrix(ncol = 2, nrow = v.count, data = rnorm(v.count *
                                                           2, 10, 2))
    colnames(lay) <- c("x", "y")
  } else {
    lay <- cbind(x = igraph::V(G)$x, y = igraph::V(G)$y)
    w <- is.na(lay[, "x"])
    if (any(w))
      lay[w, ] <- rnorm(sum(w) * 2, 10, 2)
  }
  if (is.null(igraph::get.vertex.attribute(G, "size"))) {
    igraph::V(G)$size <- rep(10, v.count)
  }
  mass <- 1 + igraph::degree(G)
  F_att <- (igraph::E(G)$weight^ew.influence)
  edge_list <- igraph::as_edgelist(G, names = F) - 1
  avg_displ <- numeric(iter)
  max_displ <- numeric(iter)
  if (barnes.hut) {
    message("Using Barnes-Hut approximation\n")
  }
  message(sprintf("Stopping tolerance: %f\n", stopping.tolerance))
  flush.console()
  ANPELA_layout_forceatlas2Cpp(lay, F_att, mass, igraph::V(G)$size, edge_list,
                               avg_displ, kgrav, iter, prevent.overlap, fixed, max_displ,
                               stopping.tolerance, barnes.hut)
  return(list(lay = lay, avg_displ = avg_displ, max_displ = max_displ))
}

#ANPELA_RemodelFLOWMAPClusterList---------------------------------------------------
ANPELA_RemodelFLOWMAPClusterList <- function (list.of.FLOWMAP.clusters, label.key)
{
  full.clusters <- data.frame()
  table.breaks <- c()
  table.lengths <- c()
  cluster.medians <- list()
  cluster.counts <- list()
  cell.assgn <- list()
  for (t in 1:length(list.of.FLOWMAP.clusters)) {
    temp.medians <- data.frame()
    temp.cell.assgn <- data.frame()
    temp.counts <- data.frame()
    for (c in 1:length(list.of.FLOWMAP.clusters[[t]]$cluster.medians)) {
      temp.medians <- rbind(temp.medians, list.of.FLOWMAP.clusters[[t]]$cluster.medians[[c]])
      temp.cell.assgn <- rbind(temp.cell.assgn, list.of.FLOWMAP.clusters[[t]]$cell.assgn[[c]])
      temp.counts <- c(unlist(temp.counts), unlist(list.of.FLOWMAP.clusters[[t]]$cluster.counts[[c]]))
      temp.counts <- as.data.frame(temp.counts)
      colnames(temp.counts) <- c("Counts")
    }
    cluster.medians[[t]] <- temp.medians
    rownames(temp.counts) <- 1:nrow(temp.counts)
    cluster.counts[[t]] <- temp.counts
    cell.assgn[[t]] <- temp.cell.assgn
    table.lengths <- c(table.lengths, dim(temp.medians)[1])
    table.breaks <- c(table.breaks, sum(table.lengths))
    full.clusters <- rbind(full.clusters, temp.medians)
  }
  full.clusters <- ANPELA_ConvertCharacterLabel(full.clusters, label.key)
  remodeled.FLOWMAP.clusters <- FLOWMAPR:::FLOWMAPcluster(full.clusters,
                                               table.breaks, table.lengths, cluster.medians, cluster.counts,
                                               cell.assgn)
  return(remodeled.FLOWMAP.clusters)
}

#ANPELA_ConvertCharacterLabel---------------------------------------------------
ANPELA_ConvertCharacterLabel <- function (data.frame.with.numeric.labels, label.key)
{
  data.frame.with.character.labels <- data.frame.with.numeric.labels
  times <- unique(data.frame.with.numeric.labels[, "Time"])
  for (t in 1:length(times)) {
    this.label <- label.key[[t]]
    this.ind <- which(data.frame.with.numeric.labels[, "Time"] ==
                        times[t])
    for (i in 1:length(this.label)) {
      fix.ind <- which(data.frame.with.numeric.labels[,
                                                      "Condition"] == i)
      use.ind <- intersect(fix.ind, this.ind)
      data.frame.with.character.labels[use.ind, "Condition"] <- this.label[i]
    }
  }
  return(data.frame.with.character.labels)
}

#ANPELA_AnnotateMultiGraph------------------------------------------------------------
ANPELA_AnnotateMultiGraph <- function (output.graph, list.of.FLOWMAP.clusters, label.key)
{
  anno <- list()
  times <- 1:length(list.of.FLOWMAP.clusters$cluster.medians)
  anno$medians <- data.frame()
  anno$count <- data.frame()
  anno$percent.total <- data.frame()
  for (t in times) {
    cat("Annotating graph for file", t, "\n")
    anno$medians <- rbind(anno$medians, list.of.FLOWMAP.clusters$cluster.medians[[t]])
    anno$medians[, "Condition"]
    anno$count <- rbind(anno$count, list.of.FLOWMAP.clusters$cluster.counts[[t]])
    total.cell <- sum(list.of.FLOWMAP.clusters$cluster.counts[[t]])
    percent.total <- data.frame(list.of.FLOWMAP.clusters$cluster.counts[[t]]/total.cell)
    colnames(percent.total) <- "percent.total"
    anno$percent.total <- rbind(anno$percent.total, percent.total)
  }
  output.anno <- cbind(anno$medians, anno$count, anno$percent.total)
  global.output.anno.pre <<- output.anno
  print("label.key")
  print(label.key)
  output.anno <- ANPELA_ConvertCharacterLabel(output.anno, label.key)
  global.output.anno.post <<- output.anno
  for (c in colnames(output.anno)) {
    output.graph <- igraph::set.vertex.attribute(output.graph, c,
                                         index = as.numeric(1:dim(output.anno)[1]), value = output.anno[,
                                                                                                        c])
  }
  igraph::V(output.graph)$name <- 1:length(igraph::V(output.graph))
  global.output.graph <<- output.graph
  return(output.graph)
}

