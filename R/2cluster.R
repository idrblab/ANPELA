data_cluster <- function (data = NULL, method = c("FlowSOM", "PhenoGraph","Mclust"), 
                          Phenograph_k = 30, FlowSOM_k = 20, FlowSeed = NULL) {
  method = match.arg(method)
  switch(method, PhenoGraph = {
    set.seed(123)
    clusters <- as.numeric(Rphenograph::Rphenograph(data, k = Phenograph_k)[[2]]$membership)
  }, FlowSOM = {
    set.seed(FlowSeed)
    clusters <- FlowSOM_integrate2cytofkit2(xdata = data, k = FlowSOM_k, FlowSeed = FlowSeed)
  },Mclust = { clusters <- as.numeric(mclust::Mclust(data)$classification)})
  if ((length(clusters) != length(data)) && (length(clusters) != nrow(data))) {
    # message("Cluster is not complete, cluster failed, try another cluster method!")
    # return(NULL)
    return("Cluster is not complete, cluster failed, try another cluster method!")
  } else {
    names(clusters) <- row.names(data)
    return(clusters)
  }
}


FlowSOM_integrate2cytofkit2 <- function(xdata, k, FlowSeed = NULL, ...){
  xdata <- as.matrix(xdata)
  
  ord <- tryCatch({
    map <- FlowSOM::SOM(xdata, silent = TRUE, ...)
    metaClusters <- suppressMessages(metaClustering_consensus2(map$codes, k = k, seed = FlowSeed))
    cluster <- metaClusters[map$mapping[,1]]
  }, error=function(cond) {
    message("Run FlowSOM failed")
    message("Here's the error message:")
    # message(cond)
    return(NULL)
  }) 
  
  if(is.null(ord)){
    cluster <- NULL
  }else{
    if(length(ord) != nrow(xdata)){
      message("Run FlowSOM failed!")
      return(NULL)
    }
    cluster <- ord
  }
  
  return(cluster)
}

metaClustering_consensus2 <- function (data, k = 7, seed = NULL) {
  results <- suppressMessages(ConsensusClusterPlus2(t(data), 
                                                    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(), 
                                                    plot = "pdf", verbose = FALSE, clusterAlg = "hc", 
                                                    distance = "euclidean", seed = seed))
  results[[k]]$consensusClass
}

ConsensusClusterPlus2 <- function (d = NULL, maxK = 3, reps = 10, pItem = 0.8, pFeature = 1, 
                                   clusterAlg = "hc", title = "untitled_consensus_cluster", 
                                   innerLinkage = "average", finalLinkage = "average", 
                                   distance = "pearson", ml = NULL, tmyPal = NULL, seed = NULL, 
                                   plot = NULL, writeTable = FALSE, weightsItem = NULL, weightsFeature = NULL, 
                                   verbose = F, corUse = "everything") {
  if (is.null(seed) == TRUE) {
    seed = timeSeed = as.numeric(Sys.time())
  }
  set.seed(seed)
  if (is.null(ml) == TRUE) {
    if (all(!class(d) %in% c("dist", "matrix", "ExpressionSet"))) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }
    if (inherits(d, "dist")) {
      if (is.null(attr(d, "method"))) {
        attr(d, "method") <- distance <- "unknown - user-specified"
      }
      if (is.null(distance) || (distance != attr(d, "method"))) {
        distance <- attr(d, "method")
      }
      if ((!is.null(pFeature)) && (pFeature < 1)) {
        message("Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n")
        pFeature <- 1
      }
      if (!is.null(weightsFeature)) {
        message("Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n")
        weightsFeature <- NULL
      }
      if (clusterAlg == "km") {
        message("Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix")
      }
    }
    else {
      if (is.null(distance)) {
        distance <- "pearson"
      }
    }
    if ((clusterAlg == "km") && inherits(distance, 
                                         "character") && (distance != "euclidean")) {
      message("Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the 'kmdist' option, or use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
      distance <- "euclidean"
    }
    if (inherits(d, "ExpressionSet")) {
      d <- exprs(d)
    }
    ml <- ccRun2(d = d, maxK = maxK, repCount = reps, diss = inherits(d, 
                                                                      "dist"), pItem = pItem, pFeature = pFeature, 
                 innerLinkage = innerLinkage, clusterAlg = clusterAlg, 
                 weightsFeature = weightsFeature, weightsItem = weightsItem, 
                 distance = distance, verbose = verbose, corUse = corUse)
  }
  res = list()
  if ((is.null(plot) == FALSE | writeTable) & !file.exists(paste(title, 
                                                                 sep = ""))) {
    dir.create(paste(title, sep = ""))
  }
  log <- matrix(ncol = 2, byrow = T, c("title", title, 
                                       "maxK", maxK, "input matrix rows", ifelse(inherits(d, 
                                                                                          "matrix"), nrow(d), "dist-mat"), "input matrix columns", 
                                       ifelse(inherits(d, "matrix"), ncol(d), ncol(as.matrix(d))), 
                                       "number of bootstraps", reps, "item subsampling proportion", 
                                       pItem, "feature subsampling proportion", ifelse(is.null(pFeature), 
                                                                                       1, pFeature), "cluster algorithm", clusterAlg, 
                                       "inner linkage type", innerLinkage, "final linkage type", 
                                       finalLinkage, "correlation method", distance, "plot", 
                                       if (is.null(plot)) NA else plot, "seed", if (is.null(seed)) NA else seed))
  colnames(log) = c("argument", "value")
  if (writeTable) {
    write.csv(file = paste(title, "/", title, ".log.csv", 
                           sep = ""), log, row.names = F)
  }
  if (is.null(plot)) {
  }
  else if (plot == "pngBMP") {
    bitmap(paste(title, "/", "consensus%03d.png", 
                 sep = ""))
  }
  else if (plot == "png") {
    png(paste(title, "/", "consensus%03d.png", 
              sep = ""))
  }
  else if (plot == "pdf") {
    pdf(onefile = TRUE, paste(title, "/", "consensus.pdf", 
                              sep = ""))
  }
  else if (plot == "ps") {
    postscript(onefile = TRUE, paste(title, "/", "consensus.ps", 
                                     sep = ""))
  }
  colorList = list()
  colorM = rbind()
  thisPal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", 
               "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", 
               "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", 
               "#B15928", "#bd18ea", "#2ef4ca", "#f4cced", 
               "#f4cc03", "#05188a", "#e5a25a", "#06f106", 
               "#85848f", "#000000", "#076f25", "#93cd7f", 
               "#4d0776", "#ffffff")
  colBreaks = NA
  if (is.null(tmyPal) == TRUE) {
    colBreaks = 10
    tmyPal = myPal2(colBreaks)
  }
  else {
    colBreaks = length(tmyPal)
  }
  sc = cbind(seq(0, 1, by = 1/(colBreaks)))
  rownames(sc) = sc[, 1]
  sc = cbind(sc, sc)
  heatmap(sc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", 
          col = tmyPal, na.rm = TRUE, labRow = rownames(sc), labCol = F, 
          main = "consensus matrix legend")
  for (tk in 2:maxK) {
    if (verbose) {
      message(paste("consensus ", tk))
    }
    fm = ml[[tk]]
    hc = hclust(as.dist(1 - fm), method = finalLinkage)
    message("clustered")
    ct = cutree(hc, tk)
    names(ct) = colnames(d)
    if (all(class(d) == "dist")) {
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    colorList = setClusterColors2(res[[tk - 1]][[3]], ct, 
                                  thisPal, colorList)
    pc = c
    pc = pc[hc$order, ]
    if (!is.null(plot) && plot == "pngBMP") {
      pc = pc[, hc$order]
      pc = rbind(pc, 0)
      oc = colorList[[1]][hc$order]
      heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", 
              col = tmyPal, na.rm = TRUE, labRow = F, labCol = F, 
              mar = c(5, 5), main = paste("consensus matrix k=", 
                                          tk, sep = ""), ColSideCol = oc)
    }
    else {
      pc = rbind(pc, 0)
      heatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, 
              symm = FALSE, scale = "none", col = tmyPal, 
              na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 
                                                            5), main = paste("consensus matrix k=", 
                                                                             tk, sep = ""), ColSideCol = colorList[[1]])
    }
    legend("topright", legend = unique(ct), fill = unique(colorList[[1]]), 
           horiz = FALSE)
    res[[tk]] = list(consensusMatrix = c, consensusTree = hc, 
                     consensusClass = ct, ml = ml[[tk]], clrs = colorList)
    colorM = rbind(colorM, colorList[[1]])
  }
  CDF2(ml)
  # clusterTrackingPlot(colorM[, res[[length(res)]]$consensusTree$order])
  if (is.null(plot) == FALSE) {
    dev.off()
  }
  res[[1]] = colorM
  if (writeTable) {
    for (i in 2:length(res)) {
      write.csv(file = paste(title, "/", title, ".k=", 
                             i, ".consensusMatrix.csv", sep = ""), 
                res[[i]]$consensusMatrix)
      write.table(file = paste(title, "/", title, 
                               ".k=", i, ".consensusClass.csv", 
                               sep = ""), res[[i]]$consensusClass, col.names = F, 
                  sep = ",")
    }
  }
  return(res)
}

ccRun2 <- function( d=d,
                    maxK=NULL,
                    repCount=NULL,
                    diss=inherits( d, "dist" ),
                    pItem=NULL,
                    pFeature=NULL,
                    innerLinkage=NULL,
                    distance=NULL, #ifelse( inherits(d,"dist"), attr( d, "method" ), "euclidean" ),@@@@@
                    clusterAlg=NULL,
                    weightsItem=NULL,
                    weightsFeature=NULL,
                    verbose=NULL,
                    corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);
  
  if (is.null( distance ) ) distance <- 'euclidean'  ## necessary if d is a dist object and attr( d, "method" ) == NULL
  
  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                            "pearson", "spearman" )
  
  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
         ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
         ( is.null( pFeature ) ||
           ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ) stop("unsupported distance.")
        
        if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        # }else if( class(try(get(distance),silent=T))=="function"){
          # main.dist.obj <- get(distance)( t( d )   )
        }else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance  
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }
  
  
  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols2( d, pItem, pFeature, weightsItem, weightsFeature )
    
    this_dist = NA
    if ( ! is.null( main.dist.obj ) ) {
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        this_dist <- as.dist( this_dist )
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {
      ## if main.dist.obj is NULL, then d is a data matrix, and either:
      ##   1) clusterAlg is 'km'
      ##   2) pFeatures < 1 or weightsFeatures have been specified, or
      ##   3) both
      ## so we can't use a main distance object and for every iteration, we will have to re-calculate either
      ##   1) the distance matrix (because we're also sampling the features as well), or
      ##   2) the submat (if using km) 
      
      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=T))!="function")  ) stop("unsupported distance.")
        if( ( class(try(get(distance),silent=T))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
        }else{
          if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,method=distance) )
          }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
          }
        }
        attr( this_dist, "method" ) <- distance  
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
             ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }
          
          this_dist <- sample_x$submat
        } 
      }
    }
    
    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix2( rep( 1,length(sample_x[[3]])),
                                   mCount,
                                   sample_x[[3]] ) 
    
    ##use samples for each k		
    for (k in 2:maxK){
      if(verbose){
        message(paste("  k =",k))
      }
      if (i==1){
        ml[[k]] = mConsist #initialize
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        ##prune to k for hc
        this_assignment = cutree(this_cluster,k)
        
        
      }else if(clusterAlg=="km"){
        ##this_dist should now be a matrix corresponding to the result from sampleCols
        this_assignment <- kmeans( t( this_dist ),
                                   k,
                                   iter.max = 10,
                                   nstart = 1,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- pam( x=this_dist,
                                k,
                                diss=TRUE,
                                metric=distance, 
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }
      ##add to tally				
      ml[[k]] <- connectivityMatrix2( this_assignment,
                                      ml[[k]],
                                      sample_x[[3]] )
    }
  }
  
  
  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    ##fill in other half of matrix for tally and count.
    tmp = triangle2(ml[[k]],mode=3)
    tmpCount = triangle2(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  message("end fraction")
  return(res)
}

sampleCols2 <- function( d,
                         pSamp=NULL,
                         pRow=NULL,
                         weightsItem=NULL,
                         weightsFeature=NULL ){
  ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
  ##  if no sampling over the features is performed, the submatrix & sample features are returned as NAs
  ##  to reduce memory overhead
  
  
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )
  
  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
         ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}

connectivityMatrix2 <- function( clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey 
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId
  
  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}

triangle2 <- function(m,mode=1){
  #mode=1 for CDF2, vector of lower triangle2.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary 
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m
  
  
  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  
  fm = t(nm)+nm
  diag(fm) = diag(m)
  
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  
  if(mode==1){
    return(vm) #vector 		
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle2 and no diagonal. no double counts.
  }
  
}

myPal2 <- function(n=10){
  #returns n colors
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
}

setClusterColors2 <-  function(past_ct,ct,colorU,colorList){
  #description: sets common color of clusters between different K
  newColors = c()
  if(length(colorList)==0){
    #k==2
    newColors = colorU[ct]
    colori=2
  }else{
    newColors = rep(NULL,length(ct))
    colori = colorList[[2]]
    mo=table(past_ct,ct)
    m=mo/apply(mo,1,sum)
    for(tci in 1:ncol(m)){ # for each cluster
      maxC = max(m[,tci])
      pci = which(m[,tci] == maxC)				
      if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
        #if new column maximum is unique, same cell is row maximum and is also unique
        ##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
        newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
      }else{ #add new color.
        colori=colori+1
        newColors[which(ct==tci)] = colorU[colori]
      }
    }
  }
  return(list(newColors,colori,unique(newColors) ))
}

CDF2 <- function(ml,breaks=100){
  #plot CDF2 distribution
  plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="consensus index",ylab="CDF2",main="consensus CDF2", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle2(ml[[i]],mode=1)
    
    #empirical CDF2 distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    #calculate area under CDF2 curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }
  legend(0.8,0.5,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),fill=this_colors)
  
  #plot area under CDF2 change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF2 curve",main="Delta area",type="b")
}

