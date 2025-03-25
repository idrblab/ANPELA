##
## change the starting point of the inferred pseudotemporal ordering
## (we assume that the experimental time is available)
##
## INPUT: - exp.time:       the initial experimental time vector e.g. exp.time = c(0,0,0,...,1,1,1,...)
##        - pseudotime:     a vector with the inferred pseudotimes
##        - ps.ordered:     flag, whether the exp.time and pseudotime are already in ascending order
##        - circular:       flag, whether the trajectory inference method used for extracting the pseudotime outputs circular topologies
##        - randomize_t0:   flag, whether to shift randomly according to the location of a cell @ t = 0
##
## OUTPUT: - the new pseudotime
##
## plot the location of the turnpoints
## plot(trnp$points)
## abline(v = trnp$tppos, col = "red" )

shift_start <- function(exp.time, pseudotime, ps.ordered = T, circular = F, randomize_t0 = T){
  
  # sort pseudotimes in ascending order
  if (!ps.ordered){
    cell_position =  order(pseudotime)
    pseudotime = pseudotime[cell_position]
    pseudo.ordered.exp.time = exp.time[cell_position]
    
    # find the indexes of cells of t = 0
    t0 = which(pseudo.ordered.exp.time == 0)
  } else {
    # find the indexes of cells of t = 0
    # t0 = which(exp.time == 0)
    t0 = which(exp.time == min(exp.time))
  }
  
  # find the turning points of the density of cells @ t = 0 on the pseudo-ordered vector
  trnp = pastecs::turnpoints(density(t0,kernel = "gaussian")$y)
  
  # if the trajectory inference method does not output a circular topology
  if (!circular){
    l = length(trnp$points)
    middle.point = round(l/2)
    area.left = sum(trnp$points[1:middle.point])
    area.right = sum(trnp$points[middle.point:l])
    if (area.right > area.left){
      new.pseudotime = reverse_pseudotime(pseudotime)
    } else { 
      new.pseudotime = pseudotime
    }
  } else {
    # if we want the new starting point to be where the density of cells from exp.time = 0 starts increasing
    if (randomize_t0) {
      # if the first turnpoint is a peak and this peak is the maximum one, don't shift
      if ((which.max(trnp$info) == 1) & (trnp$firstispeak)) {
        new.pseudotime = pseudotime
      } 
      # if the last turnpoint is a peak and this peak is the maximum one, reverse
      else if (which.max(trnp$info) == trnp$nturns){
        new.pseudotime = reverse_pseudotime(pseudotime)
      }
      # else, find the position of the pit with the maximum information
      else {
        sorted.trnp.info = order(trnp$info,decreasing = T)
        i <- 1
        is.pit = F
        while(is.pit != T){
          test.pit = trnp$pits[trnp$tppos]
          if (test.pit[sorted.trnp.info[i]]) { is.pit = T}
          else { i <- i + 1 }
        }
        max.info.pit.ix <- sorted.trnp.info[i]
        
        # recover the binning of the time vector
        bins = round(seq(1,length(exp.time),length(exp.time)/trnp$n))
    
        # the optimal starting point should be between the pit and the next best peak
        pit.pos = bins[trnp$tppos[max.info.pit.ix]]
        
        # if the next best peak is on the right
        if (trnp$info[max.info.pit.ix-1] < trnp$info[max.info.pit.ix+1]){
          peak.pos = bins[trnp$tppos[max.info.pit.ix+1]]  
          pit.peak.points = t0[ (t0 >= pit.pos) & (t0 <= peak.pos)]
          
          # set the starting point at 5% of the quantile range where t = 0 after the pit
          t0.init = t0[which.max(t0 > quantile(pit.peak.points, probs = .05))]
          
          # cycle shift pseudotime according to the inferred t0.init
          new.pseudotime = cycle.pseudotime(pseudotime,t0.init)
        } else 
          # if the next best peak is on the left we must both reverse and shift
          {
          peak.pos = bins[trnp$tppos[max.info.pit.ix-1]]  
          pit.peak.points = t0[ (t0 <= pit.pos) & (t0 >= peak.pos)]
          
          # set the starting point at 95% of the quantile range where t = 0 after the peak
          t0.init = t0[which.max(t0 > quantile(pit.peak.points, probs = .95))]
          
          # cycle shift pseudotime according to the inferred t0.init
          new.pseudotime = cycle.pseudotime(reverse_pseudotime(pseudotime),t0.init)
          
          # because the new.pseudotime had to be reversed its resulting range is [-1,0]
          new.pseudotime = new.pseudotime + 1
        }
      
        
      }
    } else {
      # # if we want the new starting point to be chosen randomly
      t0.init = 1
      totally_random == 1
      while (t0.init == 1) {
        if (totally_random == 1){ t0.init = sample(seq(1:length(exp.time)),1) }
        else { t0.init = t0[sample(seq(1:length(t0)),1)] }
      }
      
      # cycle shift pseudotime according to the inferred t0.init
      new.pseudotime = cycle.pseudotime(pseudotime,t0.init)
    }
  }
  return(new.pseudotime)
}