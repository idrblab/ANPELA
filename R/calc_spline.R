calc_spline <- function(dat,tm){
  dfo <- 3
  knts <- quantile(tm, probs = seq(0, 1, length = (dfo + 1)))[-c(1, (dfo + 1))]
  
  # remove zeros and ones , but keep only one of each
  if (length(which(knts == 0))>1 | length(which(knts == 1))>1 ) {
    knts <-  knts[-c(1:(length(which(knts == 0))-1) , (length(knts)-length(which(knts==1))+1) :(length(knts)-1))] 
  }
  names(dat) <- make.names(names(dat))
  # build the model
  fit1 <- lm(formula(paste("cbind(",paste(names(dat), collapse = ","),") ~ splines::ns(tm, knots = knts)")), data = dat, drop=TRUE)
  
  # the time points to predict
  pred.time <- seq(min(tm), max(tm), length.out = nrow(dat))
  
  # data.frame with the results of the fit
  df <- data.frame(time = pred.time, predict(fit1, data.frame(tm = pred.time)))
  
  return(df)
}