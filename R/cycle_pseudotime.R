cycle.pseudotime <- function(pst,new.init){
  # calculate the difference between consequtive pseudotime values
  d = diff(pst)
  
  # split the differences. tail(t): the ones from t0.init and on, head(h): the ones until t0.init
  t = d[new.init : length(d)]
  h = d[1 : new.init-1]
  
  # calculate the new pseudotime values
  new.t = c(min(pst),min(pst) + cumsum(t))
  new.h = tail(new.t,1) + cumsum(h)
  
  new.pst = c(new.h, new.t)
  return(new.pst)
}