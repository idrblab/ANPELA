reverse_pseudotime <- function(pst){
  # pst = pst[order(pst)]
  d = diff(pst)
  d = rev(d)
  reverse_pseudotime = rev(cumsum(c(min(pst),d)))
  return(reverse_pseudotime)
}