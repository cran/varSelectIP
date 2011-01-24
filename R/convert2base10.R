convert2base10 <- function(x) {
  return(sum(2 ^ (which(as.logical(rev(x))) - 1)))
}
