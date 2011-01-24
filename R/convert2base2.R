# helper fn
convert2base2 <- function(x, q) {
  ndigits <- 1 + floor(log(x, 2)) 
  r <- rep(0, times=ndigits)
  if (ndigits >= 1) 
      for (i in ndigits:1) {
          r[i] <- x%%2
          if (i > 1) 
              x <- x%/%2
      }
  if (ndigits < q) {
    r <- c(rep(0, times=q-ndigits), r)
  }
  r
}
