# function to write to file
writeToFile <- function(models.visited, fname="models.csv", keep, p, 
                        m0.z) {
  tmp1 <- models.visited[complete.cases(models.visited),] # get rid of NA's
  tmp1 <- unique(tmp1)
  tmp1 <- tmp1[order(tmp1[,"BF"], decreasing=TRUE),][1:keep, ]
  # Compute m_gamma(z) / m_0(z)
  if (!is.null(m0.z)) {
    tmp1[ ,"BF"] <- tmp1[, "BF"]/m0.z
  }
  tmp2 <- t(apply(as.matrix(tmp1[,1]), 1, convert2base2, p))
  tmp2 <- cbind(tmp2, tmp1[,"BF"])

  colnames(tmp2) <- c(sprintf("cov%02d", 1:p), "Bayes.Factor")

  write.csv(tmp2, fname, quote=FALSE)
  tmp2
}
