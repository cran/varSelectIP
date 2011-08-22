
# Function to update the table of models                       #
################################################################
updateTable <- function(cur.table, model, BF) {
  keep <- dim(cur.table)[1]
  BFcol <- dim(cur.table)[2]
  min.cur.BF <- min(cur.table[,BFcol])
  which.min.cur.BF <- which.min(cur.table[,BFcol])

# If updating is needed
  if(BF>min.cur.BF) {
    matching.BF <- which(BF == cur.table[,BFcol])
    if (length(matching.BF)==0) {
      cur.table[which.min.cur.BF,] <- c(model, BF)
    }
    else{
      dup <- FALSE
      for (i in matching.BF) {
        if (sum((cur.table[i,-BFcol] - model)^2)==0) {
          dup <- TRUE
        }
      }
      if (dup==FALSE) {
        cur.table[which.min.cur.BF,] <- c(model, BF)
      }
    }
  }

  cur.table
}
