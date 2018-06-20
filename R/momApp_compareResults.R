#======== todo =================================================================
#t1 documentation for compareMomentClosure

compareMomentClosure <- function(model, l = 4, m = 2){
  closureMethods <- c("setZero", "reduceDegree")
  print(paste("l = ", m*l))
  results <- lapply(closureMethods, function(i) {
    print(i)
    momApp(model, l = m*l, closure = i)
  })
  
  n <-  length(model@init) - 1
  s <- NULL
  for(i in 1:n){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+length(closureMethods))
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  colnames(s) <- c(names(model@init)[-(n+1)], closureMethods)
  s <- as.data.frame(s)
  
  for(i in 1:(n*l)){ # fill the moments
    for (j in 1:length(closureMethods))
      s[i, n+j] <- results[[j]]$contRes[nrow(results[[j]]$contRes), 
                                        prodlim::row.match(as.vector(s[i,1:n]), 
                                                           as.matrix(results[[j]]$contInd))]
  }
  s <- cbind(s, "|difference|" = abs(s[,n+1]-s[,n+2]))
  s
}

# für Modell 7 wird m = 3 benötigt (dann max Unterschied 7e-03)
# bei Modell 8 funktioniert es nicht (Differenzen sind bei m=3 größer als bei m=2)
# bei Modell 8 muss "times" viel größer sein (nichtlineare Dgln)! Tipp: fromtoby(c(from = 0, to = 400, by = 1))
