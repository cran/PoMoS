labelPoly <- function(nVar,dMax,toFind=NULL) {

	pMax <- choose(dMax+nVar,nVar)
	labels <- vector("character",pMax)
	labels[1] <- "Ct"
	j <- vector("numeric",nVar)
	for (i in 2:pMax) {
		j[nVar] <- j[nVar] + 1
		while (sum(j) > dMax) {
			l <- ((1:nVar)[j>0])
			j[l[length(l)]] <- 0
			j[l[length(l)]-1] <- j[l[length(l)]-1] + 1
		}	
		for (k in 1:nVar) {
			if (j[k] > 1) {
				labels[i] <- paste(labels[i],"X",k,"^",j[k]," ",sep="")
			}
			if (j[k] == 1) {
				labels[i] <- paste(labels[i],"X",k," ",sep="")
			}
		}	
	}
	if (!is.null(toFind)) {
		labels <- labels[toFind]
	}
	labels

}