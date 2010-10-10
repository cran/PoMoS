polyFilter <- function(nVar,dMax) {

	pMax <- choose(dMax+nVar,nVar)
	res <- matrix(0,ncol=pMax,nrow=nVar)
	j <- vector("numeric",nVar)
	for (i in 2:pMax) {
		j[nVar] <- j[nVar] + 1
		while (sum(j) > dMax) {
			l <- ((1:nVar)[j>0])
			j[l[length(l)]] <- 0
			j[l[length(l)]-1] <- j[l[length(l)]-1] + 1
		}	
		res[,i] <- j 	
	}
	res

}