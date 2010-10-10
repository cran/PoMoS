bpRegress <- function(models,labels=NULL,plot=TRUE,modByNReg=1,logPlot=FALSE,nVar=1,groupReg=NULL) {

	#Trace les boxplots des valeurs des critères des modèles testés
	nReg <- dim(models$def)[1]
	if (nReg %% nVar != 0 & is.null(groupReg)) {
		warning("Incorrect number of variables.")
		nVar <- 1
	}
	indicesClasses <- list()
	indicesMoinsUnReg <- NA*vector(length=nReg)
	tailleEchant <- vector(length=nReg)
	indices <- models$crit
	minCrit <- min(indices)
	maxCrit <- max(indices)	
	if (logPlot) {
		modif <- function(x) {
			res <- log(x - minCrit + 1)
			res
		}
	} else {
		modif <- function(x) { x }
	}
	for (i in 1:nReg) {
		cache <- c(rep(0,i-1),1,rep(0,nReg-i))
		temp <- indices[models$def[i,]==1]
		indicesClasses[[i]] <- unique(modif(temp))
		b <- findRegData(models,1-cache)
		if (length(b$crit) > 0) {
			indicesMoinsUnReg[i] <- b$crit
		}
		tailleEchant[i] <- length(indicesClasses[[i]])
	}
	boxplot(indicesClasses,names=labels,col=rgb(1-tailleEchant/max(tailleEchant),1-tailleEchant/max(tailleEchant),1))
	title(gettext("Criterion values by regressor"),paste(gettext("Total number of observations :"),length(unique(indices))),
		 xlab=gettext("Regressors"))
	if (logPlot) {
		title(ylab=gettext("Criteria logarithm"))
	} else {
		title(ylab=gettext("Criteria"))
	}
	for (j in modByNReg:1) {
		for (i in 1:nReg) {
			cache <- c(rep(0,i-1),1,rep(0,nReg-i))
			a <- findRegData(models,cache,neigh=j-1)
			a$crit <- a$crit[t(rep(1,nReg)) %*% a$def == j]
			if (length(a$crit) > 0) {
				points(rep(i,length(a$crit)),modif(a$crit),col=j+1,pch=19)
			}
		}
	}
	points(1:nReg,modif(indicesMoinsUnReg),col=1,pch=19)
	if (is.null(groupReg)) {
		if (nVar > 1) {
			for (i in 1:(nVar-1)) {
				lines(rep((nReg/nVar)*i+0.5,2),c(modif(minCrit),modif(maxCrit)))
			}
		}
	} else {
		for (i in 1:length(groupReg)) {
			lines(rep(groupReg[i]+0.5,2),c(modif(minCrit),modif(maxCrit)))
		}
	}
	indicesClasses
}