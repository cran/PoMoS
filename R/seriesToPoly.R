seriesToPoly <- function(series, dMax=NULL, pFilter=NULL) {

	if (is.vector(series)) {
		series <- t(series)
	}
	if (is.data.frame(series)) {
		series <- as.matrix(series)
	}
	nVar <- dim(series)[2]
	N <- dim(series)[1]

	#D�termination des puissances auxquelles �lever les s�ries
	if (is.null(pFilter)) {
		if (is.null(dMax)) {
			stop("'dMax' or 'pFilter' is required.")
		}
		pFilter <- polyFilter(nVar,dMax)
	} else {
		if (is.vector(pFilter)) {
			pFilter <- as.matrix(pFilter)
		}
	}

	#Calcul des s�ries de polyn�mes des s�ries initiales
	s <- (-1) ^ ((sign(series) < 0) %*% pFilter)
	seriesPoly <- s * exp(log(abs(series)) %*% pFilter)
	seriesPoly
}