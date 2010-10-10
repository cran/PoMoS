convGraph <- function(data,titleText=gettext("New graph"),
			     subTitleText="",
			     adj="remAdd",limit=1) {
	#Représente les modèles testés par recherche d'AIC. renvoie les objets 
	# graph et layout d'igraphs associés

	library(igraph)

	#Traitement des données
	nReg <- dim(data)[1] - 1
	data <- unique(data, MARGIN=2)
	def <- data[1:nReg,]
	crit <- data[nReg+1,]
	colInd <- (crit - min(crit))/(max(crit)-min(crit))
	col <- rgb(1,colInd,0)
	opti <- crit == min(crit)
	
	#Création de la matrice d'adjacence
	if (adj=="remAdd") {
		M <- (t(def) %*% def + t(1-def) %*% (1-def))
		M[M < nReg - limit | M == nReg] <- 0
		M[M >= nReg - limit & M != nReg] <- M[M >= nReg - limit & M != nReg] / (nReg - 1)
		#M <- ((t(def) %*% def + t(1-def) %*% (1-def)) == nReg-1)
	} else if(adj=="permut") {
		sumDef <- as.matrix(rep(1,dim(def)[2])) %*% t(t(rep(1,nReg) %*% def))
		M <- ((t(def) %*% def + t(1-def) %*% (1-def)) == nReg-2) & (sumDef == t(sumDef))
	} else {
		stop("The computation method of the adjacency matrix ('adj') is not correct.")
	}
	

	#Création du graphe
	gr <- graph.adjacency(M,weighted=T)
		
	#Positionnement des points dans l'espace
	lay <- layout.fruchterman.reingold(gr,weights=E(gr)$weight)

	#Affichage du graphe
	plot(gr,layout=lay,vertex.size=2 + 2*opti,vertex.color = col,vertex.frame.color=col,vertex.label=NA,edge.arrow.mode=0)
	title(titleText,subTitleText)	

	resultat <- list(graph=gr,layout=lay)

}