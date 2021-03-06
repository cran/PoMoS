findRegData <- function(models,toFind,neigh=FALSE) {

	#Renvoie les informations (d�finition, valeur des param�tres, valeur du crit�re)
	#d'un mod�le dont la d�finition est fournie. Possibilit� d'obtenir �galement
	#les informations des mod�les proches

	nReg <- dim(models$def)[1]
	resultat <- list()
	resultat$coeff <- unique(as.matrix(models$coeff[,t(toFind) %*% models$def + t(1-toFind) %*% (1-models$def) >= nReg - neigh]),MARGIN=2)
	resultat$def <- unique(as.matrix(models$def[,t(toFind) %*% models$def + t(1-toFind) %*% (1-models$def) >= nReg - neigh]),MARGIN=2)
	resultat$crit <- unique(models$crit[t(toFind) %*% models$def + t(1-toFind) %*% (1-models$def) >= nReg - neigh])
	resultat

}