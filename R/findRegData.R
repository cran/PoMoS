findRegData <- function(models,toFind,neigh=FALSE) {

	#Renvoie les informations (définition, valeur des paramètres, valeur du critère)
	#d'un modèle dont la définition est fournie. Possibilité d'obtenir également
	#les informations des modèles proches

	nReg <- dim(models$def)[1]
	resultat <- list()
	resultat$coeff <- unique(as.matrix(models$coeff[,t(toFind) %*% models$def + t(1-toFind) %*% (1-models$def) >= nReg - neigh]),MARGIN=2)
	resultat$def <- unique(as.matrix(models$def[,t(toFind) %*% models$def + t(1-toFind) %*% (1-models$def) >= nReg - neigh]),MARGIN=2)
	resultat$crit <- unique(models$crit[t(toFind) %*% models$def + t(1-toFind) %*% (1-models$def) >= nReg - neigh])
	resultat

}