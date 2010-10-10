poMoS <- function(series,nIter,samplFreq,nModAdd,nModSelec,nModStart=NULL,dt=1,show=1,toAnalyse=NULL,
			filterToExpl=NULL,filterReg=NULL,dMax=1,initSingle=F,initLog=NULL, critCalc="aic") {
	#########################################################################################
	#            										#
	# poMoS                                                      	     	                #
	# Recherche, en utilisant un critère proche de l'aic et une heuristique évolutionnaire,	#
	# un modèle de régression optimal pour les dérivées de la série étudiée			#
	#                                                                                       #
	#########################################################################################

	################################
	#Fonctions d'affichage diverses#
	################################
	marge<-1/5
	maxFenetre <- 0
	
	#Affiche l'évolution de la taille du panel de modèles à tester, en redimensionnant le graphique
	rePlotEvolNbMod <- function() {
		maxFenetre <<- max(vectNbMod)*(1+marge)
		plot(c(0,nIter),c(0,maxFenetre),type="n",xlab=gettext("Number of iterations"),ylab=gettext("Number of models"))
		title(gettext("Size evolution of the basket of model to test"))
		index <- as.vector(matrix(rep(1:length(vectNbMod[1,]),2),nrow=2,byrow=T))
		points(index,as.vector(vectNbMod),pch=19,col=rep(1:2,length(vectNbMod[1,])))
	}
	
	#Ajoute au graphique en cours les nouvelles données de taille du panel
	plotEvolNbMod <- function() {
		l <- length(vectNbMod[1,])
		if ( l == 1 | maxFenetre < vectNbMod[1,l] ) {
			rePlotEvolNbMod()
		} else {
			points(l,vectNbMod[1,l],pch=19)
			points(l,vectNbMod[2,l],pch=19,col=2)
		}
	}
	
	#Affiche la liste de régresseurs utilisée pour en retirer ou en ajouter
	showCheckBox <- function() {
		cacheR <- filterReg[,varBox$active+1]
		nRegModif <<- nReg * (!buttonCompact$active) + sum(cacheR) * buttonCompact$active
		try(for (i in 1:length(buttonReg)) {
			buttonReg[[i]]$destroy()
		},T)
		buttonReg <<- vector("list",nRegModif)
		for (i in 1:nRegModif) {
			buttonReg[[i]] <<- gtkCheckButtonNew()
			buttonReg[[i]]$label <<- paste(nomPoly[cacheR==1][i][buttonCompact$active],nomPoly[i][!buttonCompact$active],sep="")
			buttonReg[[i]]$active <<- buttonCompact$active | (cacheR[i]==1)
			regBox$add(buttonReg[[i]])
			gSignalConnect(buttonReg[[i]],"clicked",function(button) {
				for (k in 1:nReg) {
					if (button$label == nomPoly[k]) {
						break
					}
				}
				gtkWidgetSetSensitive(window,F)
				filterRegTemp[k,varBox$active+1] <<- 1 - filterReg[k,varBox$active+1]
				if(filterRegTemp[k,varBox$active+1]==0) {
					retireModDef <<- k + nReg*varBox$active
					log <<- paste(log,nomPoly[k],sep=", ")
				}
			})
		}
	}

	###########################
	#Autres fonctions diverses#
	###########################

	#Ajoute 'nModStart' modèles dans le panel des modèles à tester
	ajouterModeles <- function(nModStart) {
		if (nModStart > (2^(sum(filterReg)) - 1)^nVar) {
			stop("Too many models to add to the basket comparing to the number of regressors.")
		}
		if (!is.null(modelesATester$def)) {
			d <- dim(modelesATester$def)[2]
		} else {
			d <- 0
		}
		nModRestant <- nModStart
		while (0 < nModRestant) {	
			M <- matrix(c(modelesATester$def,round(runif(nVar*nReg*nModRestant,0,1))),nrow=nVar*nReg)
			M[filterReg==0,] <- 0
			M <- unique(M,MARGIN=2)
			cache <- t(rep(1,nVar) %*% (testZero %*% M^2 == 0))
			modelesATester$def <<- matrix(M[,cache==0],nrow=nVar*nReg)
			nModRestant <- nModStart - dim(modelesATester$def)[2] + d
		}
	}

	#Calcul le log de la vraissemblance d'un modèle
	vraisReg <- function(cacheMod) {		
		## Régressions linéaires ##
		res <- vector("numeric",length=NDecim*nVar)
		coeffStack <- NULL
		for (j in 1:nVar) {
				cacheVar <- 1:nReg + (j-1)*nReg
				cachePred <- 1:NDecim + (j-1)*NDecim
				regAUtiliser <- as.matrix(regDecim[,cacheMod[cacheVar]*filterReg[,j]==1])
				oldWarn <- options(warn=2)
				temp <- try(lsfit(regAUtiliser,derivesDecim[,j],intercept=FALSE))
				options(warn=oldWarn$warn)
				if (inherits(temp,"try-error")) {
					print(dim(regAUtiliser))
					print(j)
					print(labelPoly(nVarReg,dMax)[cacheMod[cacheVar]*filterReg[,j]==1])
					temp <- suppressWarnings(lsfit(regAUtiliser,derivesDecim[,j],intercept=FALSE))
				}
				coeffStack <- c(coeffStack,temp$coefficients)
				res[cachePred] = derivesDecim[,j] - (regAUtiliser %*% as.matrix(temp$coefficients))

		}
		## Test de normalité ##
		if (critCalc != "modif") {
			normTest <- ks.test(res[1:NDecim], mean=mean(res[1:NDecim]), sd = sqrt(var(res[1:NDecim])), "pnorm", alternative="two.sided")
			if(normTest$p.value < 0.01) {
				warning("The residuals of the linear regression are not gaussian at 99%, you should not use Akaike's or bayesian criterion.")
			}
		}
			
		## Matrice de variance-covariance ##
		matCov <- det(cov(matrix(res,nrow=NDecim)))
		list(vrais=log(matCov),param=coeffStack)

	}

	#Recheche la vraissemblance d'un modèle dans ceux déjà testés,
	#ou la calcul le cas échéant
	vraisSearch <- function(cacheMod) {
		cacheCrit <- t(cacheMod) %*% tested$def + t(1 - cacheMod) %*% (1 - tested$def) == nReg*nVar
		vrais <- tested$crit[cacheCrit] - log(sum(cacheMod))
		if (length(vrais)==0) {
			cacheCrit <- t(cacheMod) %*% modelesATester$def + t(1 - cacheMod) %*% (1 - modelesATester$def) == nReg*nVar
			vrais <- modelesATester$crit[cacheCrit] - log(sum(cacheMod))
			if (length(vrais)==0) {					
				vrais <- vraisReg(cacheMod)$vrais
			} else {
				if (is.na(vrais)) {
					listVrais <- vraisReg(cacheMod)
					modelesATester$coeff[cacheMod==1,cacheCrit] <<- listVrais$param
					modelesATester$crit[cacheCrit] <<- listVrais$vrais + log(sum(cacheMod)) 
					vrais <- listVrais$vrais
				}
			}				
		}
		vrais
	}
		
	##########################################
	#Initialisation des dérivées / decimation#
	##########################################
	if (show > 1) { message("Derivatives initialisation / decimation...")	} 
	nVarReg <- dim(series)[2]
	N <- dim(series)[1]
	nModTot <- c(0,0)
	derives <- (series[2:N,] -  series[1:(N-1),])/dt
	derivesDecim <- derives[(0:floor((N-2)/samplFreq))*samplFreq + 1,]
	seriesDecim <- series[(0:floor((N-2)/samplFreq))*samplFreq + 1,]	

	#Calcul des séries de régresseurs, décimées ou non	
	if (show > 1) { message("Regressor series creation...") }
	NDecim <- dim(derivesDecim)[1]
	dMaxP <- 0
	while ( choose(dMaxP+1+nVarReg,nVarReg) < NDecim-1) {
		dMaxP <- dMaxP + 1
	}
	if (dMax > dMaxP) {
		if (show > 0) { message("Too high polynomial degree comparing to sample size. To go on like this would mean to use more regressors than data.") }
		dMax <- dMaxP
	}
	reg <- seriesToPoly(series[1:(N-1),],dMax)
	reg <- reg[, rep(1,dim(reg)[1]) %*% reg^2 > 0]
	regDecim <- reg[(0:floor((N-2)/samplFreq))*samplFreq + 1,]
	nomPoly <- labelPoly(nVarReg,dMax)
	
	#Vérification de la cohérence des entrées
	if (!is.null(filterToExpl)) {
		derives <- as.matrix(derives[,filterToExpl])
		derivesDecim <- as.matrix(derivesDecim[,filterToExpl])
	} else {
		filterToExpl <- 1:nVarReg
	}
	nVar <- dim(derives)[2]
	nReg <- dim(regDecim)[2]
	filterToExplTemp <- filterToExpl
	if (is.vector(filterReg)) {
		filterReg <- as.matrix(filterReg)
	}
		if (is.null(filterReg)) {
		filterReg <- matrix(1,ncol=nVar,nrow=nReg)
	}
	if (dim(filterReg)[1] != nReg) {
		stop("The maximum polynomial degree ('dMax') do not match the dimensions of the regressors filter ('filterReg').")
	}
	if (dim(filterReg)[2] != nVar) {
		stop("The number of studied variables do not match the dimensions of the regressors filter ('filterReg').")
	}
	filterRegTemp <- filterReg
	log <- initLog
	if (!is.null(initLog)) {
		log <- paste(log,"|")
	}
	
	############################
	#Choix des modèles initiaux#
	############################
	if (show > 0) { message("Setting the intial models...") }
	testZero <- matrix( c(rep(0,(nVar - 1)*nReg), rep( c(rep(1,nReg), rep(0,(max(nVar-2,0))*nReg)), nVar-1 ),
				    rep(1,nReg), rep(0,(nVar - 1)*nReg)),nrow=nVar,byrow=T)		
	tested <- list()	
	tested$def <- NULL
	tested$crit <- NULL
	tested$coeff <- NULL
	modelesATester <- list()
	modelesATester$def <- toAnalyse;
	if (initSingle) {
		modelesATester$def <- cbind(modelesATester$def, diag(rep(1,nVar*nReg)))		
	} else {
		if (is.null(nModStart)) {
			stop("Enter the number of initial models or choose to take every models with only one regressor, using the option 'initSingle=T'.")
		}
		ajouterModeles(nModStart)		
	}
	modelesATester$crit <- vector("numeric",dim(modelesATester$def)[2])*NA
	modelesATester$coeff <- modelesATester$def*NA

	##########################
	#Initialisation graphique#
	##########################
	if (show > 1) { message("Graphical intialisation...") }
	continue <- T
	pause <- F
	changerDMaxEn <- Inf
	retireModDef <- 0
	finSel <- F	
	ajoutModUnReg <- F
	ajoutMod <- 0
	vectNbMod <- NULL
	nIterText <- as.character(nIter)
	if (show > 0) {	
		library(RGtk2)
		library(cairoDevice)
		
		### Import du fichier XML de l'interface graphique ###
		builder <- gtkBuilder()
		res <- builder$addFromFile(.uiFile("polyModInterface","PoMoS"))
		builder$connectSignals(NULL)
		window <- builder$getObject("windowMain")
		window$showAll()
		nRegModif <- nReg
		pageEnCours <- 0
		nSansAff <- 0
		freqAff <- Inf
		logBp <- F
		flagBoutonReg <- T
		epaisseurBut <- 30
		regBox <- builder$getObject("regBox")
		regView <- builder$getObject("regView")
		nModAddScale <- builder$getObject("adjustmentAj")
		nModAddScale$upper <- 2*nReg*nVar
		nModAddScale$lower <- nModSelec
		nModAddScale$value <- nModAdd
		nModSelScale <- builder$getObject("adjustmentSel")
		nModSelScale$upper <- nModAdd
		nModSelScale$value <- nModSelec
		freqAffScale <- builder$getObject("adjustmentBp")
		freqAffScale$upper <- nIter
		freqAffScale$value <- 0
		dMaxScale <- builder$getObject("adjustmentDMax")
		dMaxScale$upper <- dMaxP
		dMaxScale$value <- dMax
		dessAdj <- builder$getObject("adjustmentDess")
		dessAdj$upper <- sum(filterReg)
		dessAdj$lower <- 0
		dessAdj$value <- 1
		dessBpAdj <- builder$getObject("adjustmentDessBp")
		dessBpAdj$upper <- sum(filterReg)
		dessBpAdj$lower <- 1
		dessBpAdj$value <- sum(filterReg)	
		adjAjMod <- builder$getObject("adjustmentAjMod")
		acpAx1 <- builder$getObject("adjustmentAx1")
		acpAx1$value <- 1
		acpAx2 <- builder$getObject("adjustmentAx2")
		acpAx2$value <- 2		
		adjSsMod <- builder$getObject("adjustmentSsMod")
		adjGraph <- builder$getObject("adjustmentGraph")
		adjDist <- builder$getObject("adjustmentDist")
		adjDist$upper <- nReg
		adjDist$value <- 1
		nbModLab <- builder$getObject("labelNbModel")
		varLab <- builder$getObject("varLab")
		varLab$setText(paste(nVar,gettext("variable(s) :")))
		spinAx1 <- builder$getObject("spinAx1")
		spinAx2 <- builder$getObject("spinAx2")
		onglets <- builder$getObject("onglets")
		asCairoDevice(builder$getObject("evolNbModPlot"))
		idCairo <- vector("numeric",4)
		idCairo[1] <- max(dev.list())
		dev.set(idCairo[1])
		asCairoDevice(builder$getObject("boxPlotReg"))
		asCairoDevice(builder$getObject("convGraphPlot"))
		asCairoDevice(builder$getObject("ACPPlot"))
		acpList <- builder$getObject("ACPList")
		acpTreeMod <- acpList$model
		innerProdACP <- builder$getObject("innerProdACP")
		buttonReg <- NULL
		buttonPause <- builder$getObject("buttonPause")
		buttonStop <- builder$getObject("buttonStop")
		buttonIter <- builder$getObject("buttonIter")
		buttonBp <- builder$getObject("buttonBp")
		buttonRafrBp <- builder$getObject("buttonRafrBp")
		buttonModUnReg <- builder$getObject("buttonModUnReg")
		buttonAjMod <- builder$getObject("buttonAjMod")
		buttonCompact <- builder$getObject("buttonCompact")
		buttonRetAj <- builder$getObject("adjRetAj")
		buttonRetAj$active <- TRUE
		distBox <- builder$getObject("distBox")
		nbRegLab <- builder$getObject("nbRegLab")
		nbRegLab$setText(paste(gettext("Number of regressors :"),sum(filterReg)))
		entryIter <- builder$getObject("entryIter")
		entryIter$text <- nIterText
		spinBp <- builder$getObject("spinBp")
		spinDMax <- builder$getObject("spinDMax")
		varBoxBox <- builder$getObject("varBoxBox")
		varBox <- gtkComboBoxNewText()
		for (i in 1:nVar) {
			varBox$insertText(i,paste("X",filterToExpl[i],sep=""))
		}
		varBox$active <- 0
		varBox$show()
		varBoxBox$packStart(varBox,TRUE,FALSE,0)
		showCheckBox()
		regView$hide()
		regView$show()
		if (dim(derives)[2] != 1) {
			gtkWidgetSetSensitive(buttonModUnReg,F)
		}
		
		### Connexion des évènements aux fonctions à exécuter ###
		gSignalConnect(onglets,"switch-page",function(o,page,pageNum,userData) {
			pageEnCours <<- pageNum
			if (pageEnCours == 0) {
				dev.set(idCairo[1])
				rePlotEvolNbMod()		
			}
			if (pageEnCours == 1) {
				if (idCairo[2] == 0) {
					idCairo[2] <<- max(idCairo)+1
				}
				dev.set(idCairo[2])
			}
			if (pageEnCours == 3) {
				if (idCairo[3] == 0) {
					idCairo[3] <<- max(idCairo)+1
				}
				dev.set(idCairo[3])
			}
			if (pageEnCours == 4) {
				if (idCairo[4] == 0) {
					idCairo[4] <<- max(idCairo)+1
				}
				dev.set(idCairo[4])
				acpAx1$upper <<- length(c(modelesATester$crit[!is.na(modelesATester$crit)],tested$crit))
				acpAx2$upper <<- acpAx1$upper
			}
		})
		gSignalConnect(varBox,"changed",function(cBox) {
			regBox$visible <<- F
			showCheckBox()
			regBox$visible <<- T
		})
		gSignalConnect(buttonStop,"clicked",function(button) {
			continue <<- F
			if (show > 0) { message ("The user stops the algorithm.") }
		})
		gSignalConnect(buttonPause,"clicked",function(button) {
			pause <<- !pause
			if (pause) {
				buttonPause$label <<- gettext("Continue")
				if (show > 1) { message ("The user restarts the algorithm.") }
			} else {
				buttonPause$label <<- gettext("Pause")
				if (show > 1) { message ("The user pauses the algorithm.") }
			}
		})
		gSignalConnect(builder$getObject("buttonGraph"),"clicked",function(button) {
			status <- pause
			pause <<- T
			gtkWidgetSetSensitive(window,F)
			panelDef <- cbind(modelesATester$def[,!is.na(modelesATester$crit)],tested$def)
			sortedCrit <- sort(c(modelesATester$crit[!is.na(modelesATester$crit)],tested$crit),index.return=T)
			panelDefSel <- panelDef[,sortedCrit$ix[1:adjGraph$value]]
			if (buttonRetAj$getActive()) {
				methAdj <- "remAdd"
			} else {
				methAdj <- "permut"
			}
			trash <- convGraph(rbind(panelDefSel,sortedCrit$x[1:adjGraph$value]),titleText=gettext("Models neighbourhood graph"),
						 subTitleText=gettext("Red : low criterium (optimum), Yellow : high criterium"),adj=methAdj,limit=adjDist$value)
			pause <<- status
			gtkWidgetSetSensitive(window,T)
		})
		gSignalConnect(builder$getObject("buttonLog"),"clicked",function(button) {
			logBp <<- !logBp
			buttonBp$clicked()
		})
		gSignalConnect(buttonBp,"clicked",function(button) {
			status <- pause
			pause <<- T
			gtkWidgetSetSensitive(window,F)
			cacheMod <- !is.na(modelesATester$crit)
			cacheNbMod <- t(rep(1,nReg*nVar)) %*% modelesATester$def <= dessBpAdj$value
			cacheTot <- as.vector(cacheMod) & as.vector(cacheNbMod)
			cacheNbModTest <- t(rep(1,nReg*nVar)) %*% tested$def <= dessBpAdj$value
			filterRegBp <- rep(T,nReg*nVar) & ((!buttonCompact$active) | as.vector(filterReg))
			aAfficher <- try(list(def=cbind(modelesATester$def[filterRegBp,cacheTot],tested$def[filterRegBp,cacheNbModTest]),
						crit=c(modelesATester$crit[cacheTot],tested$crit[cacheNbModTest]),
						coeff=cbind(modelesATester$coeff[filterRegBp,cacheTot],tested$coeff[filterRegBp,cacheNbModTest])))
			if (inherits(aAfficher,"try-error")) {
				print(as.vector(cacheMod) & as.vector(cacheNbMod))
				print(dim(modelesATester$def))
				stop("Error creating boxplots")
			}
			groupe <- NULL
			if (buttonCompact$active && nVar > 1) {
				groupe <- sum(filterReg[,1])
				if (nVar > 2) {
					for (i in 1:(nVar-1)) {
						groupe <- c(groupe,groupe[i-1] + sum(filterReg[,i]))
					}
				}
			}
			bpRegress(aAfficher,plot=T,labels=rep(nomPoly,nVar)[filterRegBp],logPlot=logBp,
				    nVar=nVar,groupReg=groupe,modByNReg=dessAdj$value)
			pause <<- status
			gtkWidgetSetSensitive(window,T)
		})
		gSignalConnect(buttonModUnReg,"clicked",function(button) {
			ajoutModUnReg <<- T
		})
		gSignalConnect(buttonIter,"clicked",function(button) {
			if (suppressWarnings(!is.na(as.integer(nIterText)))) {
				nIter <<- as.integer(nIterText)	
				rePlotEvolNbMod()
			}
		})
		gSignalConnect(buttonRafrBp,"clicked",function(button) {
			if (freqAff == Inf) {
				freqAff <<- as.integer(spinBp$text)
				buttonRafrBp$label <<- gettext("Stop")
			} else {
				freqAff <<- Inf
				buttonRafrBp$label <<- gettext("Start")
			}
		})
		gSignalConnect(buttonCompact,"clicked",function(button) {
			showCheckBox()
		})	
		gSignalConnect(builder$getObject("buttonACP"),"clicked",function(button) {
			status <- pause
			pause <<- T
			gtkWidgetSetSensitive(window,F)
			panelDef <- cbind(modelesATester$def[,!is.na(modelesATester$crit)],tested$def)
			filterDouble <- !duplicated(panelDef,MARGIN=2)
			panelDef <- panelDef[,filterDouble]
			nMod <- dim(panelDef)[2]
			filterACP <- (panelDef %*% as.matrix(rep(1,nMod)) > 0) & (panelDef %*% as.matrix(rep(1,nMod)) < nMod)
			panelDef <- panelDef[filterACP,]
			panelCrit <- c(modelesATester$crit[!is.na(modelesATester$crit)],tested$crit)[filterDouble]
			aDessiner <- rbind(panelDef,-panelCrit)
			acp <- prcomp(t(aDessiner))
			if (nVar > 1) {
				labACP <- NULL 
				for (i in 1:nVar) {
					labACP <- c(labACP, paste(labelPoly(nVarReg,dMax)," -> X",i,sep=""))
				}
			} else {
				labACP <- labelPoly(nVarReg,dMax)
			}
			labACP <- labACP[filterACP]
			biplot(acp,choices=c(spinAx1$value,spinAx2$value),xlabs=rep("",dim(aDessiner)[2]),
				 ylabs=c(labACP,"Crit"))
			l <- length(labACP)
			coordVect <- acp$rotation %*% sqrt(diag(acp$sdev))[,c(spinAx1$value,spinAx2$value)]

			innerProd <- coordVect[1:l,] %*% as.matrix(coordVect[l+1,])
			longVect <- sqrt(coordVect[1:l,1]^2 + coordVect[1:l,2]^2)
			corre <- cor(as.matrix(-panelCrit),t(panelDef))
			nbMod <- panelDef %*% as.matrix(rep(1,dim(panelDef)[2]))
			modPrinc <- vraisSearch(as.vector(filterReg))
			aAfficher <- vector("list",l)
			for (i in 1:l) {
				aAfficher[[i]] <- list(labACP=labACP[i],innerProd=innerProd[i],longVect=longVect[i],corre=corre[i],nbMod=nbMod[i])
				if (adjSsMod$value == 0) {
					aAfficher[[i]]$test <- TRUE
				} else {
					filterTest <- filterACP
					filterTest[(1:(nReg*nVar))[filterACP][i]] <- F
					aAfficher[[i]]$test <- vraisSearch(filterTest) - modPrinc < qchisq(1 - adjSsMod$value,1)
					#print(modPrinc - vraisSearch(filterTest))
				}
			}
			modif <- T
			while (modif) {
				modif <- F
				for(i in 2:l) {
					if (aAfficher[[i]]$corre < aAfficher[[i-1]]$corre) {
						modif <- T
						temp <- aAfficher[[i]]
						aAfficher[[i]] <- aAfficher[[i-1]]
						aAfficher[[i-1]] <- temp
					}
				}
			}
			innerProdACP$clear()
			for (i in 1:l){
				innerProdPtr <- innerProdACP$append()$iter
				innerProdACP$set(innerProdPtr,0,aAfficher[[i]]$labACP,
							1, aAfficher[[i]]$innerProd,
							2, aAfficher[[i]]$longVect,
							3, aAfficher[[i]]$corre,
							4, aAfficher[[i]]$nbMod,
							5, aAfficher[[i]]$test)
			}
			acpList$getColumn(5)$setVisible(adjSsMod$value != 0)
			pause <<- status
			gtkWidgetSetSensitive(window,T)
		})	
		gSignalConnect(buttonAjMod,"clicked",function(button) {
			ajoutMod <<- adjAjMod$value	
		})
		gSignalConnect(builder$getObject("adjPerm"),"clicked",function(button) {
			distBox$sensitive <<- FALSE
		})
		gSignalConnect(buttonRetAj,"clicked",function(button) {
			distBox$sensitive <<- TRUE
		})
		gSignalConnect(builder$getObject("hscaleAj"),"value-changed",function(range) {
			nModAdd <<- range$getValue()
			nModSelScale$upper <<- nModAdd
		})
		gSignalConnect(builder$getObject("hscaleSel"),"value-changed",function(range) {
			nModSelec <<- range$getValue()
			nModAddScale$lower <<- nModSelec
		})
		gSignalConnect(entryIter,"changed",function(entry) {
			nIterText <<- entry$getText()
		})
		gSignalConnect(entryIter,"activate",function(entry) {
			buttonIter$clicked()
		})
		gSignalConnect(spinDMax,"value-changed",function(spin) {
			sVal <- suppressWarnings(as.integer(spin$value))
			if (!is.na(sVal)) {
				gtkWidgetSetSensitive(window,F)
				changerDMaxEn <<- sVal
			}
		})
		gSignalConnect(window,"delete-event",function(win,user.data) {			
			if (continue) {
				continue <<- F
				if (show > 0) { message ("The user stops the algorithm.") }
				win$hideOnDelete()
			} else {
				win$destroy()
			}
		})
	}

	#######################
	#Sélection des modèles#
	#######################
	if (show > 0) { 
		message("Comparing models...") 
	}
	k <- 1
	while (k < nIter & continue) {

		if (show > 1) { message(paste(gettext("Memory size filled (Mb) : "),memory.size(),sep="")) }
		gc()
		if (show > 1) { message(paste(gettext("Memory size filled (Mb) after gc :  "),memory.size(),sep=""))	}
		
		doEvent <- T
		while((pause || doEvent) &  continue & (show > 0)) {
			while (gtkEventsPending()) {
				gtkMainIterationDo(FALSE)
			}

			### Modification éventuelle de la liste des régresseurs ###
			if (sum(retireModDef) != 0) {
				cache <- modelesATester$def[retireModDef,]
				modelesATester$def <- modelesATester$def[,cache==0]
				modelesATester$coeff <- modelesATester$coeff[,cache==0]
				modelesATester$crit <- modelesATester$crit[cache==0]
				cache <- tested$def[retireModDef,]
				tested$def <- tested$def[,cache==0]
				tested$coeff <- tested$coeff[,cache==0]
				tested$crit <- tested$crit[cache==0]
				retireModDef <- 0
			}
			if (sum(filterReg - filterRegTemp) != 0) {
				filterReg <- filterRegTemp
				if (buttonCompact$active) {
					showCheckBox()
				}
			}
			if (changerDMaxEn != Inf) {
				if(dMax < changerDMaxEn) {
					pCache <- polyFilter(nVarReg,changerDMaxEn)
					nReg <- dim(pCache)[2]
					cacheD <- t(rep(1,nVarReg)) %*% pCache <= dMax
					oldModelesATester <- modelesATester
					nMod <- dim(oldModelesATester$def)[2]
					modelesATester <- list(def=matrix(0,nrow=nReg*nVar,ncol=nMod),
								     crit=oldModelesATester$crit,
								     coeff=matrix(NA,nrow=nReg*nVar,ncol=nMod))
					modelesATester$def[rep(cacheD,nVar),] <- oldModelesATester$def
					modelesATester$coeff[rep(cacheD,nVar),] <- oldModelesATester$coeff
					oldTested <- tested
					nMod <- dim(oldTested$def)[2]
					tested <- list(def=matrix(0,nrow=nReg*nVar,ncol=nMod),
							   crit=oldTested$crit,
							   coeff=matrix(NA,nrow=nReg*nVar,ncol=nMod))
					tested$def[rep(cacheD,nVar),] <- oldTested$def
					tested$coeff[rep(cacheD,nVar),] <- oldTested$coeff
					dMax <- changerDMaxEn
					oldReg <- reg
					reg <- matrix(0,nrow=N-1,ncol=nReg)
					reg[,cacheD] <- oldReg
					reg[,!cacheD] <- seriesToPoly(series[1:(N-1),],dMax,pCache[,!cacheD])
					reg <- reg[, rep(1,dim(reg)[1]) %*% reg^2 > 0]
					regDecim <- reg[(0:floor((N-2)/samplFreq))*samplFreq + 1,]
					oldCacheReg <- filterReg
					filterReg <- matrix(1,nrow=nReg,ncol=nVar)
					filterReg[cacheD,] <- oldCacheReg
					filterRegTemp <- filterReg
				} else {
					pCache <- polyFilter(nVarReg,dMax)
					cacheD <- t(rep(1,nVarReg)) %*% pCache <= changerDMaxEn
					nReg <- sum(cacheD)
					cacheMod <- try(t(rep(1,(dim(pCache)[2] - nReg)*nVar)) %*% modelesATester$def[!rep(cacheD,nVar),] == 0)
					if (inherits(cacheMod,"try-error")) {
						print((dim(pCache)[2] - nReg)*nVar)
						print(dim(modelesATester$def[!rep(cacheD,nVar),]))
					}
					modelesATester <- list(def=modelesATester$def[rep(cacheD,nVar),cacheMod],
								     crit=modelesATester$crit[cacheMod],
								     coeff=modelesATester$coeff[rep(cacheD,nVar),cacheMod])
					cacheMod <- try(t(rep(1,(dim(pCache)[2] - nReg)*nVar)) %*% tested$def[!rep(cacheD,nVar),] == 0)
					if (inherits(cacheMod,"try-error")) {
						print(dim(pCache)[2] - nReg)
						print(dim(tested$def[!cacheD,]))
					}
					tested <- list(def=tested$def[rep(cacheD,nVar),cacheMod],
							   crit=tested$crit[cacheMod],
							   coeff=tested$coeff[rep(cacheD,nVar),cacheMod])
					dMax <- changerDMaxEn
					reg <- reg[,cacheD]


					reg <- reg[, rep(1,dim(reg)[1]) %*% reg^2 > 0]
					regDecim <- reg[(0:floor((N-2)/samplFreq))*samplFreq + 1,]
					filterReg <- as.matrix(filterReg[cacheD,])
					filterRegTemp <- filterReg
				}
				nomPoly <- labelPoly(nVarReg,dMax)	
				showCheckBox()
				changerDMaxEn <- Inf
				testZero <- matrix( c(rep(0,(nVar - 1)*nReg), rep( c(rep(1,nReg), rep(0,(max(nVar-2,0))*nReg)), nVar-1 ),
					 		    rep(1,nReg), rep(0,(nVar - 1)*nReg)),nrow=nVar,byrow=T)
			}
			if (ajoutModUnReg) {
				for (i in 1:(nReg*nVar)) {
					if (as.vector(filterReg)[i] == 1) {
						modelesATester$def <- cbind(modelesATester$def, c(rep(0,i-1),1,rep(0,nReg*nVar-i)))
					}
				}
				modelesATester$crit <- c(modelesATester$crit, rep(NA,sum(filterReg)))
				modelesATester$coeff <- cbind(modelesATester$coeff, matrix(NA,ncol=sum(filterReg),nrow=nReg*nVar))			
				ajoutModUnReg <- F
			}	
			if (ajoutMod > 0) {
				ajouterModeles(ajoutMod)
				modelesATester$crit <- c(modelesATester$crit,vector("numeric",ajoutMod)*NA)
				modelesATester$coeff <- cbind(modelesATester$coeff,matrix(NA,ncol=ajoutMod,nrow=nVar*nReg))
				if (finSel) {
					finSel <- F
					message("Algorithm restarted.")
				}
				ajoutMod <- 0
			}
			nModAddScale$lower <- min(2*sum(filterReg),nModAddScale$lower)
			nModAddScale$upper <- 2*sum(filterReg)
			nModAddScale$value <- min(nModAddScale$value,nModAddScale$upper)
			nModSelScale$upper <- min(nModSelScale$upper,nModAddScale$value)
			nModSelScale$value <- min(nModSelScale$upper,nModSelScale$value)
			dessBpAdj$upper <- sum(filterReg)
			dessBpAdj$value <- min(dessBpAdj$value, dessBpAdj$upper)
			adjAjMod$upper <- 2^(sum(filterReg) - 1)^nVar
			adjAjMod$value <- min(adjAjMod$upper,adjAjMod$value)
			nbRegLab$setText(paste(gettext("Number of regressors :"),sum(filterReg)))
			adjGraph$upper <- length(tested$crit) + sum(!is.na(modelesATester$crit))
			adjGraph$value <- min(adjGraph$value,adjGraph$upper)
			nbModLab$setText(paste(gettext("Number of models (max : "),adjGraph$upper,") :",sep=""))
			adjDist$upper <- sum(filterReg)

			nATester <- dim(modelesATester$def)[2]
			nUniATester <- dim(unique(modelesATester$def,MARGIN=2))[2]						
			#Affichage de l'évolution de l'algorithme
			if (doEvent) {
				if (!finSel) {vectNbMod <- cbind(vectNbMod,c(nATester,nUniATester))}
				nSansAff <- nSansAff + 1
				if (pageEnCours == 0 & !finSel) {
					plotEvolNbMod()
				}
				if (pageEnCours == 1) {
					if (nSansAff >= freqAff) {
						nSansAff <- 0
						buttonBp$clicked()
					}
				}
			}
			doEvent <- F
			gtkWidgetSetSensitive(window,T)	
		}
		if (show == 0) {
			nATester <- dim(modelesATester$def)[2]
			nUniATester <- dim(unique(modelesATester$def,MARGIN=2))[2]						
		}

		### Critére d'arrêt ###
		if (show > 1) { 	message("Test to stop the algorithm...") }
		if (nUniATester < nModSelec) {
			if (show > 0 && !finSel) { message("No more models to test. End of the selection.") }
			if (show==0) { continue <- F }
			finSel <- T
		} else {
		
			### Choix des modèles qui vont s'affronter ###
			if (show > 1) { message("Choice of the models that will fight...") }
			nModRestant <- nModSelec
			compteSelec <- 0
			while (0 < nModRestant & compteSelec < 100*nModSelec) {	
				modelesSelec <- unique(ceiling(runif(nModSelec,0,nATester)))
				cache <- (as.matrix(modelesSelec) == 0)
				modelesSelec <- modelesSelec[cache==0]
				nModRestant <- nModSelec - length(modelesSelec)
				compteSelec <- compteSelec + 1
			}
			if (nModRestant > 0) {
				warning("Cannot find models to be tested.")
			}

			### Calcul du critère ###
			if (show > 1) { message("Criterium computation...") }
			cache <- is.na(modelesATester$crit[modelesSelec])
			indices <- modelesSelec[(cache * 1:length(cache))[(cache * 1:length(cache))>0]]
			cacheMod <- as.matrix(modelesATester$def[,indices]) 
			if (length(indices) != 0) {
				for (i in 1:dim(cacheMod)[2]) {
					listVrais <- vraisReg(cacheMod[,i])
					modelesATester$coeff[as.vector(filterReg)*cacheMod[,i]==1,indices[i]] <- listVrais$param
					if (critCalc == "aic") {
						modelesATester$crit[indices[i]] <- (dim(derivesDecim)[1])*listVrais$vrais + 2*sum(cacheMod[,i])
					} else if (critCalc == "bic") {
						modelesATester$crit[indices[i]] <- (dim(derivesDecim)[1])*listVrais$vrais + sum(cacheMod[,i])*log(dim(derivesDecim)[1])
					} else {
						modelesATester$crit[indices[i]] <- listVrais$vrais + log(sum(cacheMod[,i])) 
					}
				}
			}

			### Sélection du (des) meilleur(s) modèle(s) et modification des ensembles de modèles en conséquence ###
			if(show > 1) { message("Finding the best criterium...") }
			bestIndice <- modelesSelec[min(modelesATester$crit[modelesSelec]) == modelesATester$crit[modelesSelec]][1]
			mauvaisIndices <- modelesSelec[modelesSelec != bestIndice]
			tested$def <- cbind(tested$def, modelesATester$def[,mauvaisIndices])
			tested$crit <- c(tested$crit, modelesATester$crit[mauvaisIndices])
			tested$coeff <- cbind(tested$coeff, modelesATester$coeff[,mauvaisIndices])
			
				## Ajout des nouveaux modèles ##
			if (show > 1) { message("Adding new models to the basket...") }
			modifsModels <- ceiling(runif(nModAdd,0,sum(filterReg))) + 0:(nModAdd-1)*sum(filterReg)
			modelesAAjouter <- as.matrix(modelesATester$def[filterReg==1,bestIndice]) %*% t(rep(1,nModAdd))
			modelesAAjouter[modifsModels] <- 1 - modelesAAjouter[modifsModels]
			modelesAAjFin <- matrix(0,ncol=nModAdd,nrow=nVar*nReg)
			modelesAAjFin[filterReg==1,] <- modelesAAjouter
			M <- cbind(tested$def,modelesAAjFin)
			M <- unique(M,MARGIN=2)
			nTested <- dim(unique(tested$def,MARGIN=2))[2]
			cache <- t(rep(1,nVar) %*% (testZero %*% M^2 == 0))
			cache[1:nTested] = 1
			modelesAAjFin <- M[,cache==0]
			nModAddVrai <- dim(as.matrix(modelesAAjFin))[2]
		
				## Actualisation des modèles à tester ##
			if (show > 1) { message("Refreshing the basket of models...") }
			temp <- modelesATester$def
			temp[,mauvaisIndices] <- 0
			indConserves <- t(rep(1,nVar*nReg)) %*% temp > 0
			if (nModAddVrai > 0) {
				modelesATester$def <- cbind(modelesATester$def[,indConserves], as.matrix(modelesAAjFin))
				modelesATester$crit <- c(modelesATester$crit[indConserves], as.matrix(rep(NA,nModAddVrai)))
				modelesATester$coeff <- cbind(modelesATester$coeff[,indConserves], matrix(NA,ncol=nModAddVrai,nrow=nReg*nVar))			
			} else {
				modelesATester$def <- as.matrix(modelesATester$def[,indConserves])
				modelesATester$crit <- as.matrix(modelesATester$crit[indConserves])
				modelesATester$coeff <- as.matrix(modelesATester$coeff[,indConserves])
			}
			k <- k + 1
		}
	}
	continue <- 0
	testDestruct <- try(window$visible,TRUE)
	if (!inherits(testDestruct,"try-error")) {
		if (!testDestruct) {
			window$destroy()
		}
	}
	resultat <-  list(selectedMod=modelesATester,rejectedMod=tested,dMax=dMax,filterReg=filterReg,log=log)
}
