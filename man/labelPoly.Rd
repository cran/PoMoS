\name{labelPoly}
\alias{labelPoly}
\title{Polynomial labels systematic ordering in a set of N variables}

\description{Generate a systematic correspondence of polynomial terms for a polynomial 
of degree \code{dMax} with nVar variables, and displays the terms in a formal way. This correspondence is used in the package \code{poMoS}}

\usage{labelPoly(nVar, dMax, toFind=NULL)}

\arguments{
\item{nVar}{The number of time series from which the polynomial series are built.}
\item{dMax}{The maximum degree of the polynomial function. The same \code{dMax} must be provided to \code{poMoS} when the functions are used together.}
\item{toFind}{It is allowed to select only one or several polynomial terms among those returned automatically. 
This selection can be achieved by using toFind that filters the results before the end of the function. 
It can be a vector of subscripts or booleans : \code{TRUE}  for the selected labels and FALSE elsewhere. 
It is usefull if someone only wants the regressors of a model and is not interested in all the tested ones.}
}

\value{ A vector of characters. Each element is the expression of a regressor such as \eqn{X_1^2 X_3 X_4}
}

\author{Mangiarotti S., Coudret R., Drapeau L.}

\seealso{\link{poMoS}}

\examples{
	data(sprottK,package="PoMoS")
	res <- poMoS(as.matrix(sprottK),200,20,5,4,nModStart=200,filterToExpl=1,dMax=2,show=0,critCalc="modif")
	testedMod <- list()
	testedMod$def <- res$selectedMod$def[,!is.na(res$selectedMod$crit)]
	testedMod$crit <- res$selectedMod$crit[!is.na(res$selectedMod$crit)]
	bestModDef <- unique(testedMod$def[,testedMod$crit == min(testedMod$crit)],MARGIN=2)
	labelPoly(3,2,toFind=bestModDef==1)
	#Note : The data in 'x' represent a Sprott H attractor. The returned regressors may not be
	#the initial ones because this use of poMoS gives a lot of false positives.
	#For more accuracy, this function should be launched in graphical mode, with boxplots.
	#In our tests we find with this method :
	# "X3 " "X3^2 " "X2 " "X2 X3 " "X1 X3 " "X1 X2 " "X1^2 " 
	#And the initial regressors are :
	# "X2 " "X1 X3 "
}
