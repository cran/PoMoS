\name{bpRegress}
\alias{bpRegress}
\title{Boxplot representation of models criterion.}

\description{Creates boxplots of the criterions, sorted by regressors, obtained for the whole set of tested models. 
This function is used by \code{poMoS} in graphical mode.}

\usage{bpRegress(models,labels=NULL,plot=TRUE,modByNReg=1,logPlot=FALSE,nVar=1,groupReg=NULL)}

\arguments{
	\item{models}{A list of elements providing information about tested models such as their criterion or their structure. A detailed description
	of this list can be found in the documentation of \code{\link{poMoS}}}
	\item{labels}{Text to be shown below each column of the graphic. As an example, it can be the output of labelPoly \code{\link{labelPoly}} }
	\item{plot}{A boolean which indicate if boxplots must be represented graphically or not.}
	\item{modByNReg}{When a graphical output is required and if this type of model exists, \code{bpRegress} also draws 
	black points to highlight models that contain every regressors but the one of the column. Models with the regressor 
	of the column plus a fixed number of others can also be shown (with color points); the maximum number of regressors 
	to be represented by points can be interactively fixed. Thus, to be able to see the criterion of models with only 
	the regressor of the column, those with it and one other, and those with it and two others, the user should enter \code{modByNReg=3}}
	\item{logPlot}{Determine if criteria must be shown as they are or in logarithm scale. It is usefull to
	choose the last one if data are packed at the bottom of the graphic.}
	\item{nVar}{In order to separate regressors of different variables by a line, nVar should be set to 1. The columns of boxplots are then divided into \code{nVar} groups of equal size. 
	If there is not the same number of regressors for each variable, use \code{groupeReg}.}
	\item{groupReg}{The vector of cumulated number of regressors for each variable. It can also be defined as the vector of the choosen
	boxplots indexes to draw separation lines on their right. It is related to the code{nVar} parameter, but has priority over it.}
}

\value{A list with one element by regressor, which is a vector of criteria of the models within this regressor is included.}

\note{Boxplot colours depend on the relative size of the sample plotted. The larger data length, the darker, the smaller data length, the lighter.}

\author{Mangiarotti S., Coudret R., Drapeau L.}

\seealso{\link{poMoS}, \link{labelPoly} }

\examples{
	#Sprott K attractor loading
	data(sprottK,package="PoMoS")
	sprottK <- as.matrix(sprottK)
	
	# We only select regressors found in 'labelPoly' example for the first variable
	filter <- c(FALSE, rep(TRUE,4), rep(FALSE,2), rep(TRUE,3))
	filterTot <- cbind(filter, rep(TRUE,10))
	labelsForBp <- c(labelPoly(3,2,toFind=filter),labelPoly(3,2))
	
	#We search for the best model and then draw boxplots
	res <- poMoS(sprottK,200,20,20,10,nModStart=100,filterToExpl=c(1,2),show=0,dMax=2,filterReg=filterTot,critCalc="modif")
	toPlot <- list(def=res$rejectedMod$def[as.vector(filterTot),],
				   crit=res$rejectedMod$crit,
				   coeff=res$rejectedMod$coeff[as.vector(filterTot),])
	trash <- bpRegress(toPlot,labels=labelsForBp,groupReg=sum(filter))	
}
