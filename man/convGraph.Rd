\name{convGraph}
\alias{convGraph}
\title{Graphical plot of neighboring models}

\description{Draws models as points and links them if their structure is almost the same. 
The nearness of models connected can be set. This function can be launched graphically with \link{poMoS}.}

\usage{convGraph(data,titleText=gettext("New graph"), subTitleText="", adj="remAdd",limit=1)}

\arguments{
	\item{data}{A matrix in which informations to be represented are stored. Each column correpsonds to a model.
The first rows follow the same pattern as the \code{def} element of the list presented in the note of the documentation of \link{poMoS}.
The last one contain criteria.}
	\item{titleText}{The title of the graph.}
	\item{subTitleText}{The subtitle of the graph.}
	\item{adj}{A string to choose the distance that determines neighbourhood relations. With \code{"remAdd"}, the distance refers to the number 
	of regressors we have to add or remove to transform a model into another. With \code{"permut"}, the distance refers to the minimum number of permutation. 
	A permutation is defined as the action of remove a regressor  and then add another one. As a result, when adj="permut", two models of different 
	sizes have an infinite distance.}
	\item{limit}{Available if \code{adj="remAdd"}, only. Sets the maximum distance at which models are linked together. Anyway, models separated with small
	distances will appear closer from each other, because the representation of the relations depends on distances.}
}

\value{ A list of two elements :
		\item{\code{graph}}{An graph object as it is defined in the \code{igraph} package. It contains vertices that represents models
and edges which are links between them.}
		\item{\code{layout}}{For the user to understand which models are close to the others, the function calls an iterative algorithm to
display them the most clearly it is possible. This object sends the models position in a 2D space after this optimisation.} 
}

\note{The graph produced shows coloured points. The red ones are those with lowest criterion, which are consquently the optimum ones. As for the
colour yellow, it indicates highest criterion. The best model is shown by a big point.}

\author{
Mangiarotti S., Coudret R., Drapeau L.
}

\seealso{\link{poMoS}}

\examples{
	#Rossler attractor loading
	data(rossler,package="PoMoS")
	rossler <- as.matrix(rossler)
	
	#We try to explain the third series
	res <- poMoS(rossler,100,20,50,5,nModStart=100,filterToExpl=3,dMax=3,show=0,critCalc="modif")
	toPlot <- rbind(res$rejectedMod$def,res$rejectedMod$crit)
	trash <- convGraph(toPlot,"Tested models graph","Maximum Distance : 3",limit=3)
}
