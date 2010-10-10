\name{findRegData}
\alias{findRegData}
\title{Information searching from poMoS output}

\description{Find a model structure from a list of results produced by \code{poMoS}}

\usage{findRegData(models, toFind, neigh=FALSE) }

\arguments{
	\item{models}{A list of information as provided by \code{\link{poMoS}} function. See the documentation of this function for more information.}
	\item{toFind}{A vector of model definition, as long as there are 
	regressors to take into consideration, with 1 if the regressor is used
	and 0 otherwise. To see which element of the vector corresponds to
	which regressor, \code{\link{labelPoly}} can be used.
	}
	\item{neigh}{The maximal distance between the model to be found \code{toFind} and the
	obtained from the model search. To enter \code{FALSE} is equivalent to use 0 for this
	distance. That means that only the informations concerning \code{toFind}
	is wanted. The distance is defined as the minimum number of modification to be applied to a model 
	(which structure is defined thanks to a vector) to obtain the other.}
}

\value{A list of information like \code{\link{poMoS} output but containing the selected models, only. }
}

\author{
Mangiarotti S., Coudret R., Drapeau L.
}

\seealso{
	\code{\link{poMoS}}.
}

\examples{
	data(sprottK,package="PoMoS")
	res <- poMoS(as.matrix(sprottK),100,20,4,3,nModStart=15,filterToExpl=3,dt=1/20,dMax=1,show=0,critCalc="modif")
	findRegData(res$rejectedMod,c(0,1,1,0))
}
