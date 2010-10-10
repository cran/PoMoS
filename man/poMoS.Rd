\name{poMoS}
\alias{poMoS}
\title{Polynomial Model Search}

\description{poMoS function aims at determining, from a set of N time series, the optimal polynomial structure for a model built 
on first-order ordinary differential equations. An evolutionary algorithm combined with a least square fitting is 
used for this purpose. The optimality is estimated with AIC or AIC-like criterions.
Methodology: A simple evolutionary algorithm is used to generate sets of models structures and to test their efficiency 
to model the time series provided in input. Models are firstly identified through a least square fitting (the \code{lsfit} function is used for this purpose). Their performance is thus estimated with AIC-like criterions.}

\usage{poMoS(series, nIter, samplFreq, nModAdd, nModSelec, nModStart=NULL,
dt=1, show=1, toAnalyse=NULL, filterToExpl=NULL, filterReg=NULL, dMax=1,
initSingle=F, initLog=NULL, critCalc="aic")}

\arguments{
\item{series}{The matrix of time series data. Each column contains a time series.}
\item{nIter}{An integer that determines the maximum number of iterations for the algorithm. 
The algorithm stops earlier if in the pool of models \code{nModSelec} is greater than the number of models available for test. }
\item{samplFreq}{Subsampling frequency of the input series. If \code{samplFreq} is too small, there can be noise in the analysis that 
involves false positives in the regressors selection. If it is too big, the effect of the useful regressors is reduced and may become undetectable.}
\item{nModAdd}{How many models the algorithm adds to the pool at the end of each iteration. In practice, the number of added number 
can be less than \code{nModAdd} if not enough satisfying models are available.}
\item{nModSelec}{Determines how many models are compared during each iteration. Afterwards it eliminates the \code{nModSelec-1} worst models from the pool of models to test.}
\item{nModStart}{Initial number of models in the pool of models to be tested.}
\item{dt}{Time sampling of the input series. Used for the estimation of the regression coefficients.}
\item{show}{The depth of the interaction between the user and the function. It can be one of the following values :
	\describe{
		\item{0: Silent mode. No user interface, no messages.}{}
		\item{1: Normal mode. User interface, a few messages.}{}
		\item{2: Debug mode. User interface, a lot of messages.}{}
	}	
}
\item{toAnalyse}{Defines the more general model's structures to be considered in the analysis. 
These models are automatically included in the testing pool when starting the algorithm. This 
structure is given by a matrix: each column corresponds to an equation, each line corresponds 
to an algebraic regressor. Value 1 indicates that the regressor is selected; value 0 indicates 
that it is not. The correspondence between lines and regressors can be obtained from labelPoly 
providing the number nVar of time series and the maximum degree\code{dMax} of the polynomial.
}
\item{filterToExpl}{A vector setting that defines which series modeling to be analyzed. 
For example, with filterToExpl=c(2,3), the modeling of the second and third time series is considered. 
The derivatives of the time series 2 and 3 is thus computed and if no restriction is given (filterReg=NULL),
the analysis of these two series is carried on with the whole set of polynomial regressors. If the user 
does not fill this field, all the time series are analyzed.}
\item{filterReg}{A matrix that defines which regressors to be considered in the analysis. Each column 
corresponds to one series as defined in filterToExpl; the number of columns and the length of filterToExpl
 must thus be identical. Each line corresponds to one algebraic regressor. Value is 1 if the regressor
 is available, 0 if it is not. In order to ensure an iterative use of the algorithm, note that the same
 format is used for the output variable def. To set filterReg=NULL is equivalent to setting a matrix full of 1.}
\item{dMax}{The maximum degree of the polynomials. For example, let call x_1, x_2 and x_3 our 3 times series. 
The regressor x_1 x_2^3 x_3^2 is not part of the models if dMax = 6.}
\item{initSingle}{A boolean that indicates if the pool of models has to be initialized by all the models with only one
regressor. If \code{initSingle=T}, \code{nModStart} does not have to be set.}

\item{initLog}{It can be interesting for the user to keep a log of regressors removed after several launches of \code{poMoS}. To do so
, he or she can set \code{initLog} to the previous version of the element log of the list returned by the last execution of \code{poMoS}.}
\item{critCalc}{A word that defines which criteria to consider when estimating the models skill. Criteria "aic" and "bic" are available, 
as well as "modif". The latter is better adapted for deterministic series with non gaussian residuals arising from the derivative 
approximation. It can be called by setting \code{critCalc="modif"}.}}

\value{ A list of five elements :
		\item{\code{rejected}}{ Information about rejected models.}
		\item{\code{selected}}{ Information about both selected and untested models.} 
		\item{\code{dMax}}{The value of dMax at the end of the run. It can be different from the input value 
		if it was interactively changed by the user through the interface.}
		\item{\code{filterReg}}{The value of \code{filterReg} at the end of the run. It can be different from the input value for the same reason as for \code{dMax}.}
		\item{\code{log}}{The labels of the regressors removed thanks to the interactive interface during the run.}}
		
\note{The two first elemente of the returned list are lists that contain :
		\describe{
		\item{\code{def}}{ The output models structure. Each column corresponds to one selected model. Lines inform about the presence (value 1)
		of the absence (value 0) of the corresponding regressors, as in \code{toAnalyse}.}
		\item{\code{crit}}{A vector defining the models skill. The lower the value, the better 
		the models skill. If crit==NA the model was not tested by the function.}
		\item{\code{coeff}}{ The values of the regression coefficients. A matrix which has the same size as variable def and where 
		the value of the regression coefficient is given when def==1 whereas NA is given when def==0. More information about 
		these coefficients can be found in the lsfit documentation}
		}
}

\seealso{\link{labelPoly} \link{lsfit}}

\author{
Mangiarotti S., Coudret R., Drapeau L.
}

\examples{
	
	#Example 1
	data(NDVI,package="PoMoS")
	reg <- poMoS(as.matrix(NDVI)[,1:3],1000,20,20,10,nModStart=1000,dMax=3,show=0,filterToExpl=3)
        
        #Example 1
	data(rossler,package="PoMoS")
	reg <- poMoS(as.matrix(rossler),100,20,50,40,nModStart=100,dMax=2,show=0,filterToExpl=NULL)

}
