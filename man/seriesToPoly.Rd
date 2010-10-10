\name{seriesToPoly}
\alias{seriesToPoly}
\title{Polynomial series construction}

\description{Creates time series corresponding to polynomial terms by multiplying time series among them. One time series can be multiplied several times.}

\usage{seriesToPoly(series, dMax=NULL, pFilter=NULL)}

\arguments{
\item{series}{A matrix containing the time series from which the polynomial series are built. One time series by column.}
\item{dMax}{The maximum degree of the polynomial. The same value must be used in \code{\link{poMoS}}}
\item{pFilter}{A matrix filled representing one polynomial term in each column and the power of each variable 
		in each line (same structure as the output of \code{polyFilter}). If this parameter is not set, \code{polyFilter} is called, instead. In this case, 
		\code{dMax} is required, whereas \code{nVar} is deduced from the size of series. }
}

\value{A matrix which have as number of lines as series and which comprehend the polynomial series measures based on it. 
		Each column is relative to one series. To know the expression of a series, use \code{\link{labelPoly}}.}

\seealso{\link{polyFilter} \link{labelPoly}}

\author{Mangiarotti S., Coudret R., Drapeau L.}

\examples{
	#Both examples give the first values of series 'X_1 X_2^3 X_3'
	data(sprottK,package="PoMoS")
	sprottK <- as.matrix(sprottK)
	dMax <- 5
	nVar <- dim(sprottK)[2]
	
	#Example 1
	polySeries <- seriesToPoly(sprottK,dMax)
	lab <- labelPoly(nVar,dMax)
	polySeries[1:100,lab=="X1 X2^3 X3 "]
	
	#Example 2
	p <- c(1,3,1)
	polySeries <- seriesToPoly(sprottK,pFilter=p)
	polySeries[1:100]
}
