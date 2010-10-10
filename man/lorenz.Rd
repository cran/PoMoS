\name{lorenz}
\alias{lorenz}
\docType{data}

\title{Time series from Lorenz (1963) attractor}

\description{
This dataset contains simulated variables of the differential equations of the Lorenz (1963) attractor, that is:\cr
\eqn{dx/dt = s( y - x ) } \cr
\eqn{dy/dt = x(r - z) - y} \cr
\eqn{dz/dt = x y - b z}
}

\usage{lorenz}

\format{A data frame with 10000 observations on 3 variables. The columns correspond to \eqn{x}, \eqn{y} and \eqn{z}, respectively.}

\details{The integration was made thanks to the ordinary differential equations solver of \code{deSolve} package. The following parameter values were used:\cr 
\describe{
	\item{\eqn{s}}{\eqn{10}}
	\item{\eqn{r}}{\eqn{28}}
	\item{\eqn{b}}{\eqn{8/3}}
}
Besides, the time between two points of the attractor is given by :
\describe{
	\item{\eqn{dt}}{\eqn{0.005}}
}
}

\author{Mangiarotti S., Coudret R., Drapeau L.}

\references{Lorenz E.N. (1963). Deterministic Nonperiodic Flow. Journal of the atmospheric sciences, 20, pp. 130141.}

\examples{
	#Selecting a single side of the Lorenz attractor
	data(lorenz,package="PoMoS")
	k <- kmeans(lorenz,2)
	plot(lorenz[k$cluster==1,1],lorenz[k$cluster==1,2])}

\keyword{datasets}