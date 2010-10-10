\name{rossler}
\alias{rossler}
\docType{data}

\title{Time series from Rossler (1976) attractor}

\description{
This dataset contains simulated variables of the differential equations of the Rossler (1976) attractor, that is: \cr
\eqn{dx/dt = - y - z} \cr
\eqn{dy/dt = x + a y} \cr
\eqn{dz/dt = b + z ( x - c )}
}

\usage{rossler}

\format{A data frame with 199999 observations on 3 variables. The columns correspond to \eqn{x}, \eqn{y} and \eqn{z}, respectively.}

\details{
The integration was made with a Runge-Kutta 4 algorithm. The following parameter values were used:
\describe{
	\item{\eqn{a}}{\eqn{0.52}}
	\item{\eqn{b}}{\eqn{2}}
	\item{\eqn{c}}{\eqn{4}}}

	Besides, the time between two points of the attractor is given by :
\describe{
	\item{\eqn{dt}}{\eqn{0.05}}}}

\author{Mangiarotti S., Coudret R., Drapeau L.}

\references{Rossler O. E. (1976). An equation for contiuous chaos. Physics Letters, 57A(5), pp. 397398.}

\examples{
	#Computation of the derivative of z
	data(rossler,package="PoMoS")
	N <- dim(rossler)[1]
	dt <- 0.05
	derivRossZ <- (rossler[2:N,3] - rossler[1:(N-1),3])/ dt 
	
	#Linear regression between dz/dt and 1, z, x z
	linRegResults <- lsfit(cbind(rossler[1:(N-1),3],rossler[1:(N-1),3]*rossler[1:(N-1),1]),
						   derivRossZ)
	linRegResults$coefficients
	#Our results were 2.1242690, -4.0691719, 0.9951174, which is very close to the parameters
	#2, -4 and 1 entered in the initial system
}

\keyword{datasets}