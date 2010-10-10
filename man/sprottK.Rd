\name{sprottK}
\alias{sprottK}
\docType{data}

\title{Time series from Sprott K (1994) attractor}

\description{This dataset contains simulated variables of the differential equations of the SprottK attractor, that is:  : \cr
\eqn{dx/dt = -a y + x z} \cr
\eqn{dy/dt = x + b y} \cr
\eqn{dz/dt = x - z}
}

\usage{sprottK}

\format{A data frame with 199999 observations on 3 variables. The columns correspond to \eqn{x}, \eqn{y} and \eqn{z}, respectively.}

\details{
The integration was made with a Runge-Kutta 4 algorithm. The following parameter values were used: 
\describe{
	\item{\eqn{a}}{\eqn{0.432}}
	\item{\eqn{b}}{\eqn{0.2}}
}
Besides, the time step between two points of the attractor is given by:
\describe{
	\item{\eqn{dt}}{\eqn{0.05}}}}

\author{Mangiarotti S., Coudret R., Drapeau L. }
\references{Sprott J.C., 1994. Some simple chaotic flows. Physical Review E, 50(2), pp. 647-650.}

\examples{
	#Computation of the derivative of x
	data(sprottK,package="PoMoS")
	N <- dim(sprottK)[1]
	dt <- 0.05
	derivSproX <- (sprottK[2:N,1] - sprottK[1:(N-1),1])/ dt 
	
	#Linear regression between dx/dt and a model with regressors 1, x and y
	#which does not corresponds to the initial equation
	linRegResults <- lsfit(cbind(sprottK[1:(N-1),1],sprottK[1:(N-1),2]),derivSproX)
	normalizedRes <- (linRegResults$residuals - mean(linRegResults$residuals))/sqrt(var(linRegResults$residuals)) 
	ks.test(normalizedRes,"pnorm")
	#We reject the normality test, so we cannot use criterion like AIC to find
	#the good model, for these series.
}

\keyword{datasets}