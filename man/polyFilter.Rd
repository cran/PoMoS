\name{polyFilter}
\alias{polyFilter}
\title{Polynomial terms with the exponents of the variables}

\description{Generate a systematic correspondence of polynomial terms for a polynomial of degree \code{dMax} 
with \code{nVar} variables, and store it in a matrix form. Each column corresponds to one algebraic term; each line gives the exponent of the variables.}

\usage{polyFilter(nVar,dMax)}

\arguments{\item{nVar}{The number of variables for which polynomial terms are built}
	       \item{dMax}{The maximum degree of the polynomial terms. The same \code{dMax} value must be provided to \code{\link{poMoS}}}}

\value{ A matrix of degrees. Each column corresponds to one polynomial term. Each line represents the degrees of the variables. 
The polynomial term is the product each variable contribution: for example, the column c(1,2,0,2) leads to the polynomial term \eqn{X_1 X_2^2 X_4^2}}

\section{Note}{ This function can be used in complement to \code{labelPoly} to formally display the polynomial terms. 
Each expression produced by \code{labelPoly} corresponds to a vector of degree of polyFilter.}

\author{Mangiarotti S., Coudret R., Drapeau L.}

\seealso{\link{poMoS} \link{labelPoly}}

\examples{
	#Two different ways to represent the 8th polynomial series
	dMax <- 4
	nVar <- 3
	p <- polyFilter(nVar,dMax)
	p[,8]
	labelPoly(nVar,dMax,8)
}
