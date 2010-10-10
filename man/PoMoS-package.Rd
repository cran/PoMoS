\name{PoMoS}
\alias{PoMoS}
\docType{package}
\title{Polynomial (ordinary differential equation) Model Search}
\description{This package aims to determine from a set of N time series the optimal polynomial structure of a model built on first-order ordinary differential equations. The core of the package is based on the poMoS function: an evolutionary algorithm combined with a least square fitting. Optimality is estimated with AIC (Akaike, 1974) or AIC-like criterions. Although efficient in its selection, the identification of the optimal structure cannot be guaranteed. Therefore, both selected and rejected models are reconsidered after optimal solutions are obtained from the evolutionary algorithm for another analysis. This analysis is based on a statistical evaluation of the regressors quality.
An important level of handling can be kept during the model search process: regressors can thus be added and removed dynamically, combining both the evolutionary model search with the statistical analysis of the regressors. Parameters of the evolutionary algorithm can also be modified during the searching process etc. A blind run of the algorithm is also made technically possible.
}
\details{
\tabular{ll}{
Package: \tab PoMoS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-10-10\cr
Depends: \tab R (>= 2.11.1), RGtk2, cairoDevice, igraph\cr
SystemRequirements: \tab Gtk+ (>= 2.16.0)\cr
License: \tab CeCILL-2\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
bpRegress               Boxplot representation of models criterion
convGraph              	Graphical plots of neighboring models
findRegData             Information searching from poMoS output
labelPoly               Polynomial labels systematic ordering in a set of N variables
polyFilter              Polynomial terms with the exponent of the variables
poMoS                   Polynomial Model Search
seriesToPoly            Polynomial series construction
}

The main function is \code{PoMoS} from which major part of the other functions are called. Gtk+ and RGtk2 provide a link between \code{PoMoS}  and R and allows for an interactive monitoring of the algorithm search process.
Several criterions can be chosen to evaluate the models skill, including Akaike's and the BIC (Bayesian Information Criterion). These criterions are based on the sum of two terms, one term corresponding to the likelihood, another term corresponding to a penalty on larger length models.
\code{bpRegress} and \code{convGraph} are representation functions. \code{findRegData} helps in browsing the results obtained from \code{poMoS} algorithm. \code{labelPoly}, \code{polyFilter} and \code{seriesToPoly} manage the polynomial series and their incorporation into models. These functions are also required in the GloMo-package which, therefore, depends on PoMoS-package.

}

\author{
Mangiarotti S., Coudret R., Drapeau L.

Maintainer: Drapeau L., : laurent.drapeau@ird.fr
}

\references{
Akaike, Hirotugu (1974). "A new look at the statistical model identification". IEEE Transactions on Automatic Control 19 (6): 716-723.
}

\keyword{package}
\keyword{datasets}
\keyword{regression}
\keyword{ts}

\seealso{\code{\link{poMoS}}}
