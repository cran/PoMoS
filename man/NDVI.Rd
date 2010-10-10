\name{NDVI}
\alias{NDVI}
\docType{data}

\title{Normalized Difference Vegetation Index (NDVI) time series}

\description{Time series of Normalized Difference Vegetation Index (NDVI) derived from the AVHRR (Advanced Very High Resolution Radiometer) 
sensors of the NOAA (National Oceanic and Atmospheric Administration) space missions. 
These values have been spatially averaged over the province of Settat, Morocco (see hereafter for details).}

\usage{NDVI}

\format{A data frame of 4 variables with 9618 interpolated values each. The first column represents the vegetation index while the others are the successive derivatives.}

\details{NDVI data from AVHRR sensor are available since 1981. The present time series is derived from the GIMMS product of the 
Global Land Cover Facility. This product is 10-day sampling with an 8x8km2 resolution. The present time series has been spatially 
averaged over the province of Settat (Morocco), representing 150 pixels of 8x8km2. The time series has been smoothed and resampled 
at a higher time resolution (roughly daily sampling) with a Savitzky-Golay method.
\describe{\item{\code{dt}}{1/365}}
}


\author{Mangiarotti S., Coudret R., Drapeau L.}

\source{Global Land Cover Facility\cr
University of Maryland, Department of Geography\cr
\url{http://www.landcover.org}
}

\references{
Pinzon, J., Brown, M.E. and Tucker, C.J., 2005.\cr
 Satellite time series correction of orbital drift artifacts using empirical mode decomposition.\cr
 In: N. Huang (Editor), Hilbert-Huang Transform: Introduction and Applications, pp. 167-186. 

Tucker C.J., Pinzon J.E., Brown M.E., Slayback D.A., Pak E.W., Mahoney R., Vermote E.F. and Saleous N.E., 2005. \cr
An extended AVHRR 8-km NDVI dataset compatible with MODIS and SPOT vegetation NDVI data.  \cr
International Journal of Remote Sensing, 26:20, 4485-4498. http://dx.doi.org/10.1080/01431160500168686
}

\examples{
#2D plots like those drawn by the gloMo function of the GloMo package
data(NDVI,package="PoMoS")
par(mfrow=c(2,2))
for (i in 1:2) {
	for (j in (i+1):3) {
		plot(NDVI[,i],NDVI[,j],type="l",main="NDVI Series",
			xlab=paste("X",i,sep=""),ylab=paste("X",j,sep=""))
	}
}


}

\keyword{datasets}
