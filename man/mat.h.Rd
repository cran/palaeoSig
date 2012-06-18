\name{mat.h}
\alias{mat.h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MAT with minimum geographic distance for analogues }
\description{
h-block cross-validation for MAT. 
}
\usage{
mat.h(y, x, noanalogues, geodist, thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Community data }
  \item{x}{ Environmental variable }
  \item{noanalogues}{ Number of analogues }
  \item{geodist}{ Geographical distances }
  \item{thresh}{ Threshold distance below which sites are excluded during cross validation }
}
\details{
Analogues closer than the geographic distance specified by thresh are not excluded, reducing the effect of spatial autocorrelation on the performance statistics. 
  
}
\value{
  Predictions for each site.
}
\references{ Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. DOI: 10.1016/j.quascirev.2008.12.020 }
\author{ Richard Telford \email{Richard.Telford@bio.uib.no} }

\examples{

require(mapproj)
require(gstat)
data(arctic.env)
data(arctic.pollen)
try({
sinproj <- mapproject(x = arctic.env$Longitude, y = arctic.env$Latitude, proj = "sinusoidal", ori = rotate(arctic.env$Lat, arctic.env$Long))
radius <- 6372
arctic.pos <- data.frame(x=sinproj$x*radius, y=sinproj$y*radius)
arctic.dist <- dist(arctic.pos)

recon <- wajack(sqrt(arctic.pollen), arctic.env$tjul)
recon.resid <- recon-arctic.env$tjul
recon.resid.detrend <- resid(loess(recon.resid ~ arctic.env$tjul, span = .4))

resid.variogram <- variogram(recon.resid.detrend ~ 1, ~x+y, data = arctic.pos)
(resid.vModel <- fit.variogram(resid.variogram, model = vgm(2, "Sph", 300, 1)))
plot(resid.variogram, resid.vModel)

mat.recon <- mat.h(arctic.pollen, arctic.env$tjul, noanalogues = 10, arctic.dist, -1)
mat.reconH <- mat.h(arctic.pollen, arctic.env$tjul, noanalogues = 10, arctic.dist, resid.vModel$range[2])

cor(mat.recon, arctic.env$tjul)
cor(mat.reconH, arctic.env$tjul)

sqrt(mean((mat.recon - arctic.env$tjul)^2))
sqrt(mean((mat.reconH - arctic.env$tjul)^2))
})

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
