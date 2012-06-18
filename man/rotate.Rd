\name{rotate}
\alias{rotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rotate}
\description{
  Provides the rotation required to orientate the dataset across the equator at zero longitude. This makes the map projection give a good correlation between projected distances and great circle distances. Useful for local to hemispheric datasets - not global ones.
}
\usage{
rotate(lats, longs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lats}{ Latitudes of sites}
  \item{longs}{ Longitudes of sites }
}

\value{
  Orientation as required for mapproject
}

\author{Richard Telford \email{Richard.Telford@bio.uib.no} }

\seealso{  \code{\link[mapproj]{mapproject}} }
\examples{
require(mapproj)
data(arctic.env)
map(xlim = c(-180, -15), ylim = c(45, 90), projection = "sinusoidal") #default orientation
try({
points(mapproject(x = arctic.env$Longitude, y = arctic.env$Latitude), pch = 20, col = 2)
})
map.grid()

ori<-rotate(arctic.env$Lat, arctic.env$Long)

map(xlim = c(-180, -15), ylim = c(45, 90), projection = "sinusoidal", orientation = ori)
try({
points(mapproject(x = arctic.env$Longitude, y = arctic.env$Latitude, proj = "sinusoidal", orientation = ori), pch = 20, col = 2)
map.grid()
})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

