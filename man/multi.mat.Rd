\name{multi.mat}
\alias{multi.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MAT for multiple variables }
\description{
  MAT for many environmental variables simultaneously. More efficient than calculating them seperately for each variable.
}
\usage{
multi.mat(training.spp, envs, core.spp, noanalogues = 10, method = "sq-chord", run = "both")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training.spp}{ Community data }
  \item{envs}{ Environmental variables - or simulations }
  \item{core.spp}{ Optional fossil data to make predictions for}
  \item{noanalogues}{ Number of analogues to use }
  \item{method}{ distance metric to use }
  \item{run}{ Return LOO predictions or predictions for fossil data  }
}
\value{
  Matrix of predictions
}
\references{ Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. DOI: 10.1016/j.quascirev.2008.12.020  }
\author{ Richard Telford \email{Richard.Telford@bio.uib.no}  }
\examples{

require(mapproj)
require(gstat)
data(arctic.env)
data(arctic.pollen)
try({
sinproj <- mapproject(x = arctic.env$Longitude, y = arctic.env$Latitude, proj = "sinusoidal", ori = rotate(arctic.env$Lat, arctic.env$Long))
radius <- 6372
arctic.data <- data.frame(x = sinproj$x*radius, y = sinproj$y*radius, tjul = arctic.env$tjul)

arctic.variogram <- variogram(tjul~1, ~x+y, data = arctic.data, cutoff = 4000)
(arctic.vModel <- fit.variogram(arctic.variogram, model = vgm(20, "Sph", 1500)))
plot(arctic.variogram, arctic.vModel)

arctic.data$x <- jitter(arctic.data$x)
sims <- simulate(arctic.data, nsim=99, arctic.pollen, arctic.vModel)

sims <- cbind(arctic.env$tjul, sims)

sim.MAT.R2 <- mapply(cor, as.data.frame(multi.mat(arctic.pollen, sims, noanalogues = 5)[[1]]), as.data.frame(sims))^2
mean(sim.MAT.R2[1] <= sim.MAT.R2)
hist(sim.MAT.R2)
abline(v = sim.MAT.R2[1])
})

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
