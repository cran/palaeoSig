\name{simulateSpatial}
\alias{simulateSpatial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Unconditional Gaussian Simulation }
\description{
  A wrapper for krige
}
\usage{
simulateSpatial(geo.xy, nsim = 999, y, varmod, centre = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geo.xy}{ data.frame with columns x & y }
  \item{nsim}{Number of simulations}
  \item{y}{The environmental variable being simulated (required just to get mean)}
  \item{varmod}{ Variogram model}
  \item{centre}{ Should the simulations be centred (useful if using hermite transformations with large range) }
}
\details{
Simulates environmental variables with same spatial structure as the observed environmental variable.  
}
\value{
  Data.frame of simulations.
}
\references{ Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. DOI: 10.1016/j.quascirev.2008.12.020   }
\author{ Richard Telford \email{Richard.Telford@bio.uib.no}  }
\seealso{ \code{\link[gstat]{krige}}, \code{\link[gstat]{fit.variogram}} }
\examples{

  require(mapproj)
  require(gstat)
  data(arctic.env)
  data(arctic.pollen)

  sinproj <- mapproject(x = arctic.env$Longitude, y = arctic.env$Latitude, proj = "sinusoidal", ori = rotate(arctic.env$Lat, arctic.env$Long))
  radius <- 6372
  arctic.data <- data.frame(x = sinproj$x*radius, y = sinproj$y*radius, tjul = arctic.env$tjul)
  
  arctic.variogram <- variogram(tjul~1, ~x+y, data = arctic.data, cutoff = 4000)
  (arctic.vModel <- fit.variogram(arctic.variogram, model = vgm(20, "Sph", 1500)))
  plot(arctic.variogram, arctic.vModel)
  
  arctic.data$x <- jitter(arctic.data$x)

  \dontshow{
    try({
    arctic.data2<-arctic.data[1:100,]#just use first 100 observation for speed
    arctic.pollen2<-arctic.pollen[1:100,]
    
    sims <- simulateSpatial(arctic.data2, nsim=99, arctic.data2$tjul, arctic.vModel) #99 simulations is way too few except as a demonstration
    
    sims <- cbind(arctic.data2$tjul, sims)
    
    sim.MAT.R2 <- mapply(cor, as.data.frame(multi.mat(arctic.pollen2, sims,noanalogues = 5)[[1]]), as.data.frame(sims))^2
    mean(sim.MAT.R2[1] <= sim.MAT.R2)
    hist(sim.MAT.R2)
    abline(v = sim.MAT.R2[1])
    })
  }
  \donttest{
    try({ #slow   
    sims <- simulateSpatial(arctic.data, nsim=99, arctic.data$tjul, arctic.vModel) #99 simulations is way too few except as a demonstration
    
    sims <- cbind(arctic.data$tjul, sims)
    
    sim.MAT.R2 <- mapply(cor, as.data.frame(multi.mat(arctic.pollen, sims,noanalogues = 5)[[1]]), as.data.frame(sims))^2
    mean(sim.MAT.R2[1] <= sim.MAT.R2)
    hist(sim.MAT.R2)
    abline(v = sim.MAT.R2[1])
    })
  }



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
