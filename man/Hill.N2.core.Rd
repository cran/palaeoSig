\name{Hill.N2.core}
\alias{Hill.N2.core}

\title{
Calculate the effective number of species in the fossil data
}
\description{
Gives a measure of the species divesity in the fossil data.
}
\usage{
Hill.N2.core(spp)
}
\arguments{
  \item{spp}{
Species data
}
}
\details{
Uses \code{\link{Hill.N2}} from the rioja package 

}
\value{
Effective number of species
}
\references{

}
\author{
Richard Telford
}
\note{
If the effective number of species is small, WA based reconstructions are unlikely to be significant, and MAT based reconstructions should be tested instead.
}


\seealso{
\code{\link{Hill.N2}}
}
\examples{
data(RLGH)

Hill.N2.core(RLGH$spec)
}
\keyword{ attribute }

