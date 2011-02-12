\name{obs.cor}
\alias{obs.cor}
\alias{sim.cor}
\alias{plot.obscor}
\alias{plot.simcor}
\alias{obscor.sig}

\title{
Weighted correlation between weighted averaging optima and constrained ordination species scores.
}
\description{
obs.cor calculates the weighted correlation between the species weighted average optima and the axis one species scores of an ordination constrained by the WA reconstruction. 
sim.cor calculates a null distribution of this statistic.
plot.simcor plots a histogram showing the null distribution and a vertical line showing the observed weighted correlation.
plot.obscor plots the species weighted average optima against the axis 1 constrained ordination scores
obscor.sig calculates the p value of the observed correlation.
}
\usage{
obs.cor(spp, env, fos, ord = rda)
sim.cor(spp, env, fos, n = 99, ord = rda) 
plot.obscor(x, xlab = "WA optima", ylab = "RDA scores", f = 1, ...) 
plot.simcor(x, sim,...) 
obscor.sig(x, sim) 
}


\arguments{
  \item{spp}{
Data frame of modern training set species data, transformed if required, for example with sqrt  
}
  \item{env}{
Vector of a single environmental variable
}
  \item{fos}{
Data frame of fossil species data. Species codes and tranformations should match those in spp.
}
  \item{ord}{
Constrained ordination. The default is rda. cca will work with the present code. capscale would require the code to be altered.
}
  \item{n}{
Number of random trials
}
  \item{x}{
Output of obs.cor
}
  \item{xlab}{
X-axis label
}
  \item{ylab}{
Y-axis label
}
  \item{\dots}{
Other arguments to plot
}
  \item{f}{
Scale factor for the points
}
  \item{sim}{
Output of sim.cor
}

}
\details{
A square root transformation of the species data is often useful. n=99 is too small in practice to give a smooth histogram of the null model. n=999 is better.
}
\value{

obs.cor returns a list
 \item{x}{Data.frame with columns Optima, RDA1 and abun containing the species optima, ordination axis 1 scores, and mean abundance respectively}
 \item{res}{List with components wc and cc containing the weighted correlation and unweighted correlation between species optima and ordination scores respectively}

sim.cor returns a matrix with the correlation between species weighted average optima and ordination scores in the first column and the weighted correlation in the second column. Each row represents a different random environmental variable.

obscor.sig returns the p-value for the observed correlation between species weighted average optima and ordination scores.
}

\references{
Telford, R.J. and Birks, H.J.B. (In Revision) A novel method for assessing the statistical significance of quantitative reconstructions inferred from biotic assemblages. Quaternary Science Reviews.
}
\author{
Richard Telford \email{richard.telford@bio.uib.no}
}
\note{
The test of the weighted correlation between species optima and ordination axis scores is more powerful, especially with a small number of fossil observations, that the test of variance explained in  \code{\link{randomTF}} but is only applicable to WA and will have a large type II error if there are few species. 
}

\seealso{
\code{\link{randomTF}}\code{\link{WA}},\code{\link{MAT}},\code{\link{WAPLS}},\code{\link{rda}},\code{\link{cca}}
}
\examples{
data(SWAP)
data(RLGH)


rlgh.obs<-obs.cor(spp=sqrt(SWAP$spec),env=SWAP$pH,fos=sqrt(RLGH$spec))
rlgh.sim<-sim.cor(spp=sqrt(SWAP$spec),env=SWAP$pH,fos=sqrt(RLGH$spec))  #slow

plot.obscor(rlgh.obs)
plot.simcor(rlgh.obs,rlgh.sim)

obscor.sig(rlgh.obs,rlgh.sim)
}

\keyword{multivariate}
\keyword{htest}
\keyword{hplot}