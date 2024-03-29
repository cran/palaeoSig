% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.fittedAgelme.R
\name{plot.fittedAgelme}
\alias{plot.fittedAgelme}
\title{Plots fitted agelme model and dates}
\usage{
\method{plot}{fittedAgelme}(x, main, xlab = "Depth", ylab = "Calibrated Age", ...)
}
\arguments{
\item{x}{Fitted agelme model.}

\item{main}{Title of the plot.}

\item{xlab}{x axis label of the plot.}

\item{ylab}{y axis label of the plot.}

\item{\dots}{Other arguments to plot.}
}
\description{
Plots fitted agelme model and dates
}
\examples{
data(STOR)
fit.mod <- with(STOR, agelme(depthup, depthdo, cageup, cagedo))
#Predicting using the constant variance model,
#for each cm between 70 and 400 cm.
fit.pre <- predict(fit.mod, 1, 70:400)
plot(fit.pre)
}
