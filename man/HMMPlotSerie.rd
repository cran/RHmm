\name{HMMPlotSerie}
\alias{HMMPlotSerie}
\title{Plot univariates series in each estimated states}
\description{This function plots the time series in each hidden state.}
\usage{
    HMMPlotSerie(obs, states, dis="NORMAL", color="green")
    }
\arguments{
    \item{obs}{the vector, list of vectors of observations}
    \item{states}{a ViterbiClass object which gives the hidden states or a vector or a lis of vectors of integer 1 to the number of hidden states}
    \item{dis}{Distribution name = 'NORMAL', 'DISCRETE', 'MIXTURE'. Default 'NORMAL'.}
    \item{color}{color for the kernel density plot}
}

\value{none}
\note{
 HMMPlotSerie is not implemented for multivariate observations.\cr

 The time series of observations for each hidden states of the model are plotted using:\cr
    plot(obs[states=i])) and i in 1..max(States)

  }
\examples{
data(geyser)
obs <- geyser$duration
#Fits an 3 states gaussian model for geyser duration
ResFitGeyser <- HMMFit(obs, nStates=3)
VitGeyser <- viterbi(ResFitGeyser, obs)
#plot the series
HMMPlotSerie(obs, VitGeyser)
}
\seealso{viterbi, plot}
\keyword{htest}
