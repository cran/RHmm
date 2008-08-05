\name{HMMGraphicDiag}
\alias{HMMGraphicDiag}
\title{Graphic diagnostic of the HMM estimation}
\description{This function plots the kernel density of the observations and the normal (mixture of normal, discrete) density
    with estimated parameters for each hidden states. The hidden states}
\usage{
HMMGraphicDiag(vit, HMM, obs, color="green")
}
\arguments{
    \item{vit}{a ViterbiClass object which gives the hidden states}
    \item{HMM}{a HMMClass or a HMMFitClass object which describes the model}
    \item{obs}{the vector, list of vectors of observations}
    \item{color}{color for the kernel density plot}
}

\value{none}
\note{
 HMMGraphicDiag is not implemented for multivariate distributions.\cr

 The kernel densities of observations for each hidden states of the model are plotting using:\cr
    plot(density(obs[vit$states=i])) and i in 1..HMM$nStates (or HMM$HMM$nStates)

  }
\examples{
data(geyser)
obs <- as.matrix(geyser)
#Fits an 3 states gaussian model for geyser duration
ResFitGeyser <- HMMFit(obs, dis='MIXTURE', nStates=3, nMixt=2)
VitGeyser <- viterbi(ResFitGeyser, obs)
m1Duration <- c(ResFitGeyser$HMM$distribution$mean[[1]][[1]][1], ResFitGeyser$HMM$distribution$mean[[1]][[2]][1])
m2Duration <- c(ResFitGeyser$HMM$distribution$mean[[2]][[1]][1], ResFitGeyser$HMM$distribution$mean[[2]][[2]][1])
m3Duration <- c(ResFitGeyser$HMM$distribution$mean[[3]][[1]][1], ResFitGeyser$HMM$distribution$mean[[3]][[2]][1])
v1Duration <- c(ResFitGeyser$HMM$distribution$cov[[1]][[1]][1,1], ResFitGeyser$HMM$distribution$cov[[1]][[2]][1,1])
v2Duration <- c(ResFitGeyser$HMM$distribution$cov[[2]][[1]][1,1], ResFitGeyser$HMM$distribution$cov[[2]][[2]][1,1])
v3Duration <- c(ResFitGeyser$HMM$distribution$cov[[3]][[1]][1,1], ResFitGeyser$HMM$distribution$cov[[3]][[2]][1,1])
prop <- list(ResFitGeyser$HMM$distribution$proportion[[1]], ResFitGeyser$HMM$distribution$proportion[[2]], 
    ResFitGeyser$HMM$distribution$proportion[[3]])
prop[[1]] <- prop[[1]] / sum(prop[[1]])
prop[[2]] <- prop[[2]] / sum(prop[[2]])

HMMDuration <- HMMSet(dis='MIXTURE', transMat=ResFitGeyser$HMM$transMat, initProb=ResFitGeyser$HMM$initProb, 
    mean=list(m1Duration, m2Duration, m3Duration), var=list(v1Duration, v2Duration, v3Duration), proportion=prop) 
# Graphic diagnostic
HMMGraphicDiag(VitGeyser, HMMDuration, obs[,1])
}
\seealso{HMMFit, viterbi}
\keyword{htest}
