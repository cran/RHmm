\name{forwardbackward}
\alias{forwardbackward}
\title{forward-backward procedure}
\description{The forward-backward procedure is used to compute quantities used in the Baum-Welch algorithm.}
\usage{
forwardbackward(HMM, obs)
}
\arguments{
    \item{HMM}{a HMMClass or a HMMFitClass object}
    \item{obs}{a vector (matrix) of observations, or a list of vectors (or matrices) if there are more than one samples}
    }

\value{ If obs is one sample, a list of following elements, if obs is a list of samples, a
    list of list of following elements. See \bold{note} for mathematical definitions.   
    \item{Alpha}{The matrix of 'forward' probabilities (size: number of obs. times number of hidden states)}
    \item{Beta}{The matrix of 'backward' probabilities (size: number of obs. times number of hidden states)}
    \item{Gamma}{The matrix of probabilities of being at time t in state i (size: number of obs. times number of hidden states)}
    \item{Xsi}{The matrix of probabilities of being in state i at time t and being in state j at time t + 1 (size: number of obs. times number of hidden states)}
    \item{Rho}{The vector of probabilities of seeing the partial sequence obs[1] \ldots obs[t] (size number of obs.)}
    \item{LLH}{Log-likelihood}
     }

\note{
 Let \eqn{o=(o(1),\,\ldots,\,o(T))}{obs=(obs[1], \ldots obs[T])} be the 
 vector of observations, and \eqn{O=(O(t), t=1,\,\ldots,\,T)}{O=(O[t], 
 1, \ldots, T)}, the corresponding random variables. Let \eqn{Q=(Q(t), t=1,\,\ldots,\,T)}{(Q[t], t=1, \ldots, T)} 
 be the hidden Markov chain whose values are in \eqn{\left\{1,\,\ldots,\,nStates\right\}}{{1, \ldots, nStates}} 
 We have the 
 following definitions:\cr
 
 \eqn{\alpha_i(t) = 
 P(O_1=o(1),\,\ldots,\,O(t)=o(t),\,Q(t)=i\,|\,HMM)}{Alpha[i][t] = 
 P(O[1]=obs[1],\,\ldots,\,O[t]=obs[t],\,Q[t]=i | HMM)} which is 
 the probability of seeing the partial sequence 
 \eqn{o(1),\,\ldots,\,o(t)}{obs[1], \ldots, obs[t]} and ending up 
 in state i at time t.\cr
 
 \eqn{\beta_i(t) = P(O_{t+1}=o(t+1),\,\ldots,\,O(T)=o(T),\,Q(t)=i 
| HMM)}{Beta[i][t] = 
 P(O[t+1]=obs[t+1],\,\ldots,\,O[T]=obs[T],\,Q[t]=i | HMM)} which 
 is the probability of the ending partial sequence \eqn{o(t+1),\,\ldots,\,o(T)}{obs[t+1], \ldots, obs[T]} 
 given that we started at state i at time t.\cr
 
 \eqn{\Gamma_i(t) = P(Q(t) = i\,|\,O=o,\,HMM)}{Gamma[i][t] = P(Q[t]=i | O=obs, HMM)} which is the probability of being in state i 
 at time t for the state sequence \eqn{O=o}{O=obs}. \cr
 \eqn{\xi_i(t)=P(Q(t)=i,\,Q(t+1)=j\,|\,O=o,\,HMM)}{Xsi[i][t]=P(Q[t]=i, Q[t+1]=j | O=obs, HMM)} which is the probability of being 
 in state i at time t and being in state j at time t + 1.\cr
 
 \eqn{\rho(t) = P(O_1=o(1),\,\ldots,\,O_t=o(t)\,|\, HMM)}{Rho[t] = P(O[1]=obs[1], \ldots, O[t]=obs(t) | HMM)} witch is probabilities of seeing 
 the partial sequence \eqn{o(1),\,\ldots,\,o(t)}{obs[1] \ldots obs[t]}.\cr
 
 \eqn{LLH=\ln\rho[T]}{LLH=ln(Rho[T])}
}
\references{
    Jeff A. Bilmes (1997) \emph{ A Gentle Tutorial of the EM Algorithm and its Application to Parameter
    Estimation for Gaussian Mixture and Hidden Markov Models} \url{http://ssli.ee.washington.edu/people/bilmes/mypapers/em.ps.gz}
}
\examples{
    data(geyser)
    obs <- geyser$duration
    #Fits an 2 states gaussian model for geyser duration
    ResFitGeyser <- HMMFit(obs)
    #Forward-backward procedure
   fb <- forwardbackward(ResFitGeyser, obs)
}    

\keyword{htest}
