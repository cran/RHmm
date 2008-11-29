\name{asymptoticCovMat}
\alias{asymptoticCovMat}
\title{Asymptotic covariance matrix of the HMM parameters}
\description{This function calculates the empirical asymptotic covariance matrix of the HMM parameters}
\usage{
asymptoticCovMat(HMM, obs, asymptMethod=c("nlme", "optim"))
}
\arguments{
    \item{HMM}{a HMMClass or a HMMFitClass object}
    \item{obs}{The vector, matrix, data frame, list of vectors or list of matrices of observations}
    \item{asymptMethod}{A string which indicates the numerical method for computing the Hessian of parameters. Default 'nlme'.}
    
}
\value{A matrix}

\section{Numerical computations}{
    The Hessian of the LLH function is computed using finite difference approximations. 
    Either the stat package  'optimhess' internal function or the nlme package 'fdHess' function is used 
    for these computations. \cr
    There are a lot of numerical difficulties in computing derivatives in such models. 'optimhess' or 'fdHess' could return non inversible Hessian matrix.
}

\examples{
  data(n1d_3s)
  Res_n1d_3s<-HMMFit(obs_n1d_3s, nStates=3)
  covMat <- asymptoticCovMat(Res_n1d_3s, obs_n1d_3s, asymptMethod='optim')
}

\references{    
    Visser Ingmar, Raijmakers Maartje E. J. and  Molenaar Peter C. M.(2000) \emph{Confidence intervals for hidden Markov
    model parameters}, British Journal of Mathematical and Statistical Psychology, 53, 317-327.
    
    Mann Tobias P. (2006) \emph{Numerically Stable Hidden Markov Model Implementation}, \url{http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf}
}

\seealso{HMMFit}
