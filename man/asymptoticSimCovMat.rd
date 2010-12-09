\name{asymptoticSimCovMat}
\alias{asymptoticSimCovMat}
\title{Compute the asymptotic covariance matrix of a fitted HMM by simulation}
\description{This function compute the empirical asymptotic covariance matrix of the fitted HMM.}
\usage{
asymptoticSimCovMat(HMM, obs, nSimul, verbose)
}
\arguments{
    \item{HMM}{a HMMClass or HMMFitClass object}
    \item{obs}{A vector, a matrix, a data frame, a list of vectors or a list of matrices of observations. See 
           HMMFit.}
    \item{nSimul}{The number of simulation}
    \item{verbose}{A boolean. if true, displays some informations. Default false.}

}
\value{the empirical covariance matrix}

\section{Numerical computations}{This is an ``experimental'' method. The HMM model is simulated nSimul times then fitted
and the emirical covariance matrix is computed.}



\seealso{asymptoticCovMat, setAsymptoticCovMat}
