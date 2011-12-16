\name{asymptoticSimCovMat}
\alias{asymptoticSimCovMat}
\title{Compute the asymptotic covariance matrix of a fitted HMM by simulation}
\description{This \sQuote{old} function computes the empirical asymptotic covariance matrix of the fitted HMM.}
\usage{
asymptoticSimCovMat(HMM, obs, nSimul, verbose=FALSE, oldCovMat=NULL, oldNSimul=0)
}
\arguments{
    \item{HMM}{a HMMClass or HMMFitClass object}
    \item{obs}{A vector, a matrix, a data frame, a list of vectors or a list of matrices of observations. See \code{\link{HMMFit}}.}
    \item{nSimul}{The number of simulation.}
    \item{verbose}{A boolean. if true, displays some informations. Default false.}
    \item{oldCovMat}{Last covariance matrix.}
    \item{oldNSimul}{Last number of simulations used to compute \sQuote{oldCovMat}. Useful to increase the size of the sample.}
}
\value{The empirical covariance matrix.}

\section{Numerical computations}{This is an ``experimental'' method. The HMM model is simulated nSimul times then fitted
and the empirical covariance matrix is computed.}
\examples{
    # Fit a 3 states 1D-gaussian model
    data(n1d_3s)
    Res <- HMMFit(obs_n1d_3s, nStates=3)
    # First 10 computations of covariance matrix
    Cov <- asymptoticSimCovMat(Res, obs_n1d_3s, 10)
    # 10 more computations of covariance matrix
    Cov <- asymptoticSimCovMat(Res, obs_n1d_3s, 10, verbose=TRUE, oldCovMat=Cov, oldNSimul=10)
    Res<-setAsymptoticCovMat(Res, Cov)
    summary(Res)
    }


\seealso{asymptoticCovMat, setAsymptoticCovMat}
