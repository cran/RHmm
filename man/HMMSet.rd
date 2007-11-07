\name{HMMSet}
\alias{HMMSet}
\title{Set the parameters for the hidden Markov models}

\description{This function is used to create a HMMClass object which contains the parameters of the HMM. An HMM is described by an initial state probability vector,
a transition matrix and a distributionClass object. }
\synopsis{
HMMSet(initProb, transMat, ...)
}
\usage{

HMMSet(initProb, transMat, distribution)
HMMSet(initProb, transMat, dis="NORMAL", mean, var)
HMMSet(initProb, transMat, dis="NORMAL", mean, cov)
HMMSet(initProb, transMat, dis="MIXTURE", mean, var, proportion)
HMMSet(initProb, transMat, dis="DISCRETE", proba, labels=NULL)
}
\arguments{
    \item{initProb}{the vector of probabilities of the initial state}
    \item{transMat}{the transition matrix of the hidden Markov chain}
    \item{distribution}{the distributionClass object of the observations}
    \item{dis}{dis parameter. See \bold{distributionSet}}
    \item{mean}{mean parameter. See \bold{distributionSet}}
    \item{var}{var parameter. See \bold{distributionSet}}
    \item{cov}{cov parameter. See \bold{distributionSet}}
    \item{proportion}{proportion parameter. See \bold{distributionSet}}
    \item{proba}{proba parameter. See \bold{distributionSet}}
    \item{labels}{labels parameter. See \bold{distributionSet}}
}
\value{ an object of class HMMClass
    \item{initProb}{initial state probabilities vector}
    \item{transMat}{transition matrix}
    \item{distribution}{distributionClass object}
}
\examples{
    # 3 hidden states Markov Model with univariate normal distributions
    # for the observations
    #   obs | hidden state = 1 are N(1, 1)
    #   obs | hidden state = 2 are N(-2, 2)
    #   obs | hidden state = 3 are N(5, 4)

        n_1d_3s <- distributionSet("NORMAL", c(1, -2, 5), c(1, 2, 4))
        initProb3 <- rep(1,3)/3
        transMat3 <- rbind(c(0.5, 0.4, 0.1), c(0.3, 0.4, 0.3),
            c(0.2, 0.1, 0.7))
        hmm1 <- HMMSet(initProb3, transMat3, n_1d_3s)
        # or directly
        hmm2 <- HMMSet(initProb3, transMat3, "NORMAL", mean=c(1, -2, 5),
            var=c(1, 2, 4))
 }
\seealso{\code{\link{distributionSet}}}
\keyword{models}
\keyword{distribution}
