\name{data_mixture}
\docType{data}
\alias{data_mixture}
\title{Simulated univariate mixture of 3 gausssian distributions}
\description{This data set contains 10 samples with 100 observations of this HMM}
\note{Parameters for the simulation
Model:\cr
------\cr
2 states HMM with gaussian mixture distribution\cr
\cr
Initial probabilities:\cr
  Pi1 Pi2\cr
  0.4 0.6\cr

Transition matrix:\cr
        State 1 State 2\cr
State 1     0.8     0.2\cr
State 2     0.4     0.6\cr
\cr
Conditionnal distribution parameters:\cr
\cr
Distribution parameters:\cr
  State 1\cr
        mean var prop\cr
mixt. 1    1   2  0.1\cr
mixt. 2    2   3  0.2\cr
mixt. 3    3   4  0.7\cr
\cr
  State 2\cr
        mean var prop\cr
mixt. 1   -1   4  0.4\cr
mixt. 2   -2   2  0.3\cr
mixt. 3   -3   1  0.3\cr
}

\usage{data(data_mixture)}
\format{R script}

\keyword{datasets}
