 ###############################################################
 #### RHmm version 1.4.3                               
 ####                                                         
 #### File: RHmm-HyperNew.R 
 ####                                                         
 #### Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr>
 #### Author: Sebastian BAUER <sebastian.bauer@charite.de>
 #### Date: 2010/12/01                                      
 ####                                                         
 ###############################################################

tolMin <- .Machine$double.eps*100
# Contraintes sur les paramètres
eProba <- 1
eVar <- 2
eCor <- 3
eNone <- 0

sumList <- function(List, n)
{
    if (is.list(List))
    {   Res <- 0
        for (i in 1:n)
            Res <- Res + List[[i]]
    }
    else
        Res <- List
    return(Res)
}

GetNParam<-function(object) UseMethod("GetNParam")
GetNAllParam<-function(object) UseMethod("GetNAllParam")
GetVectParam <- function(object) UseMethod("GetVectParam")
GetVectAllParam <- function(object) UseMethod("GetVectAllParam")
SetVectParam <- function(object, Vect) UseMethod("SetVectParam")
SetVectAllParam <- function(object, Vect) UseMethod("SetVectAllParam")
GetNConstraint <- function(object) UseMethod("GetNConstraint")
GradConstraint <- function(object) UseMethod("GradConstraint")

GetNParam.default<-function(object)
{
    return(list(nParam=0,paramConstr=NULL))
}

GetNAllParam.default<-function(object)
{
    return(list(nParam=0, paramConstr=NULL))
}

GetVectParam.default<-function(object)
{
    return(NULL)
}

GetVectAllParam.default<-function(object)
{
    return(NULL)
}

SetVectParam.default<-function(object, Vect)
{
    return(NULL)
}


SetVectAllParam.default<-function(object, Vect)
{
    return(NULL)
}

GetNConstraint.default<-function(object)
{
    return(0)
}

GradConstraint.default<-function(object)
{
    return(NULL)
}


GetNParam.HMMFitClass<-function(object)
{
    return(GetNParam(object$HMM))
}

GetNAllParam.HMMFitClass<-function(object)
{
    return(GetNAllParam(object$HMM))
}

GetVectParam.HMMFitClass<-function(object)
{
    return(GetVectParam(object$HMM))
}

GetVectAllParam.HMMFitClass<-function(object)
{
    return(GetVectAllParam(object$HMM))
}

SetVectParam.HMMFitClass<-function(object, Vect)
{
    return(SetVectParam(object=object$HMM, Vect=Vect))
}

SetVectAllParam.HMMFitClass<-function(object, Vect)
{
    return(SetVectAllParam(object=object$HMM, Vect=Vect))
}

GetNConstraint.HMMFitClass<-function(object)
{
    return(GetNConstraint(object$HMM))
}

GradConstraint.HMMFitClass<-function(object)
{
    return(GradConstraint(object$HMM))
}


GetNParam.HMMClass<-function(object)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nStates <- object$distribution$nStates
    nOtherParam <- (nStates - 1)*(nStates + 1)
    Res1<-GetNParam(object$distribution)
    nParam <- Res1$nParam * nStates + nOtherParam
    paramConstr <- c(rep(eProba, nOtherParam), rep(Res1$paramConstr, nStates))
    Res <- list(nParam=nParam, paramConstr=paramConstr)
    return(Res)
}


GetNAllParam.HMMClass<-function(object)
{   nStates <- object$distribution$nStates
    Res1<-GetNAllParam(object$distribution)
    nOtherParam <- nStates*(nStates+1)
    nParam <- Res1$nParam * nStates + nOtherParam
    paramConstr <- c(rep(eProba, nOtherParam), rep(Res1$paramConstr, nStates))
    Res <- list(nParam=nParam, paramConstr=paramConstr)
    return(Res)
}

GetVectParam.HMMClass<-function(object)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nParam <- GetNParam(object)$nParam
    nStates <- object$distribution$nStates
    Res <- rep(0, nParam)
    #initProb
    for (i in 1:(nStates-1))
    {   Res[i] <- object$initProb[i]
    }
    k <- nStates
    for (i in 1:nStates)
    {   for (j in 1:(nStates -1))
        {   Res[k] <- object$transMat[i,j]
            k <- k + 1
        }
    }
    Res[k:nParam] <- GetVectParam(object$distribution)
    return(Res)
}

GetVectAllParam.HMMClass<-function(object)
{
    nStates <- object$distribution$nStates
    Res <- object$initProb
    
    for (i in 1:nStates)
        Res <- c(Res, object$transMat[i,])
    Res1 <- GetVectAllParam(object$distribution)
    return(c(Res, Res1))
}

SetVectParam.HMMClass<-function(object, Vect)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nStates <- object$distribution$nStates
    
    #initProb
    initProb <- rep(0, nStates)
    initProb[1:(nStates-1)] <- Vect[1:(nStates-1)]
    initProb[nStates] <- 1.0 - sum(initProb)
    
    #transMat
    transMat <- matrix(0, nrow=nStates, ncol=nStates)
    k <- nStates
    for (i in 1:nStates)
    {   transMat[i,1:(nStates-1)] <- Vect[k:(k+nStates-2)]
        k <- k + nStates-1
        transMat[i,nStates] <- 1.0 - sum(transMat[i,])
    }
    
    ResDistr <- SetVectParam(object=object$distribution, Vect=Vect[k:length(Vect)])
    HMM<-HMMSet(initProb=initProb, transMat=transMat, ResDistr)
    return(HMM)
}

SetVectAllParam.HMMClass<-function(object, Vect)
# nStates - 1 probInit + nStates*(nStates-1) matTrans + GetNParam de la distribution
{
    nStates <- object$distribution$nStates
    
    #initProb
    initProb <- Vect[1:nStates]
    
    #transMat
    transMat <- matrix(0, nrow=nStates, ncol=nStates)
    k <- nStates+1
    for (i in 1:nStates)
    {   transMat[i,1:nStates] <- Vect[k:(k+nStates-1)]
        k <- k + nStates
    }
    
    ResDistr <- SetVectAllParam(object=object$distribution, Vect=Vect[k:length(Vect)])
    HMM<-HMMSet(initProb=initProb, transMat=transMat, ResDistr)
    return(HMM)
}

GetNConstraint.HMMClass<-function(object)
{
    nStates <- object$distribution$nStates
    nConstr <- 1 + nStates + GetNConstraint(object=object$distribution)
    return(nConstr)
}

GradConstraint.HMMClass<-function(object)
{
    LParam <- GetNAllParam(object)
    nConstr <- GetNConstraint(object)
    Res <- matrix(0, nrow=nConstr, ncol=LParam$nParam)
    nStates <- object$distribution$nStates
    Grad1 <- rep(1, nStates)
#   probInit
    Res[1,1:nStates] <- Grad1
#   transMat
    for (i in 1:nStates)
        Res[i+1, (i*nStates+1):((i+1)*nStates)] <- Grad1
#   Distribution    
    Aux <- GradConstraint(object$distribution)
    if (! is.null(Aux))
    {    Res[(nStates+2):nConstr,(nStates*(nStates+1) + 1):LParam$nParam] <- Aux
    }
    return(Res)
}


GetNParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetNParam", object=object))

}

GetNAllParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetNAllParam", object=object))

}

GetVectParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetVectParam", object=object))

}

GetVectAllParam.distributionClass<-function(object)
{
    return(NextMethod(generic="GetVectAllParam", object=object))

}

SetVectParam.distributionClass<-function(object, Vect)
{
    return(NextMethod(generic="SetVectParam", object=object, Vect=Vect))

}

SetVectAllParam.distributionClass<-function(object, Vect)
{
    return(NextMethod(generic="SetVectAllParam", object=object, Vect=Vect))

}

GetNConstraint.distributionClass<-function(object)
{
    return(NextMethod(generic="GetNConstraint", object=object))
}

GradConstraint.distributionClass<-function(object)
{
    return(NextMethod(generic="GradConstraint", object=object))
}


GetNParam.univariateNormalClass<-function(object)
{
    return(list(nParam=2, paramConstr=c(0,eVar)))
}

GetNAllParam.univariateNormalClass<-function(object)
{
    return(list(nParam=2, paramConstr=c(0,eVar)))
}

GetVectParam.univariateNormalClass<-function(object)
{   nStates <- object$nStates
    nParam <- nStates * 2
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   Res[k] <- object$mean[i]
        k<-k+1
        Res[k] <- object$var[i]
        k<-k+1
    }
    return(Res)    
}

GetVectAllParam.univariateNormalClass<-function(object)
{   return(GetVectParam.univariateNormalClass(object))    
}

SetVectParam.univariateNormalClass<-function(object, Vect)
{   nStates <- object$nStates
    mean<-var<-rep(0, nStates)
    k <- 1
    for (i in 1:nStates)
    {   mean[i] <- Vect[k]
        k<-k+1
        var[i] <- Vect[k]
        k<-k+1
    }
    Res<-distributionSet(dis='NORMAL', mean=mean, var=var, verif=FALSE)
    return(Res)    
}

SetVectAllParam.univariateNormalClass<-function(object, Vect)
{   return(SetVectParam.univariateNormalClass(object, Vect)) 
}

GetNConstraint.univariateNormalClass<-function(object)
{
    return(0)
}

GradConstraint.univariateNormalClass<-function(object)
{
    return(NULL)
}

GetNParam.multivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nParam <- as.integer(dimObs + dimObs*(dimObs+1)/2)
    paramConstr <- c(rep(eNone,dimObs), rep(eVar, dimObs), rep(eCor, as.integer(dimObs*(dimObs-1)/2)))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetNAllParam.multivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nParam <- as.integer(dimObs + dimObs*(dimObs+1)/2)
    paramConstr <- c(rep(eNone,dimObs), rep(eVar, dimObs), rep(eCor, as.integer(dimObs*(dimObs-1)/2)))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetVectParam.multivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nStates <- object$nStates
    nParam <- as.integer(dimObs + dimObs*(dimObs+1)/2)*nStates
    Res <- rep(0, nParam)
# Moyennes, puis variances puis corrélations
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+dimObs-1)] <- object$mean[[i]]
        k <- k + dimObs
        Res[k:(k+dimObs-1)] <- diag(object$cov[[i]])
        k <- k + dimObs
        aa<-diag(1/sqrt(diag(object$cov[[i]])))
        matCor <- aa%*%object$cov[[i]]%*%aa
        for (n in 1:(dimObs-1))
        {   for (m in (n+1):dimObs)
            {   Res[k] <- matCor[m,n]
                k <- k + 1
            }
        }
    }
    return(Res)
}

GetVectAllParam.multivariateNormalClass<-function(object)
{
    return(GetVectParam.multivariateNormalClass(object))
}

SetVectParam.multivariateNormalClass<-function(object, Vect)
{
    dimObs <- object$dimObs
    nStates <- object$nStates
# Moyennes, puis variances puis corrélations
    k <- 1
    distr <- object 
    for (i in 1:nStates)
    {   distr$mean[[i]] <- Vect[k:(k+dimObs-1)]
        k <- k + dimObs
        aa <- sqrt(diag(Vect[k:(k+dimObs-1)]))
        k <- k + dimObs
        matCor <- diag(rep(1, dimObs))
        for (n in 1:(dimObs-1))
        {   for (m in (n+1):dimObs)
            {  matCor[m,n] <- matCor[n,m] <- Vect[k]
                k <- k + 1
            }
        }
        distr$cov[[i]] <- aa %*% matCor %*% aa
    }
    return(distr)
}

SetVectAllParam.multivariateNormalClass<-function(object, Vect)
{
    return(SetVectParam.multivariateNormalClass(object, Vect))
}

GetNConstraint.multivariateNormalClass<-function(object)
{
    return(0)
}

GradConstraint.multivariateNormalClass<-function(object)
{
    return(NULL)
}
   
GetNParam.mixtureUnivariateNormalClass<-function(object)
{
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nProba <- nMixt - 1
    nParam <- nMean + nVar + nProba
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eProba, nProba))
    return(list(nParam=nParam, paramConstr=paramConstr))
}


GetNAllParam.mixtureUnivariateNormalClass<-function(object)
{
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nParam <- nMean + nVar + nMixt
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eProba, nMixt))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetVectParam.mixtureUnivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nProba <- nMixt - 1
    nParam <- (nMean + nVar + nProba)*nStates
    Res <- rep(0, nParam)
#Les NMixt moyennes, puis les NMixt var puis les NMixt-1 proba
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+nMixt-1)] <- object$mean[[i]]
        k <- k + nMixt
        Res[k:(k+nMixt-1)] <- object$var[[i]]
        k <- k + nMixt
        Res[k:(k+nMixt-2)] <- object$prop[[i]][1:(nMixt-1)]
        k <- k + nMixt - 1
    }
    return(Res)
 }

GetVectAllParam.mixtureUnivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nParam <- (nMean + nVar + nMixt)*nStates
    Res <- rep(0, nParam)
#Les NMixt moyennes, puis les NMixt var puis les NMixt-1 proba
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+nMixt-1)] <- object$mean[[i]]
        k <- k + nMixt
        Res[k:(k+nMixt-1)] <- object$var[[i]]
        k <- k + nMixt
        Res[k:(k+nMixt-1)] <- object$prop[[i]]
        k <- k + nMixt
    }
    return(Res)
 }

SetVectParam.mixtureUnivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    nProba <- nMixt - 1
    distr <- object
#Les NMixt moyennes, puis les NMixt var puis les NMixt-1 proba
    k <- 1
    for (i in 1:nStates)
    {   distr$mean[[i]] <- Vect[k:(k+nMixt-1)]
        k <- k + nMixt
        distr$var[[i]] <- Vect[k:(k+nMixt-1)]
        k <- k + nMixt
        distr$prop[[i]][1:(nMixt-1)] <- Vect[k:(k+nMixt-2)]
        distr$prop[[i]][nMixt] <- 1.0 -sum(distr$prop[[i]][1:(nMixt-1)])
        k <- k + nMixt - 1
    }
    return(distr)
 }

SetVectAllParam.mixtureUnivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt
    distr <- object
#Les NMixt moyennes, puis les NMixt var puis les NMixt-1 proba
    k <- 1
    for (i in 1:nStates)
    {   distr$mean[[i]] <- Vect[k:(k+nMixt-1)]
        k <- k + nMixt
        distr$var[[i]] <- Vect[k:(k+nMixt-1)]
        k <- k + nMixt
        distr$prop[[i]] <- Vect[k:(k+nMixt-1)]
        k <- k + nMixt
    }
    return(distr)
 }

GetNConstraint.mixtureUnivariateNormalClass<-function(object)
{
    return(object$nStates)
}

GradConstraint.mixtureUnivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    nMixt <- object$nMixt
    LParam <- GetNAllParam(object)
    Res <- matrix(0, nrow=nStates, ncol=LParam$nParam*nStates)
    k <- 2*nMixt
    Grad1 <- rep(1, nMixt)
    for ( i in 1:nStates)
    {    Res[i,(k+1):(k+nMixt)] <- Grad1
        k <- k + 3*nMixt
    }
    return(Res)
}


 
GetNParam.mixtureMultivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nStates <- object$nStates
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    nProba <- (nMixt - 1)
    nParam <- nMean + nVar + nCor + nProba
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eCor, nCor), rep(eProba, nProba))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetNAllParam.mixtureMultivariateNormalClass<-function(object)
{
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nStates <- object$nStates
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    nParam <- nMean + nVar + nCor + nMixt
    paramConstr <- c(rep(eNone, nMean), rep(eVar, nVar), rep(eCor, nCor), rep(eProba, nMixt))
    return(list(nParam=nParam, paramConstr=paramConstr))
}

GetVectParam.mixtureMultivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    nProba <- nMixt - 1
    nParam <- (nMean + nVar + nCor + nProba)*nStates
    Res <- rep(0, nParam)
    k <- 1
# Les Nmixt x dimObs mean, puis les NMixt*dimObs Var, puis les cor, puis les proba      
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Res[k:(k+dimObs-1)] <- object$mean[[i]][[j]]
            k <- k + dimObs
        }
        for (j in 1:nMixt)
        {   aa <- diag(object$cov[[i]][[j]])
            Res[k:(k+dimObs-1)] <- aa
            k <- k + dimObs
        }
        for (j in 1:nMixt)
        {   aa <- object$cov[[i]][[j]]
            bb<-diag(1/sqrt(diag(aa)))
            matCor <- bb%*%aa%*%bb
            for (n in 1:(dimObs-1))
            {   for (m in (n+1):dimObs)
                {   Res[k] <- matCor[m,n]
                    k <- k + 1
                }
            }
        }
        Res[k:(k+nMixt-2)] <- object$proportion[[i]][1:(nMixt-1)]
        k <- k + nMixt -1
    }
    return(Res)
}

GetVectAllParam.mixtureMultivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
   nParam <- (nMean + nVar + nCor + nMixt)*nStates
    Res <- rep(0, nParam)
    k <- 1
# Les Nmixt x dimObs mean, puis les NMixt*dimObs Var, puis les cor, puis les proba      
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Res[k:(k+dimObs-1)] <- object$mean[[i]][[j]]
            k <- k + dimObs
        }
        for (j in 1:nMixt)
        {   aa <- diag(object$cov[[i]][[j]])
            Res[k:(k+dimObs-1)] <- aa
            k <- k + dimObs
        }
        for (j in 1:nMixt)
        {   aa <- object$cov[[i]][[j]]
            bb<-diag(1/sqrt(diag(aa)))
            matCor <- bb%*%aa%*%bb
            for (n in 1:(dimObs-1))
            {   for (m in (n+1):dimObs)
                {   Res[k] <- matCor[m,n]
                    k <- k + 1
                }
            }
        }
        Res[k:(k+nMixt-1)] <- object$proportion[[i]]
        k <- k + nMixt
    }
    return(Res)
}

SetVectParam.mixtureMultivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    k <- 1
    distr <- object
# Les Nmixt x dimObs mean, puis les NMixt*dimObs Var, puis les cor, puis les proba      
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   distr$mean[[i]][[j]] <- Vect[k:(k+dimObs-1)]
            k <- k + dimObs
        }
        for (j in 1:nMixt)
        {   distr$cov[[i]][[j]] <- diag(Vect[k:(k+dimObs-1)])
            k <- k + dimObs
        }
        
        for (j in 1:nMixt)
        {   matCor <- diag(rep(1, dimObs))
            for (n in 1:(dimObs-1))
            {   for (m in (n+1):dimObs)
                {   matCor[m,n] <- matCor[n,m] <- Vect[k]
                    k <- k + 1
                }
            }
            aa <- sqrt(distr$cov[[i]][[j]])
            distr$cov[[i]][[j]] <- aa %*% matCor %*% aa
        }
        distr$proportion[[i]][1:(nMixt-1)] <- Vect[k:(k+nMixt-2)]
        distr$proportion[[i]][nMixt] <- 1.0 - sum(distr$proportion[[i]][1:(nMixt-1)])
        k <- k + nMixt -1
    }
    return(distr)
}

SetVectAllParam.mixtureMultivariateNormalClass<-function(object, Vect)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    k <- 1
    distr <- object
# Les Nmixt x dimObs mean, puis les NMixt*dimObs Var, puis les cor, puis les proba      
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   distr$mean[[i]][[j]] <- Vect[k:(k+dimObs-1)]
            k <- k + dimObs
        }
        for (j in 1:nMixt)
        {   distr$cov[[i]][[j]] <- diag(Vect[k:(k+dimObs-1)])
            k <- k + dimObs
        }
        
        for (j in 1:nMixt)
        {   matCor <- diag(rep(1, dimObs))
            for (n in 1:(dimObs-1))
            {   for (m in (n+1):dimObs)
                {   matCor[m,n] <- matCor[n,m] <- Vect[k]
                    k <- k + 1
                }
            }
            aa <- sqrt(distr$cov[[i]][[j]])
            distr$cov[[i]][[j]] <- aa %*% matCor %*% aa
        }
        distr$proportion[[i]] <- Vect[k:(k+nMixt-1)]
        k <- k + nMixt
    }
    return(distr)
}

GetNConstraint.mixtureMultivariateNormalClass<-function(object)
{
    return(object$nStates)
}

GradConstraint.mixtureMultivariateNormalClass<-function(object)
{
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    LParam <- GetNAllParam(object)
    nMean <- nVar <- nMixt*dimObs
    nCor <- as.integer(nMixt*(dimObs*(dimObs-1)/2))
    k <- nOtherParam <- nMean + nVar + nCor
    Res <- matrix(0, nrow=nStates, ncol=LParam$nParam*nStates)
    Grad1 <- rep(1, nMixt)
    for ( i in 1:nStates)
    {   Res[i,(k+1):(k+nMixt)] <- Grad1
        k <- k + nMixt + nOtherParam
    }
    return(Res)
}

GetNParam.discreteClass <- function(object)
{
    nParam <- object$nLevels - 1
    return(list(nParam=nParam, paramConstr=rep(eProba, nParam)))
}

GetNAllParam.discreteClass <- function(object)
{
    nParam <- object$nLevels
    return(list(nParam=nParam, paramConstr=rep(eProba, nParam)))
}

GetVectParam.discreteClass <- function(object)
{   nStates <- object$nStates
    nLevels <- object$nLevels
    nParam <- (nLevels - 1)*nStates
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+nLevels-2)] <- as.vector(object$proba[[i]][1:(nLevels-1)]) 
        k <- k+nLevels-1
    }
    return(Res)
}

GetVectAllParam.discreteClass <- function(object)
{   nStates <- object$nStates
    nLevels <- object$nLevels
    nParam <- nLevels*nStates
    Res <- rep(0, nParam)
    k <- 1
    for (i in 1:nStates)
    {   Res[k:(k+nLevels-1)] <- as.vector(object$proba[[i]]) 
        k <- k+nLevels
    }
    return(Res)
}

SetVectParam.discreteClass <- function(object, Vect)
{
    nStates <- object$nStates
    nLevels <- object$nLevels
    distr <- object
    k <- 1
    for (i in 1:nStates)
    {   distr$proba[[i]][1:(nLevels-1)] <-  Vect[k:(k+nLevels-2)] 
        distr$proba[[i]][nLevels] <- 1.0 -sum(distr$proba[[i]][1:(nLevels-1)])
        k <- k+nLevels-1
    }
    return(distr)
}

SetVectAllParam.discreteClass <- function(object, Vect)
{
    nStates <- object$nStates
    nLevels <- object$nLevels
    distr <- object
    k <- 1
    for (i in 1:nStates)
    {   distr$proba[[i]] <-  Vect[k:(k+nLevels-1)] 
        k <- k+nLevels
    }
    return(distr)
}

GetNConstraint.discreteClass<-function(object)
{
    return(object$nStates)
}

logForwardBackward<-function(HMM, obs)
{   if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM
        
    if (is.data.frame(obs))
    {   dimObs <- dim(obs)[2]
        if (dimObs == 1)
            obs <- obs[,1]
        else
            obs <- as.matrix(obs[,1:dimObs])
    }

    HMM <- setStorageMode(HMM)
    maListe <- TransformeListe(HMM$distribution, obs)
    Res1 <- .Call("RLogforwardbackward", HMM, maListe$Zt)
    names(Res1) <- c("LLH")
    if (!is.list(obs))
    {   Res1$LLH <- Res1$LLH[[1]]
    }
     return(Res1)
}

NumLLH <- function(Teta, HMM, obs, nSample)
{
    HMM1<-SetVectAllParam(HMM, Teta)
    LLH <- sumList(logForwardBackward( HMM1, obs)$LLH, nSample)
    return(-LLH)
}

NumLLH1 <- function(Teta, HMM, obs, nSample)
{
    HMM1<-SetVectParam(HMM, Teta)
    LLH <- sumList(logForwardBackward( HMM1, obs)$LLH, nSample)
    return(-LLH)
}

optimHessian <- function(par, fn, ...) 
{
    fn1 <- function(par) fn(par, ...)
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
        length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100L, 
        abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
        beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
        factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
    hess <- try(.Internal(optimhess(par, fn1, NULL, con)), silent=TRUE)
    nParam <- length(par)
    if (class(hess) == "try-error")
    {   warning("non-finite finite-difference value while computing the Hessian matrix ...")
        hess <- matrix(0, nrow=nParam, ncol=nParam)
    }
    else
    {   hess <- 0.5 * (hess + t(hess))
    }

    return(hess)
}


nlmeHessian <- function(par, fn, ...) 
{
    fn1 <- function(par) fn(par, ...)
    hess <- try(fdHess(par, fn1)$Hessian, silent=TRUE)
    nParam <- length(par)
    if (class(hess) == "try-error")
    {   warning("non-finite finite-difference value while computing the Hessian matrix ...")
        hess <- matrix(0, nrow=nParam, ncol=nParam)
    }
    return(hess)
}

ComputeHessian <- function(HMM, obs, asymptMethod)
{
    if (is.list(obs))
    {   nSample <- length(obs)
    }
    else
    {   nSample <- 1
    }
    
    listeParam <- GetNAllParam(HMM)
    Teta0 <- GetVectAllParam(HMM)

    if (asymptMethod == 'nlme')
        Hess<- try(
                     nlmeHessian(par=Teta0, fn=NumLLH, HMM=HMM, obs=obs, nSample=nSample)
                )
    else
        Hess<- try(
                    optimHessian(par=Teta0, fn=NumLLH, HMM=HMM, obs=obs, nSample=nSample)
                )

    if (class(Hess) == "try-error")
       return(matrix(NaN, nrow=listeParam$nParam, ncol=listeParam$nParam))
    else
        return(Hess)
}

ComputeHessian1 <- function(HMM, obs, asymptMethod)
{
    if (is.list(obs))
    {   nSample <- length(obs)
    }
    else
    {   nSample <- 1
    }
    
    listeParam <- GetNParam(HMM)
    Teta0 <- GetVectParam(HMM)

    if (asymptMethod == 'nlme')
        Hess<- try(
                     nlmeHessian(par=Teta0, fn=NumLLH1, HMM=HMM, obs=obs, nSample=nSample)
                )
    else
        Hess<- try(
                    optimHessian(par=Teta0, fn=NumLLH1, HMM=HMM, obs=obs, nSample=nSample)
                )

    if (class(Hess) == "try-error")
       return(matrix(NaN, nrow=listeParam$nParam, ncol=listeParam$nParam))
    else
        return(Hess)
}


asymptoticCovMat <- function(HMM, obs, asymptMethod=c('nlme', 'optim'))
{
    if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM
    Hh <- ComputeHessian(HMM, obs, asymptMethod[1])
    nParam <- dim(Hh)[1]
    K <- GradConstraint(HMM)
    Dh <- Hh + t(K) %*% K
    Dm1 <- try(solve(Dh), silent=T)
    if (class(Dm1)=="try-error")
    {   asymptMatCov <- matrix(NaN, ncol=nParam, nrow=nParam)
        colnames(asymptMatCov) <- rownames(asymptMatCov) <- NomsParamHMM(HMM)
        warning("Hessian matrix is not inversible")
        return(asymptMatCov)
    }
    Aux <- K %*% Dm1 %*% t(K)
    Auxm1 <- try(solve(Aux), silent=T)
    if (class(Auxm1)=="try-error")
    {   asymptMatCov <- matrix(NaN, ncol=nParam, nrow=nParam)
        colnames(asymptMatCov) <- rownames(asymptMatCov) <- NomsParamHMM(HMM)
        warning("Hessian matrix is not inversible")
        return(asymptMatCov)
    }
    asymptMatCov <- Dm1 - Dm1 %*% t(K) %*% Auxm1 %*% K %*% Dm1
    if (is.list(obs))
    {   nTot <- length(obs[[1]])
        nSample <- length(obs)
        if (nSample > 1)
        {   for (j in 2:nSample)
                nTot <- nTot + length(obs[[j]])
        }
    }
    else
        nTot <- length(obs)

    colnames(asymptMatCov) <- rownames(asymptMatCov) <- NomsParamHMM(HMM)
    return(asymptMatCov)
}

asymptoticSimCovMat <- function(HMM, obs, nSimul, verbose=FALSE)
{
    if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM
    nParam <- GetNAllParam(HMM)$nParam
    Teta0 <- GetVectAllParam(HMM)
    matCov <- matrix(0, nParam, nParam)
    if (is.list(obs))
    {   nList <- length(obs)
        nObs <- rep(0, nList)
        for (j in 1:nList)
        {   nObs[j] <- length(obs[[j]])
        }   
    }
    else
    {   nList <- 1
        nObs <- c(length(obs))
    }
    for (i in 1:nSimul)
    {   simObs <- NULL
        for (j in 1:nList)
        {   simObs <- c(simObs, HMMSim(nObs[j], HMM=HMM))
        }
        if (HMM$distribution$dis=="NORMAL")
        {   Res <- HMMFit(simObs$obs, dis="NORMAL", nStates=HMM$distribution$nStates, asymptCov=FALSE, control=list(init="USER", initPoint=HMM))
            Teta <- GetVectAllParam(Res$HMM)
            u <- Teta - Teta0
            matCov <- (i*matCov + u%*%t(u))/(i+1)
            if (verbose)
            {   cat(sprintf("iteration %d\n", i))
            }
        }
    }
    colnames(matCov) <- rownames(matCov) <- NomsParamHMM(HMM)
    return(matCov)
}

setAsymptoticCovMat<-function(HMMFit, asymptCovMat)
{
    if ( ( class(HMMFit) != "HMMFitClass" ) )
        stop("class(HMMFit) must be 'HMMFitClass'\n")
    if (! is_numeric_matrix(asymptCovMat) )
        stop("asymptCovMat must be a matrix\n")
 #   if (! is_positive_definite(asymptCovMat) )
 #       stop("asymptCovMat must be a definite positive matrix\n")
    colnames(asymptCovMat) <- rownames(asymptCovMat) <- NomsParamHMM(HMM)
    HMM$asymptCov <- asymptCovMat
    return(HMM)
}
    
 
NomsParamHMM <- function(object) UseMethod("NomsParamHMM")

NomsParamHMM.default <- function(object)
{   return(NULL)
}

NomsParamHMM.distributionClass <- function(object)
{
    return(NextMethod(generic="NomsParamHMM", object=object))
}

NomsParamHMM.discreteClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    nLevels <- object$nLevels
    namesProba <- names(object$proba[[1]])
    if (is.null(namesProba))
    {   for (j in 1:nLevels)
        {   Aux <- sprintf("p[%d]", j)
            namesProba <- c(Res, Aux)
        }
    }   
    for (i in 1:nStates)
    {   for (j in 1:nLevels)
        {   Aux <- sprintf("State[%d]-%s", i, namesProba[j])
            Res <- c(Res, Aux)
        }
    }
    return(Res)
}

NomsParamHMM.univariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    for (i in 1:nStates)
    {   Aux <- sprintf("mean[%d]", i)
        Res <- c(Res, Aux)
        Aux <- sprintf("var[%d]", i)
        Res <- c(Res, Aux)
    }
    return(Res)
    
}

NomsParamHMM.mixtureUnivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    nMixt <- object$nMixt
    for (i in 1:nStates)
    {   for (j in 1:nMixt)
        {   Aux <- sprintf("state[%d]-mean[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:nMixt)
        {   Aux <- sprintf("state[%d]-var[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:nMixt)
        {   Aux <- sprintf("state[%d]-prop[%d]", i, j)
            Res <- c(Res, Aux)
        }
    
    }
    return(Res)   
}

NomsParamHMM.multivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    dimObs <- object$dimObs
    for (i in 1:nStates)
    {   for (j in 1:dimObs)
        {   Aux <- sprintf("state[%d]-mean[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:dimObs)
        {   Aux <- sprintf("state[%d]-var[%d]", i, j)
            Res <- c(Res, Aux)
        }
        
        for (j in 1:(dimObs-1))
        {   for (k in (j+1):dimObs)
            {   Aux <- sprintf("state[%d]-Cor[%d,%d]", i, j, k)
                Res <- c(Res, Aux)
            }
        }
    }
    return(Res)   
}

NomsParamHMM.mixtureMultivariateNormalClass <- function(object)
{
    Res <- NULL
    nStates <- object$nStates
    dimObs <- object$dimObs
    nMixt <- object$nMixt
    for (i in 1:nStates)
    {   for (m in 1:nMixt)
        {   for (j in 1:dimObs)
            {   Aux <- sprintf("state[%d]-Mixt[%d]-mean[%d]", i, m, j)
                Res <- c(Res, Aux)
            }
            for (j in 1:dimObs)
            {   Aux <- sprintf("state[%d]-Mixt[%d]-var[%d]", i, m, j)
                Res <- c(Res, Aux)
            }
            for (j in 1:(dimObs-1))
            {   for (k in (j+1):dimObs)
                {   Aux <- sprintf("state[%d]-Mixt[%d]-Cor[%d,%d]", i, m, j, k)
                    Res <- c(Res, Aux)
                }
            }
        }
        for (j in 1:nMixt)
        {   Aux <- sprintf("state[%d]-prop[%d]", i, j)
            Res <- c(Res, Aux)
        }
                
    }
    return(Res)   
}

NomsParamHMM.HMMClass <- function(object)
{
    nStates <- object$distribution$nStates
    Noms <- NULL
    for (i in 1:nStates)
    {   Aux <- sprintf("Pi[%d]", i)
        Noms <- c(Noms, Aux)
    }

    for (i in 1:nStates)
        for (j in 1:nStates)
        {   Aux <- sprintf("transMat[%d,%d]", i, j)
            Noms <- c(Noms, Aux)
        }
    Aux <- NomsParamHMM(object$distribution)
    Noms <- c(Noms, Aux)
    return(Noms)
}

NomsParamHMM.HMMFitClass <- function(object)
{
    return(NomsParamHMM(object$HMM))
}

summary.HMMFitClass <-function (object, ...) 
{
    ans = NULL
    ans$call = object$call
    ans$nIter <- object$nIter
    ans$relVariation <- object$relVariation
    y <- object$HMM$distribution
    if (y$dis == "NORMAL")
    {    if (y$dimObs==1)
            nomloi <- "univariate gaussian"
        else
            nomloi <- sprintf("%d-d gaussian", y$dimObs)
    }
    if (y$dis == "DISCRETE")
        nomloi <- "discrete"
    if (y$dis == "MIXTURE")
    {   if (y$dimObs == 1)
            nomloi <- sprintf("mixture of %d gaussian", y$nMixt)
        else
            nomloi <- sprintf("mixture of %d %d-d gaussian", y$nMixt, y$dimObs)
    }
    Model <- sprintf("%d states HMM with %s distribution", y$nStates, nomloi)     
    ans$model <- Model
    ans$LLH = object$LLH
    ans$BIC <- object$BIC
    if (is.null(object$asymptCov))
    {   cat(sprintf("Computing the asymptotic covariance matrix of estimates\n"))
        object$asymptCov <- asymptoticCovMat(object, object$obs)     
    }

    nAllParam <- GetNAllParam(object)
    Value <- GetVectAllParam(object) 
    asymptVar <- diag(object$asymptCov)
    se.coef <- sqrt(asymptVar)
    tval = Value/se.coef
    prob = 2 * (1 - pnorm(abs(tval)))
    
    Noms <- NomsParamHMM(object$HMM)
    ans$coef = cbind(Value, se.coef, tval, prob)
    dimnames(ans$coef) = list(Noms, c(" Estimate", 
        " Std. Error", " t value", "Pr(>|t|)"))
    class(ans) <- 'summary.HMMFitClass'
    return(ans)
}

print.summary.HMMFitClass <- function(x, ...)
{  # Description:
    #   Print summary method for an x of class "HMMFitClass".
    
    # FUNCTION:

    # Call and Model:
    cat("\nCall:", sep="\n")
    cat("----", sep="\n")
    cat(deparse(x$call), "\n", sep = "")

    cat("\nModel:", sep="\n")
    cat("------", sep="\n")
    cat(x$model, "\n", sep = "")
    cat("\nBaum-Welch algorithm status:", sep="\n")
    cat("----------------------------", sep="\n")
    cat(sprintf("Number of iterations : %d\n", x$nIter), sep="")
    cat(sprintf("Last relative variation of LLH function: %f\n", x$relVariation), sep="")

    # Coef Matrix:
    cat("\nCoefficient(s):\n")
    cat("---------------\n")
    signif.stars = getOption("show.signif.stars")
    digits = max(4, getOption("digits") - 4)
    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
    
     cat("\nLog Likelihood: ", 
        format(round(x$LLH, 2)), "\n")

    cat("BIC Criterion: ", 
        format(round(x$BIC, 2)), "\n")
    
    # Return Value:
#    cat("\n")
    invisible()
}
