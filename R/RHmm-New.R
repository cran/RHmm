tolMin <- .Machine$double.eps*100
AsymptProbaVector <-function(probVector, obs, indDeb, HMM, Delta=1e-4)
{
    if (is.list(obs))
    {   nSample <- length(obs)
    }
    else
    {   nSample <- 1
    }
    
    nProba <- length(probVector)
    Res0 <- setAsymptVarProbaVector(probVector, 2*Delta)
    Res0$hh <- Res0$hh/2
    
    LLHi <- rep(0, (nProba-1))
    LLHij <- matrix(0, nrow=nProba-1, ncol=nProba-1)
    Hess <- LLHij
    nParam <- GetNParam(HMM)
    Teta0 <- GetAllParam(HMM, nParam)
    LLH0 <- sumList(forwardbackward(HMM, obs)$LLH, nSample)
     
    for (i in 1:(nProba-1))
    {   Tetai <- Teta0
        Tetai[indDeb+i-1] <- Tetai[indDeb+i-1] + Res0$hh[i]*Res0$sens[i]
        Tetai[indDeb+nProba-1] <- Tetai[indDeb+nProba-1] - Res0$hh[i]*Res0$sens[i]
        HMMi <- SetAllParam(HMM, Tetai)
        LLH1 <- forwardbackward(HMMi, obs)$LLH
        LLHi[i] <- sumList(LLH1, nSample)
    }
    
    for (i in 1:(nProba-1))
    {   for (j in i:(nProba-1))
        {   Tetaij <- Teta0
            Tetaij[indDeb+i-1] <- Tetaij[indDeb+i-1] + Res0$hh[i]*Res0$sens[i]
            Tetaij[indDeb+nProba-1] <- Tetaij[indDeb+nProba-1] - Res0$hh[i]*Res0$sens[i]
            Tetaij[indDeb+j-1] <- Tetaij[indDeb+j-1] + Res0$hh[j]*Res0$sens[j]
            Tetaij[indDeb+nProba-1] <- Tetaij[indDeb+nProba-1] - Res0$hh[j]*Res0$sens[j]
            HMMij <- SetAllParam(HMM, Tetaij)
            LLH2 <- forwardbackward(HMMij, obs)$LLH
            LLHij[i,j] <- sumList(LLH2, nSample)
            LLHij[j,i] <- LLHij[i,j]
            Hess[i,j] <- (LLHij[i,j] - LLHi[i] - LLHi[j] + LLH0)/(Res0$hh[i]*Res0$sens[i]*Res0$hh[j]*Res0$sens[j])
            Hess[j,i] <- Hess[i,j]
        }
    }
    if (abs(min(eigen(-Hess)$values)) < tolMin)
    {    varProbVector <- rep(NaN, nProba)
    }
    else
    {   matCov <- solve(-Hess)
        varProbVector <- diag(matCov)
        Un <- rep(1, (nProba-1))
        varProbVector[nProba] <- t(Un) %*% matCov %*% Un
    }
    return(varProbVector)
}



AsymptVar <- function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL) UseMethod("AsymptVar")

AsymptVar.default <-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
{   cat("object must be a HMMClass or HMMFitClass or distribution object\n")
    return(invisible(0))
}

AsymptVar.HMMFitClass<-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
{
     AsymptVar(object=object$HMM, obs, Delta, initProb, transMat)
}

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

AsymptVar.distributionClass<-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
{
    NextMethod(generic="AsymptVar", object=object)
}


AsymptVar.univariateNormalClass<-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
#Pas de contraintes sur les paramètres
{   
    HMM0 <- HMMSet(initProb=initProb, transMat=transMat, distribution=object) 
    if (is.list(obs))
        nSample <- length(obs)
    else
        nSample <- 1
    
    LLH0 <- sumList(forwardbackward(HMM0, obs)$LLH, nSample)
    nStates <- object$nStates
    nParam <- GetNParam(object)*nStates
    Teta0 <- GetDistrParam(object, nParam)
    Deriv <- rep(0, nParam)
    for (i in 1:nParam)
    {   hh <- abs(Delta * Teta0[i])
        while (hh < tolMin)
        {     hh <- 2.0 * hh
        }
        Teta1 <- Teta0
        Teta1[i] <- Teta1[i] + hh
        Teta2 <- Teta1
        Teta2[i] <- Teta2[i] + hh
        distr1 <- SetDistrParam(object, Teta1)
        distr2 <- SetDistrParam(object, Teta2)
        HMM1 <- HMMSet(initProb=initProb, transMat=transMat, distribution=distr1)
        HMM2 <- HMMSet(initProb=initProb, transMat=transMat, distribution=distr2)
        LLH2 <- sumList(forwardbackward(HMM2, obs)$LLH, nSample)
        LLH1 <- sumList(forwardbackward(HMM1, obs)$LLH, nSample)
        
        Aux <- (LLH2 - 2*LLH1 + LLH0)/(hh^2)
        Deriv[i] <- -1.0/Aux
    }
    distr <- SetDistrParam(object, Deriv)
    return(distr)
}

AsymptVar.multivariateNormalClass<-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
#Pas de contraintes sur les paramètres
{   
    HMM0 <- HMMSet(initProb=initProb, transMat=transMat, distribution=object) 
    if (is.list(obs))
        nSample <- length(obs)
    else
        nSample <- 1
    
    LLH0 <- sumList(forwardbackward(HMM0, obs)$LLH, nSample)
    nStates <- object$nStates
    nParam <- GetNParam(object)*nStates
    Teta0 <- GetDistrParam(object, nParam)
    Deriv <- rep(0, nParam)
    for (i in 1:nParam)
    {   hh <- abs(Delta * Teta0[i])
        if (hh < tolMin)
             hh <-  tolMin
        
        Teta1 <- Teta0
        Teta1[i] <- Teta1[i] + hh
        Teta2 <- Teta1
        Teta2[i] <- Teta2[i] + hh
        distr1 <- SetDistrParam(object, Teta1)
        distr2 <- SetDistrParam(object, Teta2)
        HMM1 <- HMMSet(initProb=initProb, transMat=transMat, distribution=distr1)
        HMM2 <- HMMSet(initProb=initProb, transMat=transMat, distribution=distr2)
        LLH2 <- sumList(forwardbackward(HMM2, obs)$LLH, nSample)
        LLH1 <- sumList(forwardbackward(HMM1, obs)$LLH, nSample)
        
        Aux <- (LLH2 - 2*LLH1 + LLH0)/(hh^2)
        Deriv[i] <- -1.0/Aux
    }
    distr <- SetDistrParam(object, Deriv)
    return(distr)
}

AsymptVar.discreteClass<-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
#Avec des contraintes sur les paramètres
{   
    HMM0 <- HMMSet(initProb=initProb, transMat=transMat, distribution=object) 
    nStates <- object$nStates
    k <- 1 + nStates*(nStates+1)
    varProba <- object$proba
    for (i in 1:nStates)
    {   varProba[[i]] <- AsymptProbaVector(object$proba[[i]], obs, k, HMM0, Delta)
        k <- k +object$nLevels
    }
    distr <- SetDistrParam(object, varProba)
    return(object)
}

AsymptVar.mixtureUnivariateNormalClass<-function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
#avec des contraintes sur les paramètres
{   
    HMM0 <- HMMSet(initProb=initProb, transMat=transMat, distribution=object) 
    if (is.list(obs))
        nSample <- length(obs)
    else
        nSample <- 1
    
    LLH0 <- sumList(forwardbackward(HMM0, obs)$LLH, nSample)
    nStates <- object$nStates
    nMixt <- object$nMixt
    nParam <- GetNParam(object)*nStates
    Teta0 <- GetDistrParam(object, nParam)
    Deriv <- rep(0, nParam)
    k <- 0
    for (i in 1:nStates) #nMixt mean vector + nMixt var vector + nMixt prop vector
    {   for (j in 1:(2*nMixt))
        {   k <- k + 1
            hh <- abs(Delta * Teta0[k])
            if (hh < tolMin)
                hh <- tolMin
            
            Teta1 <- Teta0
            Teta1[k] <- Teta1[k] + hh
            Teta2 <- Teta1
            Teta2[k] <- Teta2[k] + hh
            distr1 <- SetDistrParam(object, Teta1)
            distr2 <- SetDistrParam(object, Teta2)
            HMM1 <- HMMSet(initProb=initProb, transMat=transMat, distribution=distr1)
            HMM2 <- HMMSet(initProb=initProb, transMat=transMat, distribution=distr2)
            LLH2 <- sumList(forwardbackward(HMM2, obs)$LLH, nSample)
            LLH1 <- sumList(forwardbackward(HMM1, obs)$LLH, nSample)
        
            Aux <- (LLH2 - 2*LLH1 + LLH0)/(hh^2)
            Deriv[k] <- -1.0/Aux
        }
        propVector <- Teta0[(k+1):(k+nMixt)]
        Deriv[(k+1):(k+nMixt)] <- AsymptProbaVector(propVector, obs, k+1, HMM0, Delta)
        k <- k + nMixt
    }
    distr <- SetDistrParam(object, Deriv)
    return(distr)
}

AsymptVar.HMMClass <- function(object, obs, Delta=1e-4, initProb=NULL, transMat=NULL)
{
    HMM0 <- object
    LLH0 <- forwardbackward(HMM0, obs)$LLH
    if (is.list(obs))
    {
        nSample <- length(obs)
    }
    else
    {    nSample <- 1
    }
    LLH0 <- sumList(LLH0, nSample)
    nStates <- length(object$initProb)
     
# initProb
 
    varInitProb<-AsymptProbaVector(HMM0$initProb, obs, 1, HMM0, Delta)
             
# transMat
    varTransMat <- matrix(0, nStates, nStates)
    k <- nStates+1
    for (m in 1:nStates)
    {   varTransMat[m, ] <- AsymptProbaVector(HMM0$transMat[m,], obs, k, HMM0, Delta)
        k <- k +nStates   
    }
    
    distrVar <- AsymptVar(object$distribution, obs, Delta, HMM0$initProb, HMM0$transMat) 

    return(HMMSet(initProb=varInitProb, transMat=varTransMat, distribution=distrVar))
    
}



GetAllParam <- function(object, nParam) UseMethod("GetAllParam")
SetAllParam <- function(object, nParam) UseMethod("SetAllParam")

GetDistrParam <- function(object, nDistrParam) UseMethod("GetDistrParam")
SetDistrParam <- function(object, Vector) UseMethod("SetDistrParam")

GetDistrParam.default <- function(object, nDistrParam)
{
    return(NULL)
}

GetDistrParam.discreteClass <- function(object, nDistrParam)
{
    k <- 0
    Res <- rep(0, nDistrParam)
    for (i in 1:object$nStates)
    {   Res[(k+1):(k+object$nLevels)] <- object$proba[[i]]
        k <- k + object$nLevels
    }   
    return(Res)
}

SetDistrParam.discreteClass <- function(object, Vector)
{
    k <- 0
    for (i in 1:object$nStates)
    {   Vector[(k+1):(k+object$nLevels)] -> object$proba[[i]]
        k <- k + object$nLevels
    }   
    return(object)
}

GetDistrParam.univariateNormalClass<-function(object, nDistrParam)
{
    k <- 1
    Res <- rep(0, nDistrParam)
    for (i in 1:object$nStates)
    {   Res[k] <- object$mean[i]
        Res[k+1] <- object$var[i]
        k <- k+2
    }
    return(Res)
 }

SetDistrParam.univariateNormalClass<-function(object, Vector)
{
    k <- 1
    for (i in 1:object$nStates)
    {   Vector[k] -> object$mean[i]
        Vector[k+1] -> object$var[i]
        k <- k+2
    }
    return(object)
 }
 
GetDistrParam.multivariateNormalClass<-function(object, nDistrParam)
{
    k <- 0
    Res <- rep(0, nDistrParam)
    for (i in 1:object$nStates)
    {   Res[(k+1):(k+object$dimObs)] <- object$mean[[i]]
        k <- k + object$dimObs
        for (j in 1:object$dimObs)
        {   Res[(k+1):(k+object$dimObs-j+1)] <- object$cov[[i]][j,j:object$dimObs]
            k <- k + object$dimObs -j + 1
        }
     }
    return(Res)
 }
 
SetDistrParam.multivariateNormalClass<-function(object, Vector)
{
    k <- 0
    for (i in 1:object$nStates)
    {   Vector[(k+1):(k+object$dimObs)] -> object$mean[[i]]
        k <- k + object$dimObs
        for (j in 1:object$dimObs)
        {   Vector[(k+1):(k+object$dimObs-j+1)] -> object$cov[[i]][j,j:object$dimObs]
            k <- k + object$dimObs -j + 1        
        }
        for (m in 1:(object$dimObs-1))
        {    for (n in (m+1):object$dimObs)
                object$cov[[i]][n, m] <- object$cov[[i]][m, n]
        }
    }
    return(object)
 }

GetDistrParam.mixtureUnivariateNormalClass<-function(object, nDistrParam)
{
    k <- 0
    Res <- rep(0, nDistrParam)
    for (i in 1:object$nStates)
    {   Res[(k+1):(k+object$nMixt)] <- object$mean[[i]]
        k <- k + object$nMixt
        Res[(k+1):(k+object$nMixt)] <- object$var[[i]]
        k <- k + object$nMixt
        Res[(k+1):(k+object$nMixt)] <- object$proportion[[i]]
        k <- k + object$nMixt
    }
    return(Res)
 }
 
SetDistrParam.mixtureUnivariateNormalClass<-function(object, Vector)
{
    k <- 0
    nMixt <- object$nMixt
    for (i in 1:object$nStates)
    {   Vector[(k+1):(k+nMixt)] -> object$mean[[i]]
        k <- k + nMixt
        Vector[(k+1):(k+nMixt)] -> object$var[[i]]
        k <- k + nMixt
        Vector[(k+1):(k+nMixt)] -> object$proportion[[i]]
        k <- k + nMixt
    }
    return(object)
 }

GetDistrParam.mixtureMultivariateNormalClass<-function(object, nDistrParam)
{
    k <- 0
    Res <- rep(0, nDistrParam)
    for (i in 1:object$nStates)
    {   for (j in 1:object$nMixt)
        {   Res[(k+1)] <- object$proportion[[i]][[j]]
            Res[(k+1):(k+object$dimObs)] <- object$mean[[i]][[j]]
            k <- k + object$dimObs
            for (l in 1:object$dimObs)
            {    Res[(k+1):(k+object$dimObs-l+1)] <- object$cov[[i]][[j]][l,l:object$dimObs]
                k <- k + object$dimObs-l+1
            }
            k <- k + 1
        }
    }
    return(Res)
 }

SetDistrParam.mixtureMultivariateNormalClass<-function(object, Vector)
{
    k <- 0
    for (i in 1:object$nStates)
    {   for (j in 1:object$nMixt)
        {   Vector[(k+1)] -> object$proportion[[i]][[j]]
            Vector[(k+1):(k+object$dimObs)] -> object$mean[[i]][[j]]
            k <- k + object$dimObs
            for (l in 1:object$dimObs)
            {   Vector[(k+1):(k+object$dimObs-l+1)] -> object$cov[[i]][[j]][l,l:object$dimObs]
                k <- k + object$dimObs-l+1
            }
            for (m in 1:(object$dimObs-1))
            {    for (n in (m+1):object$dimObs)
                    object$cov[[i]][[j]][n, m] <- object$cov[[i]][[j]][m, n]
            }

             k <- k + 1
        }
    }
    return(object)
 }

GetNParam<-function(object) UseMethod("GetNParam")

GetNParam.default<-function(object)
{
    return(0)
}

GetNParam.HMMFitClass<-function(object)
{
    return(GetNParam(object$HMM))
}

GetNParam.HMMClass<-function(object)
{
    nStates <- object$distribution$nStates
    Res <- nStates*(1 + nStates + GetNParam(object$distribution))
    return(Res)
}

GetNParam.univariateNormalClass<-function(object)
{
    return(2)
}

GetNParam.multivariateNormalClass<-function(object)
{
    return(as.integer(object$dimObs + object$dimObs*(object$dimObs + 1)/2))
}

GetNParam.mixtureUnivariateNormalClass<-function(object)
{
    return(3*object$nMixt)
}

GetNParam.mixtureMultivariateNormalClass<-function(object)
{
    return(as.integer(object$nMixt*(1 + object$dimObs + object$dimObs*(object$DimObs+1)/2)))
}

setAsymptVarProbaVector<-function(Vector, Delta)
{
    n <- length(Vector)
    sens <- rep(0, n-1)
    hh <- sens
    for (i in 1:(n-1))
    {   hh0 <- abs(Delta*Vector[i])
        if (hh0 < tolMin)
        {    hh0 <- Delta
        }
        if (Vector[i] + hh0 > 1)
        {   sens[i] <- -1
        }
        else
        {    sens[i] <- 1
        }
        hh[i] <- hh0
    }
    return(list(hh = hh, sens=sens))
}

GetAllParam <- function(object, nParam) UseMethod("GetAllParam")
SetAllParam <- function(object, Vector) UseMethod("SetAllParam")

GetAllParam.HMMFitClass <- function(object, nParam)
{
    return(GetAllParam(object$HMM, nParam))
}

SetAllParam.HMMFitClass <- function(object, Vector)
{
    return(SetAllParam(object$HMM, Vector))
}


GetAllParam.HMMClass <- function(object, nParam)
{
    nStates <- object$distribution$nStates
    Teta <- rep(0, nStates*(nStates + 1))
    
    Teta[1:nStates] <- object$initProb
    k<-nStates
    for (i in 1:nStates)
    {   Teta[(k+1):(k+nStates)] <- object$transMat[i,]
        k <- k + nStates
    }
    Teta <- c(Teta, GetDistrParam(object$distribution,GetNParam(object$distribution)))
    return(Teta)
}

SetAllParam.HMMClass <- function(object, Vector)
{   nStates <- object$distribution$nStates
    
    Vector[1:nStates] -> object$initProb
    k<-nStates
    for (i in 1:nStates)
    {   Vector[(k+1):(k+nStates)] -> object$transMat[i,]
        k <- k + nStates
    }
    object$distribution <- SetDistrParam(object$distribution,Vector[(k+1):length(Vector)])
    return(object)
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
    
    for (i in 1:nStates)
    {   namesProba <- names(object$proba[[i]])
        if (is.null(namesProba))
        {   for (j in 1:nLevels)
            {   Aux <- sprintf("p[%d]", j)
                Res <- c(Res, Aux)
            }
        }          
        for (j in 1:nLevels)
        {   Aux <- sprintf("State %d - %s", i, namesProba[j])
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
        {   Aux <- sprintf("state %d - mean[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:nMixt)
        {   Aux <- sprintf("state %d - var[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:nMixt)
        {   Aux <- sprintf("state %d - prop[%d]", i, j)
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
        {   Aux <- sprintf("state %d - mean[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:dimObs)
        {   for (k in j:dimObs)
            {   Aux <- sprintf("state %d - Cov[%d,%d]", i, j, k)
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
    for (i in 1:nStates)
    {   for (j in 1:dimObs)
        {   Aux <- sprintf("state %d - mean[%d]", i, j)
            Res <- c(Res, Aux)
        }
        for (j in 1:dimObs)
        {   for (k in j:dimObs)
            {   Aux <- sprintf("state %d - Cov[%d,%d]", i, j, k)
                Res <- c(Res, Aux)
            }
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
        {   Aux <- sprintf("transMat[%d, %d]", i, j)
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
    
    nParam <- GetNParam(object)
    Value <- GetAllParam(object)
    if (!is.null(object$asymptVar))
    {   AsymptVariance <- GetAllParam(object$asymptVar)
        se.coef <- sqrt(AsymptVariance)
        tval = Value/se.coef
        prob = 2 * (1 - pnorm(abs(tval)))
    }
    else
    {   AsymptVariance <- rep(NaN, nParam)
        se.coef <- AsymptVariance
        tval <- AsymptVariance
        prob <- AsymptVariance
    }
    
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
