########## RHmm-Internals
setStorageMode <- function(object) UseMethod("setStorageMode")

setStorageMode.paramHMM <- function(object)
{   x <- object
    storage.mode(x$nStates) <- "integer"
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$nMixt) <- "integer"
    storage.mode(x$nLevels) <- "integer"
    storage.mode(x$noHMM) <-"integer"
    storage.mode(x) <- "list"
    class(x) <- "paramHMM"
    return(x)
}

setStorageMode.paramAlgoBW <- function(object)
{   x <- object
    storage.mode(x$iter) <- "integer"
    storage.mode(x$verbose) <- "integer"
    storage.mode(x$nInit) <- "integer"
    storage.mode(x$nIterInit) <- "integer"
    if (!is.null(x$initPoint))
       x$initPoint <- setStorageMode(x$initPoint)
    storage.mode(x) <- "list"
    class(x) <- "paramAlgoBW"
    return(x)
}    

setStorageMode.HMMClass <- function(object)
{   x <- object
    storage.mode(x$initProb) <- "double"
    storage.mode(x$transMat) <- "double"
    x$distribution <- setStorageMode(object$distribution)
    class(x$distribution) <- class(object$distribution)
    storage.mode(x) <- "list"
    class(x) <- "HMMClass"
    return(x)
}

setStorageMode.distributionClass <- function(object)
{   x<-NextMethod("setStorageMode", object)
    storage.mode(x$nStates) <- "integer"
    return(x)
}

setStorageMode.univariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$mean) <- "double"
    storage.mode(x$var) <- "double"
    storage.mode(x) <- "list"
    class(x) <- c("distributionClass", "univariateNormalClass")
    return(x)
}

setStorageMode.multivariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$mean) <- "list"
    storage.mode(x$cov) <- "list"
    storage.mode(x) <- "list"
    class(x) <- c("distributionClass", "multivariateNormalClass")
    return(x)
}

setStorageMode.mixtureUnivariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$nMixt) <- "integer"
    storage.mode(x$mean) <- "list"
    storage.mode(x$var) <- "list"
    storage.mode(x$proportion) <- "list"
    storage.mode(x) <- "list"
    return(x)        
}

setStorageMode.mixtureMultivariateNormalClass <- function(object)
{   x <- object
    storage.mode(x$dimObs) <- "integer"
    storage.mode(x$nMixt) <- "integer"
    storage.mode(x$mean) <- "list"
    storage.mode(x$cov) <- "list"
    storage.mode(x$proportion) <- "list"
    storage.mode(x) <- "list"
    return(x)        
}

setStorageMode.discreteClass <- function(object)
{  x <- object
   storage.mode(x$nLevels) <- "integer"
   storage.mode(x) <- "list"
   return(x)        
}    

make_labels <- function(factorlist)
{   l <- length(factorlist)
    Aux <- as.vector(factorlist[[1]])
    Res <- list(1)
    for (i in 2:l)
    {   x <- as.vector(factorlist[[i]])
        Aux <- c(Aux, x)
        Res <- c(Res, list(1))
    }
    
    Aux <- factor(Aux)
    labels <- levels(Aux)
    nLevels <- length(labels)
    for (i in 1:l)
        Res[[i]] <- factor(factorlist[[i]], levels=labels)
    
    return(list(labels = labels, nLevels=nLevels, obs=Res))
}

TransformeList <- function(dis, obs)
{   labels <- NULL
    nLevels <- 0
    if (is.list(obs))
    {   if (dis == "DISCRETE")
        {   Aux <- make_labels(obs)
            nLevels <- Aux$nLevels
            labels <- Aux$labels
            Zt <- list(as.double(Aux$obs[[1]])-1)
            for (i in 2:length(obs))
                Zt <- c(Zt, list(as.double(Aux$obs[[i]])-1))
        }
        else
            Zt <- obs
    }
    else
    {   if (dis == "DISCRETE")
        {   Aux <- factor(obs)
            labels <- levels(Aux)
            nLevels <- length(labels)
            Zt <- list(as.double(factor(obs))-1) 
        }
        else
            Zt <- list(obs)
    }
    return(list(Zt=Zt, nLevels = nLevels, labels = labels))
}
########## end of RHmm_internals
