tol <- .Machine$double.eps
min_eigen_value <- function(x)
{
    return(min(eigen(x)$values))
}

is_numeric_vector <- function(x)
{   if (!is.vector(x))
        return(FALSE)
    if (!is.numeric(x))
        return(FALSE)
    return(TRUE)
}

is_list_numeric_vector <- function(x)
{
    if (!is.list(x))
        return(FALSE)
    
    if (!all(sapply(x, is_numeric_vector)))
        return(FALSE)
    
    return(TRUE)
} 

is_list_list <- function(x)
{
    if (!is.list(x))
        return(FALSE)
    if (!all(sapply(x, is.list)))
        return(FALSE)
    return(TRUE)
}

is_list_list_numeric_vector <- function(x)
{
    if (!is.list(x))
        return(FALSE)
    if (!all(sapply(x, is_list_numeric_vector)))
        return(FALSE)
    return(TRUE)
}

is_numeric_matrix <- function(x)
{
   if (!is.matrix(x))
        return(FALSE)
    if (!is.numeric(x))
        return(FALSE)
    return(TRUE)    
}
 
is_list_numeric_matrix <- function(x)
{
   if (!is.list(x))
        return(FALSE)
    if (!all(sapply(x, is_numeric_matrix)))
        return(FALSE)
    return(TRUE)
}

is_var <- function(x)
{   return(all(x>=0))
}

is_list_var <- function(x)
{
    return(all(sapply(x, is_var)))
}

is_proba_vector <- function(x)
{
    if (!all(x >=0))
        return(FALSE)
    return(abs(sum(x)-1) <= 2*tol)
}

is_list_proba_vector <- function(x)
{
    return(all(sapply(x, is_proba_vector)))
}

s_list_list_proba_vector <- function(x)
{
    return(all(sapply(x, is_list_proba_vector)))
}

is_positive_definite <- function(x)
{
   return( min_eigen_value(x) > 0 )
}

is_list_positive_definite <- function(x)
{   if (is.null(x))
        return(TRUE)
    return(all(sapply(x,is_positive_definite)))
} 

is_list_list_positive_definite <- function(x)
{
    if (is.null(x))
        return(TRUE)
    if (!is.list(x))
        return(FALSE)
    if (!all(sapply(x, is_list_numeric_matrix)))
        return(FALSE)
    if(!all(sapply(x, is_list_positive_definite)))
        return(FALSE)
    return(TRUE)
}

univariateNormalSet <- function(mean, var, verif=TRUE)
{   if (verif)
    {   if (!is_numeric_vector(mean))
            return("'mean' must be a vector for univariate normal distributions.\n")
        if (!is_numeric_vector(var))
            return("'var' must be a vector for univariate normal distributions.\n")
        nStates <- length(mean)
        if (nStates != length(var))
            return("'mean' and 'var' parameters must have the same length, which is the number of hidden states.\n")
        if ( !is_var(var) )
            return("all variances must be positive\n")
    }
    else
        nStates <- length(mean)
    
    value <- list(dis="NORMAL", nStates=nStates, dimObs=as.integer(1), mean=mean, var=var)
    class(value) <- c("distributionClass", "univariateNormalClass")
    return(value)    
}

multivariateNormalSet <- function(mean, cov, verif=TRUE)
{   if (verif)
    {   if (!is_list_numeric_vector(mean))
            return("'mean' must be a list of vectors for multivariate normal distributions.\n")
        if (!is_list_numeric_matrix(cov))
            return("'cov' must be a list of matrices for multivariate normal distributions.\n")
        nStates <- length(mean)
        if (nStates != length(cov))
            return("'mean' and 'cov' parameters must have the same length, which is the number of hidden states.\n")
        dimObs <- length(mean[[1]])
        Aux <- as.integer(lapply(mean, length)) - dimObs
        if (sum(abs(Aux)) != 0)
            return("Each element of 'mean' must have the same length, which is the dimension of the observations.\n")
        Aux <- as.integer(lapply(cov, ncol)) - dimObs
        if (sum(abs(Aux)) != 0)
            return("The number of columns of each element of 'cov' must be equal to the dimension of the observations.\n")
        Aux <- as.integer(lapply(cov, nrow)) - dimObs
        if (sum(abs(Aux)) != 0)
            return("The number of row of each element of 'cov' must be equal to the dimension of the observations.\n")
        if (!is_list_positive_definite(cov))
            return("Each element of 'cov' must be a semi definite positive matrix\n")
     }
    else
    {   nStates <- length(mean)
        dimObs <- length(mean[[1]])
    }
    Res <- list(dis="NORMAL", nStates = nStates, dimObs=dimObs, mean=mean, cov=cov)
    class(Res) <- c("distributionClass", "multivariateNormalClass")
    return(Res)    
}

discreteSet <- function(proba, labels=NULL, verif = TRUE)
{   
    if (verif)
    {   if (!is_list_numeric_vector(proba))
            return("'proba' parameter for discrete distributions must be a list of vectors.\n")
        nStates <- length(proba)
       
        nLevels <- length(proba[[1]])
        Aux <- as.integer(lapply(proba, length)) - nLevels
        if (sum(abs(Aux)) != 0)
            return("Each element of 'proba' must have the same length, which is the number of the different discrete observations.\n") 
        if (!is_list_proba_vector(proba) )
            return("The sum of each element of 'proba' must be equal to 1.\n")
    }
    else
    {   nStates <- length(proba)
        nLevels <- length(proba[[1]])
    }
    
    if (!is.null(labels))
    {   Aux <- as.integer(lapply(labels, is.character))
        if ( (sum(Aux) != length(Aux)) && verif )
            return("'labels' must be a vector of characters.\n")
        if ( (length(labels) != nLevels) && verif )
            return("wrong number of labels\n")
    }
    else
        labels <- as.character(gl(nLevels, 1, label="p"))
    
    for (i in 1:nStates)
        names(proba[[i]]) <- labels
    
         
    Res <- list(dis="DISCRETE", nStates=nStates, nLevels=nLevels, proba=proba)
    class(Res) <- c("distributionClass", "discreteClass")
    return(Res)
}

univariateMixtureSet <- function(mean, var, proportion, verif = TRUE)
{   
    if (verif)
    {    if (!is_list_numeric_vector(mean))
            return("'mean' parameter must be a list of vectors.\n")
        nStates <- length(mean)
        if (!is_list_numeric_vector(var))
            return("'var' parameter must be a list of vectors.\n")
        if (!is_list_var(var))
            return("all 'var' elements must be positive\n")
        if (!is_list_numeric_vector(proportion))
            return("'proportion' parameter must be a list of vectors.\n")
        if (!is_list_proba_vector(proportion))
            return("The sum of any elements of 'proportion' must be 1\n")
        if (length(var) != nStates)
            return("length of 'var' parameter must be equal to length of 'mean' which is the number of hidden states.\n")
        if (length(proportion) != nStates)
            return("length of 'proportion' parameter must be equal to length of 'mean' which is the number of hidden states.\n")
        
        nMixt <- length(mean[[1]])
        Aux1 <- as.integer(lapply(mean, length)) - nMixt
        Aux2 <- as.integer(lapply(var, length)) - nMixt
        Aux3 <- as.integer(lapply(proportion, length)) - nMixt
        if (sum(abs(Aux1)+abs(Aux2)+abs(Aux3)) != 0)
            return("The number of mixtures must be the same for every hidden states.\n")
     }
    else
    {   nStates <- length(mean)
        Aux <- mean[[1]]
        nMixt <- length(Aux)
     }
    Res <- list(dis="MIXTURE", nStates=nStates, nMixt=nMixt, dimObs = as.integer(1), mean=mean, var=var, proportion=proportion)
    class(Res) <- c("distributionClass", "mixtureUnivariateNormalClass")
    return(Res)
}

multivariateMixtureSet <- function(mean, cov, proportion, verif = TRUE)
{   
     if (verif)
    {    if (!is_list_list_numeric_vector(mean))
            return("'mean' parameter must be a list of list of vectors.\n")
        nStates <- length(mean)
        if (!is_list_list_positive_definite(cov))
            return("'cov' parameter must be a list of list of positive definite matrices.\n")
        if (!is_list_numeric_vector(proportion))
            return("'proportion' parameter must be a list of vectors.\n")
        if (!is_list_proba_vector(proportion))
            return("The sum of any elements of 'proportion' must be 1\n")
        if (length(cov) != nStates)
            return("length of 'cov' parameter must be equal to length of 'mean' which is the number of hidden states.\n")
        if (length(proportion) != nStates)
            return("length of 'proportion' parameter must be equal to length of 'mean' which is the number of hidden states.\n")
        nMixt <- length(mean[[1]])
        Aux1 <- sum(abs(as.integer(lapply(mean, length)) - nMixt))
        #Aux2 <- sum(abs(as.integer(lapply(cov, dim)) - nMixt))
        Aux3 <- sum(abs(as.integer(lapply(proportion, length)) - nMixt))
        if (Aux1+Aux3 != 0)
            return("The number of mixtures must be the same for every hidden states.\n")
        Aux <- mean[[1]]
        dimObs = length(Aux[[1]])
    }
    else
    {   nStates <- length(mean)
        Aux <- mean[[1]]
        nMixt <- length(Aux)
        dimObs <- length(Aux[[1]])
    }
    Res <- list(dis="MIXTURE", nStates=nStates, nMixt=nMixt, dimObs=dimObs, mean=mean, cov=cov, proportion=proportion)
    class(Res) <- c("distributionClass", "mixtureMultivariateNormalClass")
    return(Res)
}
 
distributionSet <- function(dis, ...)
{   
# DEBUT de distributionSet

    if (is.na(match(dis, c('NORMAL', 'MIXTURE', 'DISCRETE'))))
        stop("dis must be in 'NORMAL', 'MIXTURE', 'DISCRETE'\n")
    args <- list(...)
    extras <- match.call(expand.dots = FALSE)$... # pour récupérer les noms
    mcall <- list(as.name("toto"))
    lextras <- length(extras)
    for (i in 1:lextras)
    {   mcall <- c(mcall, args[i])
        names(mcall[i+1]) <- names(extras[i])
    }    
    
    if (dis=="NORMAL")
    {   if ( (lextras == 2) || (lextras == 3) )
        {   if (!any(sapply(args, is.list)))
            {   mcall[[1]] <- as.name("univariateNormalSet")
                nmcall <- names(mcall)
                for (i in 1:length(nmcall))
                    if ( (!is.null(nmcall[i])) && (nmcall[i] !='') )
                        if ( !(nmcall[i] %in% c("mean", "var", "verif")) )
                        {   mess <- sprintf("'%s' parameter unknown", nmcall[i])
                            stop(mess)
                        }
            }    
            else
            {   mcall[[1]] <- as.name("multivariateNormalSet")
                nmcall <- names(mcall)
                for (i in 1:length(nmcall))
                    if ( (!is.null(nmcall[i])) && (nmcall[i] !='') )
                        if ( !(nmcall[i] %in% c("mean", "cov", "verif")) )
                        {   mess <- sprintf("'%s' parameter unknown", nmcall[i])
                            stop(mess)
                        }
            }
            value <- (eval(as.call(mcall)))
            if (is.character(value))
                stop(value)
            else 
                return(value)
        }
        else
            stop("Wrong number of parameters.\n")
    }
    
    if (dis == "MIXTURE")
    {   if ( (lextras == 3) || (lextras == 4) )
        {   if (any(sapply(args, is_list_list)))
            {   mcall[[1]] <- as.name("multivariateMixtureSet")
                nmcall <- names(mcall)
                 for (i in 1:length(nmcall))
                    if ( (!is.null(nmcall[i])) && (nmcall[i] !='') )
                        if (! (nmcall[i] %in% c("mean", "cov", "proportion", "verif")))
                        {   mess <- sprintf("'%s' parameter unknown", nmcall[i])
                            stop(mess)
                        }
                value <- (eval(as.call(mcall)))
                if (is.character(value))
                    stop(value)
                else 
                    return(value)
           }
            else
            {   mcall[[1]] <- as.name("univariateMixtureSet")
                nmcall <- names(mcall)
                for (i in 1:length(nmcall))
                    if ( (!is.null(nmcall[i])) && (nmcall[i] !='') )
                        if (! (nmcall[i] %in% c("mean", "var", "proportion", "verif")))
                        {   mess <- sprintf("'%s' parameter unknown", nmcall[i])
                            stop(mess)
                        }
                value <- (eval(as.call(mcall)))
                if (is.character(value))
                    stop(value)
                else 
                    return(value)
            }
        }
        else
            stop("wrong number of parameters.\n")
    }
    
    if (dis == "DISCRETE")
    {   if  ( (lextras == 1) || (lextras == 2) || (lextras == 3) )
        {   mcall[[1]] <- as.name("discreteSet")
            nmcall <- names(mcall)
            for (i in 1:length(nmcall))
                if ( (!is.null(nmcall[i])) && (nmcall[i] !='') )
                    if (! (nmcall[i] %in% c("proba", "labels", "verif")))
                    {   mess <- sprintf("'%s' parameter unknown", nmcall[i])
                        stop(mess)
                    }
           value <- (eval(as.call(mcall)))
            if (is.character(value))
            stop(value)
            else 
            return(value)
         }
        else
            stop("wrong number of parameters.\n")
    }     
}

print.univariateNormalClass <- function(x, ...)
{   Aux <- cbind(x$mean, x$var)
    Aux <- as.data.frame(Aux, row.names=" ")
    names(Aux) <- c("mean", "var")
    rnames <- as.character(gl(x$nStates, 1, label="State "))
    rownames(Aux) <- rnames
    print.data.frame(Aux, quote=FALSE, right=TRUE)
}

print.multivariateNormalClass <- function(x, ...)
{   
    for (i in 1:x$nStates)
    {   cat(sprintf("  State %d\n", i), sep="")
        Aux <- cbind(x$mean[[i]], x$cov[[i]])
        Aux <- as.data.frame(Aux)
        rnames <- rep("    ", x$dimObs)
        cnames <- c("    ", rnames)
        cnames[1] <- "mean"
        cnames[2] <- "cov matrix"
        names(Aux) <- cnames
        Aux <- as.matrix(Aux)
        row.names(Aux) <- rnames
        #print.matrix(Aux, quote=FALSE, right=TRUE)    
        print(Aux, quote=FALSE, right=TRUE)    
        cat("\n")
    }
}   

print.discreteClass <- function(x, ...)
{   proba <- x$proba
    Aux <- matrix(nrow=x$nStates, ncol=x$nLevels)
    for (i in 1:x$nStates)
        Aux[i, ]<- t(proba[[i]])
    rnames <- as.character(gl(x$nStates, 1, label="State "))
    Aux <- as.data.frame(Aux)
    if (is.null(names(proba[[1]])))
        cnames <- as.character(gl(x$nLevels, 1, label="p"))
    else
        cnames <- names(proba[[1]])
            
    names(Aux) <- cnames
    rownames(Aux) <- rnames
    print.data.frame(Aux, quote=FALSE, right=TRUE)
}

print.mixtureUnivariateNormalClass <- function(x, ...)
{   rnames <- as.character(gl(x$nMixt, 1, label="mixt. "))
    for (i in 1:x$nStates)
    {   cat(sprintf("  State %d\n", i), sep="")
        Aux <- matrix(c(x$mean[[i]], x$var[[i]], x$proportion[[i]]), ncol=3)
        Aux <- as.data.frame(Aux)
        names(Aux) <- c("mean", "var", "prop")
        rownames(Aux) <- rnames
        print.data.frame(Aux, quote=FALSE, right=TRUE)
        cat("\n")
    }
}

print.mixtureMultivariateNormalClass <- function(x, ...)
{
     for (i in 1:x$nStates)
    {   cat(sprintf("  State %d\n", i), sep="")
        for (j in 1:x$nMixt)
        {   cat(sprintf("Mixt. %d\n", j), sep="")
            prop <- c(x$proportion[[i]][j], rep('', x$dimObs -1))
            Aux <- cbind(x$mean[[i]][[j]], x$cov[[i]][[j]], prop)
            colnames(Aux) <- c("mean", "cov", rep(" ", x$dimObs-1), "prop.")
            Aux <- as.data.frame(Aux)
             print.data.frame(Aux, quote=FALSE, right=TRUE, check.names=FALSE)
            cat("\n")
        }
        cat("\n")
    }
}

print.distributionClass <- function(x, ...)
{        
    call <- match.call()
    if (is.null(call$doNotAffiche))
    {   if (x$dis == "NORMAL")
        {    if (x$dimObs==1)
                nomloi <- "univariate gaussian"
             else
                nomloi <- sprintf("%d-D gaussian", x$dimObs)
        }
        if (x$dis == "DISCRETE")
            nomloi <- "discrete"
        if (x$dis == "MIXTURE")
        {    if (is.null(x$dimObs))
                nomloi <- sprintf("mixture of %d gaussian", x$nMixt)
            else
                nomloi <- sprintf("mixture of %d %d-d gaussian", x$nMixt, x$dimObs)
        }
        Model <- sprintf("%d states HMM with %s distributions", x$nStates, nomloi)  
        cat("\nModel:\n", Model, "\n", sep = "")
    }
    cat("\nDistribution parameters:\n", sep="")
    NextMethod("print", x)    
}

HMMSet <- function(initProb, transMat, ...)
{
    args <- list(...)
    extras <- match.call(expand.dots = FALSE)$...
    lextras <- length(extras)
    if (lextras == 0)
        stop("Wrong number of parameters")
 
    if (lextras == 1)
        distribution <- args[[1]]
    else
    {   mcall <- list(1)
        for (i in 2:lextras)
            mcall <- c(mcall, list(1))
        for (i in 1:lextras)
            mcall[i] <- args[i]
        names(mcall) <- names(extras)
        mcall <- c(as.name("distributionSet"), mcall)
        mcall <- as.call(mcall)
        distribution <- eval(mcall)
    }

    if (is.na(match("distributionClass", class(distribution))))
        stop("'distribution' must be a distributionClass object\n")
    
    if (!is.vector(initProb))
        stop('initProb must be a vector\n')
    if (!is.matrix(transMat))
        stop('MatTrans must be a matrix\n')
    nStates <- length(initProb)
    if ( (ncol(transMat) != nStates) || (nrow(transMat) != nStates) || (distribution$nStates != nStates) )
        stop('length of initProb, dim of transMat or dim of distribution parameters do not match\n')
    
    Res <- list(initProb = initProb, transMat = transMat, distribution = distribution)
    class(Res) <- 'HMMClass'
    return(Res)
}


print.HMMClass <- function(x, ...)
{   y <- x$distribution
    if (y$dis == "NORMAL")
    {   if (y$dimObs==1)
        {   nomloi <- "univariate gaussian"
        }
        else
        {   nomloi <- sprintf("%d-d gaussian", y$dimObs)
        }
    }
    if (y$dis == "DISCRETE")
        nomloi <- "discrete"
    if (y$dis == "MIXTURE")
    {   if ( y$dimObs==1)
        {    nomloi <- "gaussian mixture"
        }
        else
        {   nomloi <- sprintf("%d-d gaussian mixture", y$dimObs)
        }
    }
    arg=list(...)
    if (!is.null(arg))
    {   if (is.null(arg$doNotAffiche))
            doNotAffiche <- FALSE
        else
            doNotAffiche <- arg$doNotAffiche
        if (is.null(arg$noHMM))
            noHMM <- FALSE
        else
            noHMM <- arg$noHMM
    }
    else
    {   doNotAffiche <- FALSE
        noHMM <- FALSE
    }
    
    if (!doNotAffiche)
    {
        if (!noHMM)
            Model <- sprintf("%d states HMM with %s distribution", y$nStates, nomloi)     
        else
            Model <- sprintf("%d mixtures of %s distributions", y$nStates, nomloi)
        cat("\nModel:", sep="\n")
        cat("------", sep="\n")
        cat(Model, "\n", sep = "")
    }
    
   cat("\nInitial probabilities:\n", sep="")
        Aux <- t(x$initProb)
        Aux <- as.data.frame(Aux)
        cnames <- as.character(gl(x$distribution$nStates, 1, labels="Pi"))
        names(Aux) <- cnames
        row.names(Aux) <- " "
        print.data.frame(Aux, quote=FALSE, right=TRUE)
      if (! noHMM)
    {      cat("\nTransition matrix:\n", sep="")
            Aux <- as.data.frame(x$transMat)
        cnames <- as.character(gl(x$distribution$nStates, 1, labels="State "))
        names(Aux) <- cnames
        rownames(Aux) <- cnames
        print.data.frame(Aux, quote=FALSE, right=TRUE)
    }
    cat("\nConditionnal distribution parameters:\n", sep="")
    print(x$distribution, quote=FALSE, right=TRUE, doNotAffiche=TRUE)
}

rdiscrete <- function(nSim, proba)
{   probaCum <- proba
    nLevels <- length(proba)
    for (i in 1:nLevels)
        probaCum[i] <- sum(proba[1:i])
    Eps <- runif(nSim)
    Res <- rep(0,nSim)
    for (i in 1:nSim)
    {   k <- 1
        while (Eps[i] > probaCum[k])
            k <- k+1
        Res[i] <- k
    }
    return(Res) 
}

sim <- function(object, nSim, ...)
UseMethod("sim")

sim.univariateNormalClass <- function(object, nSim)
{   value <- NULL
    for (i in 1:object$nStates)
        value <- cbind(value, rnorm(nSim, object$mean[i], sqrt(object$var[i])))
    return(value)
}

sim.multivariateNormalClass <- function(object, nSim)
{   value <- array(dim=c(nSim, object$dimObs, object$nStates))
    for (i in 1:object$nStates)
            value[, , i] <- mvrnorm(nSim, mu=object$mean[[i]], Sigma=object$cov[[i]])
    return(value)
}

trouve_indice <- function(x, probaCum)
{   Aux <- rank(c(x, probaCum))
    return(as.integer(Aux[1]))
}

rmixt <- function(nSim, mean, sd, prop)
{   nMixt <- length(mean)
    YAux <- matrix(ncol=nMixt, nrow=nSim)
    for (i in 1:nMixt)
        YAux[,i] <- rnorm(nSim, mean=mean[i], sd=sd[i])
    probaCum <- rep(0, nMixt)
    for (i in 1:nMixt)
        probaCum[i] <- sum(prop[1:i])
    Eps <- runif(nSim)
    value <- rep(0, nSim)
    for (i in 1:nSim)
        value[i] <- YAux[i, trouve_indice(Eps[i], probaCum)]
    return(value)
}

rmultimixt <- function(nSim, mean, cov, prop)
{   nMixt <- length(mean)
    dimObs <- length(mean[[1]])
    YAux <- array(dim=c(nSim, dimObs, nMixt))
    for (i in 1:nMixt)
        YAux[ , , i] <- mvrnorm(nSim, mu=mean[[i]], Sigma=cov[[i]])
    
    probaCum <- rep(0, nMixt)
    for (i in 1:nMixt)
        probaCum[i] <- sum(prop[1:i])
    Eps <- runif(nSim)
    value <- matrix(0, nSim, dimObs)
    for (n in 1:nSim)
    {   j <-  trouve_indice(Eps[n], probaCum)
        value[n, ] <- YAux[n, , j]
    }
    return(value)
}


dmixt <- function(x, mean, var, prop)
{   nMixt <- length(mean)
    YAux <- matrix(0, length(x), nMixt)
    for (i in 1:nMixt)
        YAux[, i] <- dnorm(x, mean=mean[i], sd=sqrt(var[i]))
    value <- YAux %*% prop
    return(value)
}
    

sim.mixtureUnivariateNormalClass <- function(object, nSim)
{   value <- NULL
    for (i in 1:object$nStates)
        value<- cbind(value, rmixt(nSim, object$mean[[i]], sqrt(object$var[[i]]), object$proportion[[i]]))
    return(value)
}

sim.mixtureMultivariateNormalClass <-  function(object, nSim)
{
    value <- array(dim=c(nSim, object$dimObs, object$nStates))
    for (i in 1:object$nStates)
            value[, , i] <- rmultimixt(nSim, object$mean[[i]], object$cov[[i]], object$proportion[[i]])
    return(value)
}

rdiscrete <- function(nSim, proba)
{   nProba <- length(proba)
    probaCum <- rep(0, nProba)
    for (i in 1:nProba)
        probaCum[i] <- sum(proba[1:i])
    Aux <- runif(nSim)
    value <- rep(0, nSim)
    for (i in 1:nSim)
        value[i] <- trouve_indice(Aux[i], probaCum)
    return(value)
}

sim.discreteClass <- function(object, nSim)
{   
    value <- NULL
    for (i in 1:object$nStates)
        value <- cbind(value, rdiscrete(nSim, object$proba[[i]]))
    
    return(value)
}

sim.distributionClass <- function(object, nSim)
{
    return(NextMethod("sim", object))
}

sim.markovChainClass <- function(object, nSim)
{
    probaCum <- rep(0, object$nStates)
    Aux <- runif(nSim)
    value <- as.integer(rep(0, nSim))
    for (i in 1:object$nStates)
        probaCum[i] <- sum(object$initProb[1:i])
    value[1] <- trouve_indice(Aux[1], probaCum)
    probaCum <- matrix(0, nrow=object$nStates, ncol=object$nStates)
    for (i in 1:object$nStates)
        for (j in 1:object$nStates)
            probaCum[i,j] <- sum(object$transMat[i,1:j])
  
    for (t in 2:nSim)
        value[t] <- trouve_indice(Aux[t], probaCum[value[t-1],])
    
    return(as.integer(value))
}

sim.HMMClass <- function(object, nSim)
{
    mc <- list(nStates = object$distribution$nStates, initProb=object$initProb, transMat=object$transMat)
    class(mc) <- "markovChainClass" 
    Eps <- sim(mc, nSim)
    obs <- sim(object$distribution, nSim)
    if (all(is.na(match(c("multivariateNormalClass","mixtureMultivariateNormalClass"), class(object$distribution)))))
    {   value <- rep(0, nSim)
        for (t in 1:nSim)
            value[t] <- obs[t,Eps[t]]
        if (!is.na(match("discreteClass", class(object$distribution))))
        {   value <- as.factor(value)
            levels(value) <- names(object$distribution$proba[[1]])
        }
    }
    else
    {   value <- matrix(0, nrow=nSim, ncol=object$distribution$dimObs)
        for (t in 1:nSim)
            value[t, ] <- obs[t, ,Eps[t]]
     }
    
    return(list(obs=value, states=Eps))    
    
}

HMMSim <- function(nSim, HMM)
{
    if (nSim %% 1 != 0)
        stop("nSim must be a positive integer\n")
    if (nSim < 0)
        stop('nSim must be a positive integer')
    if (class(HMM) != "HMMClass")
        stop("HMM be be a HMMClass object. See HMMSet")
            
    return(sim(HMM, nSim))
}

incr <- function(x)
{
    x<-x+1
    
}

HMMKMeans <- function(obs, nClass)
{  if (is.list(obs))
    {   x <- obs[[1]]
        dimObs <- dim(x)[2]
        if (!is.null(dimObs))
        {   for (i in 2:length(obs))
                x <- rbind(x, obs[[i]])
        }
        else
        {   for (i in 2:length(obs))
                x <- c(x, obs[[i]])
            dimObs <- 1
        }
    }
    else
    {    x <- obs
        dimObs <- dim(x)[2]
        if (is.null(dimObs))
            dimObs <- 1
    }
    z <- kmeans(x, centers=nClass)
     if (dimObs > 1)
    {   CovAux <- matrix(0, dimObs, dimObs)
        Cov <- list(rep(0, nClass))
         for (i in 1:nClass)
        {   Cov[[i]]  <- CovAux
        }
    }
    else
    {
        Var <- rep(0, nClass)
    }
    if (dimObs > 1)
        Mean <- list(rep(0, nClass))
    else
        Mean <- rep(0, nClass)
    for (i in 1:nClass)
    {   if (dimObs > 1)
        {   Cov[[i]] <- cov(x[z$cluster==i, ])
            Mean[[i]] <- z$centers[i,]
        }
        else
        {   Var[i] <- var(x[z$cluster==i])
            Mean[i] <- z$center[i]
        }
        
    }
    TT <- length(z$cluster)
    transMat <- matrix(0, nClass, nClass)
    for (tt in 2:TT)
        transMat[z$cluster[tt-1], z$cluster[tt]]<-transMat[z$cluster[tt-1], z$cluster[tt]]+1
    for (j in 1:nClass)
        transMat[j,] <- transMat[j, ]/sum(transMat[j,])
    initProb <- rep(1/nClass, nClass)
    
    if (dimObs > 1)
    {    
#        Res <- list(Mean=Mean, Cov=Cov, transMat=transMat, initProb=initProb)
        Res <- HMMSet(dis='NORMAL', transMat=transMat, initProb=initProb, mean=Mean, cov = Cov)
    }
    else
    {   
#    Res <- list(Mean=Mean, Var = Var, transMat=transMat, initProb=initProb)
        Res <- HMMSet(dis='NORMAL', transMat=transMat, initProb=initProb, mean=Mean, var = Var)

    }
    return(Res)
}

BaumWelch<-function(paramHMM, obs, paramAlgo)
{   
    paramAlgo <- setStorageMode(paramAlgo)
        
    maListe <- TransformeList(paramHMM$dis, obs)
    paramHMM$nLevels <- maListe$nLevels
    paramHMM <- setStorageMode(paramHMM)
    
    Res1 <- .Call("RBaumWelch", paramHMM, maListe$Zt, paramAlgo)
    
    if (paramHMM$dis=="NORMAL") 
    {   if (paramHMM$dimObs == 1)
            distribution <- distributionSet(dis="NORMAL", mean=Res1[[3]], var=Res1[[4]], verif=FALSE)
        else
            distribution <- distributionSet(dis="NORMAL", mean=Res1[[3]], cov=Res1[[4]], verif=FALSE)
        Autres <- Res1[[5]]
    }
        
    if (paramHMM$dis == "DISCRETE")
    {   distribution <- distributionSet(dis="DISCRETE", proba=Res1[[3]], labels=maListe$labels, verif=FALSE)
        Autres <- Res1[[4]]
    }
    
    if (paramHMM$dis == "MIXTURE")
    {   if (paramHMM$dimObs==1)
            distribution <- distributionSet(dis="MIXTURE", mean=Res1[[3]], var=Res1[[4]], proportion=Res1[[5]], verif=FALSE)
        else
            distribution <- distributionSet(dis="MIXTURE", mean=Res1[[3]], cov=Res1[[4]], proportion=Res1[[5]], verif=FALSE)
        Autres <- Res1[[6]]
    }
 
    Res2 <- HMMSet(initProb=Res1[[1]], transMat=Res1[[2]], distribution=distribution)
        
    Res <- list(HMM=Res2, LLH=Autres[1], BIC=Autres[2], nIter=as.integer(Autres[3]), relVariation=Autres[4])
    class(Res) <- "HMMFitClass"
    return(Res)
}

HMMFit <- function(obs, dis="NORMAL", nStates = 2,...)
{   noHMM=FALSE # not implemented yet
    nMixt <- NULL
    control <- NULL
    args <- list(...)
    extras <- match.call(expand.dots=FALSE)$...
    lextras <- 0
    if (!is.null(extras))
        lextras <- length(extras)
    
    if (dis == "MIXTURE")
    {    if (lextras == 0)
            stop("HMMFit with gaussian MIXTURE distributions needs a 'nMixt' parameter.\n")
        if (lextras > 2) 
            stop("too many parameters.\n")
        if (is.list(args[[1]]))
        {   control <- args[[1]]
            if (lextras != 2)
                stop("HMMFit with gaussian MIXTURE distributions needs an integer ('nMixt') parameter.\n")
            else
                nMix <- args[[2]]
        }
        else
        {   nMixt <- args[[1]]
            if (lextras == 2)
                control <- args[[2]]
        }
        
        nextras <- names(extras)
        for (nn in nextras)
            if ( !(nn %in% c("nMixt", "control")) && (nn != '') )
                {   nerror <- sprintf("unknown '%s' parameter.\n", nn)
                    stop(nerror)
                }
        if ( !(is.null(control)) && (!is.list(control)) )
            stop("'control' parameter must be a list.\n")
    }
    
    if (dis != "MIXTURE")
    {   if (lextras > 1)
            stop("too many parameters.\n")
        
        if (lextras == 1)
        {   control <- args[[1]]    
            nnames <- names(extras[[1]])
            if ( !is.null(nnames) )
                if ( (nnames != '') && (nnames != 'control') )
                {   nerror <- sprintf("unknown '%s' parameter.\n", nnames)
                    stop(nerror)
                }
            if (!is.list(control))
                stop("'control' parameter must be a list.\n")
        }
    }
        
    if (is.null(control))
        control <- list(init="RANDOM", iter=500, tol=1e-6, verbose=0, nInit=5, nIterInit=5, initPoint=NULL)
    
    nnames <- names(control)
    lnames <- length(nnames)
    rnames <- c("init", "iter", "tol", "verbose", "nInit", "nIterInit", "initPoint")
    for (i in 1:lnames)
    {   if ( (is.null(nnames[i])) || (nnames[i] == '') )
            names(control)[i] <- rnames[i]
        else
            if (is.na(match(nnames[i], rnames)))
            {   nerror <- sprintf("unknown '%s' element of 'control' parameter", nnames[i])
                stop(nerror)
            }
    }
       
    if (is.null(control$init))
        control$init <- "RANDOM"
    if (is.null(control$iter))
        control$iter <- 500
    if (is.null(control$tol))
        control$tol <- 1e-6
    if (is.null(control$verbose))
        control$verbose <- 0
    if (is.null(control$nInit))
        control$nInit <- 5
    if (is.null(control$nIterInit))
        control$nIterInit <- 5
   
    
    
    if (is.na(match(dis, c('NORMAL', 'MIXTURE', 'DISCRETE'))))
        stop("'dis' parameter must be in 'NORMAL', 'MIXTURE', 'DISCRETE'\n")
    
    if (!is.numeric(nStates))
        stop('nStates must be a positive number\n')
    nStates <- as.integer(nStates)
    if (nStates <= 0)
        stop("nStates must be a positive number\n")
        
    if (is.null(nMixt))
        nMixt <- as.integer(0)
    else
    {   if (!is.numeric(nMixt))
            stop('nMixt must be a a positive integer\n')
        nMixt <- as.integer(nMixt)
        if (nMixt < 0)
            stop("nMixt must be a positive integer\n")
    }
    
    if ( (dis == 'MIXTURE') && (nMixt < 2) )
    {   nMixt <- as.integer(2)
        warning("gaussian MIXTURE conditionnal distribution with nMixt < 2. nMixt = 2 used instead.\n")
    }
     
    if (!is.null(control$initPoint))
    {    if (class(control$initPoint) != "HMMClass")
            stop("control$initPoint class must be HMMClass. See HMMSet\n")
         control$init <- "USER"
         initPoint <- control$initPoint
    }
    else
    {   initPoint <- NULL
        if (control$init == 'USER')
        {   control$init <- 'RANDOM'
            warning("'USER' initialization and no initPoint. 'RANDOM' initialization used instead.\n")
        }
    }
    
    if ( (control$init == 'KMEANS') && (dis != 'NORMAL') )
    {   warning("'KMEANS' initialization only for 'NORMAL' conditionnal distribution. 'RANDOM' init used instead.\n")
        control$init <- 'RANDOM'
    }
        
    if (is.na(match(control$init, c("RANDOM", "KMEANS", "USER"))))
        stop("control$init must be in 'RANDOM', 'KMEANS', 'USER'")
        
    init <- control$init
    
    if (!is.numeric(control$iter))
        stop('control$iter must be a positive integer\n')
    iter <- as.integer(control$iter)
    if (iter <= 0)
        stop('control$iter must be a positive integer\n')
    
    if (!is.numeric(control$tol))
        stop('control$tol must be a positive real\n')
    tol <- as.double(control$tol)
    if (tol < 0)
        stop('control$tol must be a positive double\n')
   
    if (!is.numeric(control$verbose))
        stop('control$verbose must be 0 or 1\n')
    verbose <- as.integer(control$verbose)
    if (is.na(match(verbose, c(0,1,2))))
        stop('control$verbose must be 0 or 1\n')
        
    if (!is.numeric(control$nInit))
        stop('control$nInit must be a positive integer\n')
    nInit <- as.integer(control$nInit)
    if  (nInit < 0)
        stop('control$nInit must be a positive integer\n')
   
    if (!is.numeric(control$nIterInit))
        stop('control$nIterInit must be a positive integer\n')
    nIterInit <- as.integer(control$nIterInit)
    if  (nIterInit < 0)
        stop('control$nIterInit must be a positive integer\n')
     
    if (control$init == 'KMEANS')
    {   initPoint <- HMMKMeans(obs, nStates)
        init='USER'
    }
    
    if (!is.logical(noHMM))
    {   stop('noHMM must be a boolean\n')
    }
    
    paramAlgo <- list(init=init, iter=iter, tol=tol, verbose=verbose, nInit=nInit, nIterInit=nIterInit, initPoint=initPoint)

   class(paramAlgo) <- "paramAlgoBW"
    if (is.list(obs))
        dimObs <- ncol(as.matrix(obs[[1]]))
    else
        dimObs <- ncol(as.matrix(obs))
    
    noHMMInt <- as.integer(noHMM)        
    paramHMM <- list(nStates=nStates, dimObs=dimObs, nMixt = nMixt, nLevels = 0, noHMM=noHMMInt, dis=dis) 
    class(paramHMM) <- "paramHMM"
   
    Res<-BaumWelch(paramHMM, obs, paramAlgo)
    Res$call <- match.call()
    Res$noHMM <- noHMM
    return(Res)
}

print.HMMFitClass <- function(x, ...)
{
   # Call and Model:
    cat("\nCall:", sep="\n")
    cat("----", sep="\n")
    cat(deparse(x$call), "\n", sep = "")
    y <- x$HMM$distribution
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
    noHMM <- as.logical(x$noHMM)
    
    if (!noHMM)
        Model <- sprintf("%d states HMM with %s distribution", y$nStates, nomloi)     
    else
        Model <- sprintf("%d mixtures of %s distributions", y$nStates, nomloi)
    cat("\nModel:", sep="\n")
    cat("------", sep="\n")
    cat(Model, "\n", sep = "")
    # Algorithm
    
    cat("\nBaum-Welch algorithm status:", sep="\n")
    cat("----------------------------", sep="\n")
    cat(sprintf("Number of iterations : %d\n", x$nIter), sep="")
    cat(sprintf("Last relative variation of LLH function: %f\n", x$relVariation), sep="")
   
   # Estimation
   cat("\nEstimation:", sep="\n")
   cat("-----------", sep="\n")
   print(x$HMM, doNotAffiche=TRUE, noHMM=noHMM)
   cat("\nLog-likelihood: ",format(round(x$LLH, 2)), "\n", sep="")
   cat("BIC criterium: ",format(round(x$BIC, 2)), "\n", sep="")
}

forwardbackward<-function(HMM, obs)
{   if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM
        
    HMM <- setStorageMode(HMM)
    maListe <- TransformeList(HMM$distribution$dis, obs)
    Res1 <- .Call("Rforwardbackward", HMM, maListe$Zt)
    names(Res1) <- c("Alpha", "Beta", "Gamma", "Xsi", "Rho", "LLH")
    if (!is.list(obs))
    {   Res1$Alpha <- t(Res1$Alpha[[1]])
        Res1$Beta <- t(Res1$Beta[[1]])
        Res1$Gamma <- t(Res1$Gamma[[1]])
        Res1$Xsi <- t(Res1$Xsi[[1]])
        Res1$Rho <- Res1$Rho[[1]]
        Res1$LLH <- Res1$LLH[[1]]
    }
    else
    {   for (n in 1:length(obs))
        {   Res1$Alpha[[n]] <- t(Res1$Alpha[[n]])
            Res1$Beta[[n]] <- t(Res1$Beta[[n]])
            Res1$Gamma[[n]] <- t(Res1$Gamma[[n]])
            Res1$Xsi[[n]] <- t(Res1$Xsi[[n]])
         }
    }   
    return(Res1)
}

inc <- function(x)
{
    if (is.list(x))
    {   for (i in 1:length(x))
            x[[i]] <- x[[i]]+1
    }
    else
        x <- x + 1
    return(x)
}

viterbi<-function(HMM, obs)
{   if ( ( class(HMM) != "HMMFitClass" ) && (class(HMM) != "HMMClass") )
        stop("class(HMM) must be 'HMMClass' or 'HMMFitClass'\n")
    
    if (class(HMM) == "HMMFitClass")
        HMM <- HMM$HMM

    
    HMM <- setStorageMode(HMM)
   
    maListe <- TransformeList(HMM$distribution$dis, obs)
    
    Res1 <- .Call("RViterbi", HMM, maListe$Zt)
    names(Res1) <- c("states", "logViterbiScore")
    Res1$states <- inc(Res1$states)
    Aux <- forwardbackward(HMM, obs)
    LLH <- Aux$LLH
 
     if (is.list(obs))
    {    probSeq <- sapply(Res1$logViterbiScore, sum) - sapply(LLH, sum)
        Res<-list(states=Res1$states, logViterbiScore = Res1$logViterbiScore, logProbSeq=as.list(probSeq))
    }
    else
        Res <- list(states=Res1$states[[1]], logViterbiScore = Res1$logViterbiScore[[1]], logProbSeq=Res1$logViterbiScore[[1]]-LLH)
    class(Res) <- "viterbiClass"
    return(Res)
}

graphicDiag <- function(object, vit, obs, color="green", ...) UseMethod("graphicDiag")

graphicDiag.default <- function(object, vit, obs, color="green", ...)
{   cat("Not yet implemented for this kind of distribution\n")
    return(invisible(0))
}

graphicDiag.HMMFitClass <- function(object, vit, obs, color="green", ...)
{   graphicDiag(object$HMM$distribution, vit, obs, color)
}

graphicDiag.HMMClass <- function(object, vit, obs, color="green", ...)
{   graphicDiag(object$distribution, vit, obs, color)
}

graphicDiag.distributionClass <- function(object, vit, obs, color="green", ...)
{
    if (is.na(match(color, colours())))
        stop(sprintf("unknown %s color\n", color))
      #x11()
    if (is.list(obs))
    {   if (!is.list(vit$states))
            return("incompatibilty beetween 'vit' and 'obs' parameters\n")
        if (length(vit$states) != length(obs))
            return("incompatibilty beetween 'vit' and 'obs' parameters\n")
        states <- NULL
        Aux <- NULL
        for (n in 1:length(vit$states))
        {   states <- c(states, vit$states[[n]])
            Aux <- c(Aux, obs[[n]])
        }
        obs <- Aux
        vit$states <- states
    }   
      
    NextMethod(generic="graphicDiag", object=object)
}

graphicDiag.univariateNormalClass <- function(object, vit, obs, color="green")
{  
        
    #x11()
    par(bg="white")
    split.screen(c(object$nStates,1))
    for (i in 1:object$nStates)
    {   screen(i)
        y <- obs[vit$states==i]
        m <- object$mean[i]
        sig <- sqrt(object$var[i])
        from <- m-3*sig
        to <- m+ 3*sig
        x0 <- seq(from, to, (to-from)/512)
        y0 <- dnorm(x0, mean=m, sd=sig)
        z <- density(y, from=from, to=to)
        if (max(z$y) > max(y0))
        {   plot(z$x, z$y, col=color, type='l', xlab="", ylab="Density", lwd=2)
            lines(x0, y0)
        }
        else
        {   plot(x0, y0, type='l', xlab="", ylab="Density")
                lines(z$x, z$y, col=color, lwd=2)    
        }
        sub <- sprintf("Estimated Density (%s) - Normal density with estimated parameters (black)\n", color)
        titre <- sprintf("Graphic diagnostic for state %d", i)
        title(main=titre, sub=sub)
    } 
    invisible(close.screen(all.screens = TRUE))
}

graphicDiag.mixtureUnivariateNormalClass <- function(object, vit, obs, color="green")
{  
 #   x11()
 par(bg="white")
    split.screen(c(object$nStates,1))
    for (i in 1:object$nStates)
    {   screen(i)
        y <- obs[vit$states==i]
        m <- object$mean[[i]]
        var <- object$var[[i]]
        prop <- object$proportion[[i]]
        xbarre <- mean(y)
        sig <- sd(y)
        from <- xbarre-3*sig
        to <- xbarre + 3*sig
        x0 <- seq(from, to, (to-from)/512)
        y0 <- dmixt(x0, mean=m, var=var, prop=prop)
        z <- density(y, from=from, to=to)
        if (max(z$y) > max(y0))
        {   plot(z$x, z$y, col=color, type='l', xlab="", ylab="Density", lwd=2)
            lines(x0, y0)
        }
        else
        {   plot(x0, y0, type='l', xlab="", ylab="Density")
            lines(z$x, z$y, col=color, lwd=2)    
        }
        sub <- sprintf("Estimated Density (%s) - Mixture of univariate normal density with estimated parameters (black)\n", color)
        titre <- sprintf("Graphic diagnostic for state %d", i)
        title(main=titre, sub=sub)
    } 
    invisible(close.screen(all.screens = TRUE))
}

graphicDiag.discreteClass <- function(object, vit, obs, color="green")
{  
    par(bg="white")
    split.screen(c(object$nStates,1))
    for (i in 1:object$nStates)
    {   screen(i)
        y <- obs[vit$states==i]
        ny <- length(y)
        y <- table(y) / ny
        y0 <- object$proba[[i]]
            
        if (max(y) > max(y0))
        {   plot(y, col=color, xlab="", ylab="Probability", lwd=2)
            lines(y0, type="p", lwd=3)
        }
        else
        {   plot(y, col=color, xlab="", ylab="Probability", lwd=2, type="n")
            lines(y0, type='p', lwd=3)
            lines(y, col=color, lwd=2, type="h")    
        }
            
        sub <- sprintf("Frequency (%s) - Estimated parameters (black)\n", color)
        titre <- sprintf("Graphic diagnostic for state %d", i)
        title(main=titre, sub=sub) 
    }
    invisible(close.screen(all.screens = TRUE))
}

HMMGraphicDiag <- function(vit, HMM, obs, color="green")
{   if (class(vit) !="viterbiClass")
        stop("vit must be a viterbiClass object. See viterbi\n")
    if ( (class(HMM) != "HMMClass") && (class(HMM) != "HMMFitClass") && is.na(match("distributionClass", class(HMM))) )   
        stop("distribution must be a distributionClass, HMMClass or a HMMFitClass object\n")
    
    graphicDiag(HMM, vit, obs, color)
    
}

PlotSerie <- function(object, vit, obs, color="black", ...) UseMethod("PlotSerie")
PlotSerie.default <- function(object, vit, obs, color="black", ...)
{   cat("Not yet implemented for this kind of distribution\n")
    return(invisible(0))
}

plotSerie.distributionClass <- function(object, vit, obs, color="black", ...)
{
    if (is.na(match(color, colours())))
        stop(sprintf("unknown %s color\n", color))
      #x11()
    if (is.list(obs))
    {   if (!is.list(vit$states))
            return("incompatibilty beetween 'vit' and 'obs' parameters\n")
        if (length(vit$states) != length(obs))
            return("incompatibilty beetween 'vit' and 'obs' parameters\n")
        states <- NULL
        Aux <- NULL
        for (n in 1:length(vit$states))
        {   states <- c(states, vit$states[[n]])
            Aux <- c(Aux, obs[[n]])
        }
        obs <- Aux
        vit$states <- states
    }   
      
    NextMethod(generic="plotSerie", object=object)
}



HMMPlotSerie <- function(obs, states, dis = "NORMAL", color="green")
{   
    if ( (class(states) !="viterbiClass") && (!is_numeric_vector(states)) && (!is_list_numeric_vector(states)))
        stop("vit must be a viterbiClass object, a numeric vector or a list of numeric vector.\n")
    if (class(states) == "viterbiClass") 
        states <- states$states
        
    if (is.list(obs))
    {   if (!is.list(states))
            stop("incompatibilty beetween 'states' and 'obs' parameters\n")
        if (length(states) != length(obs))
            stop("incompatibilty beetween 'states' and 'obs' parameters\n")
        statesN <- NULL
        Aux <- NULL
        for (n in 1:length(states))
        {   statesN <- c(statesN, states[[n]])
            Aux <- c(Aux, obs[[n]])
        }
        states <- statesN
    }
    else
    {   Aux <- obs
    }
      
    if (!is.null(dim(obs)))
        stop("Not yet implemented for this kind of distribution\n")
   
    nStates <- max(states)
    #x11()
    par(bg="white")
    split.screen(c(nStates,1))
    for (i in 1:nStates)
    {   screen(i)
        z <- TransformeList(dis, Aux)
        y <- (z$Zt[[1]]+1)*(states==i)
        if (dis != 'DISCRETE')
            plot(y, col=color, type='l', xlab="", ylab="Serie", lwd=1)
        else
            plot(y, col=color, type='l', xlab="", ylab="Serie", lwd=1)
        titre <- sprintf("Serie for state %d", i)
        title(main=titre)
    } 
    invisible(close.screen(all.screens = TRUE))
    
}
