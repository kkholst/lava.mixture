frobnorm <- function(x,y=0,...) {
    sum((x-y)^2)^.5
}

###{{{ mixture

##' Estimate mixture latent variable model
##'
##' Estimate parameters in a mixture of latent variable models via the EM algorithm.
##' @title Estimate mixture latent variable model.
##' @param x List of \code{lvm} objects. If only a single \code{lvm}
##' object is given, then a \code{k}-mixture of this model is fitted
##' (free parameters varying between mixture components).
##' @param data \code{data.frame}
##' @param k Number of mixture components
##' @param control Optimization parameters (see details)
##' @param type Type of EM algorithm (standard, classification, stochastic)
##' @param ... Additional arguments parsed to lower-level functions
##' @author Klaus K. Holst
##' @details
##' The performance of the EM algorithm can be tuned via the \code{control}
##' argument, a list where a subset of the following members can be altered:
##' 
##' \describe{
##' \item{start}{Optional starting values}
##' \item{nstart}{Evaluate \code{nstart} different starting values and run the
##'  EM-algorithm on the parameters with largest likelihood}
##' \item{tol}{Convergence tolerance of the EM-algorithm.  The algorithm is
##'  stopped when the absolute change in likelihood and parameter (2-norm)
##'  between successive iterations is less than \code{tol}}
##' \item{iter.max}{Maximum number of iterations of the EM-algorithm}
##' \item{gamma}{Scale-down (i.e. number between 0 and 1) of the step-size
##'  of the Newton-Raphson algorithm in the M-step}
##' \item{trace}{Trace information on the EM-algorithm is printed on every
##' \code{trace}th iteration}
##'}
##'
##' Note that the algorithm can be aborted any time (C-c) and still be saved
##' (via on.exit call).
##' @seealso \code{mvnmix}
##' @keywords models, regression
##' @export
##' @examples
##' \donttest{
##' set.seed(1)
##' m0 <- lvm(list(y~x+z,x~z))
##' distribution(m0,~z) <- binomial.lvm()
##' d <- sim(m0,500,p=c("y<-z"=2,"y<-x"=1))
##' 
##' ## unmeasured confounder example
##' m <- baptize(lvm(y~x));
##' covariance(m,~x) <- "v"
##' intercept(m,~x+y) <- NA
##' 
##' M <- mixture(m,k=2,data=d,control=list(trace=1,tol=1e-6))
##' summary(M)
##' lm(y~x,d)
##' ## True slope := 1
##' }
mixture <- function(x, data, k=length(x), control=list(), type=c("standard","CEM","SEM"),...) {    
    MODEL <- "normal"    
    type <- tolower(type[1])
    if (type[1]!="standard") {
        return(mixture0(x,data=data,k=k,control=control,type=type,...))
    }   


    optim <- list(start=NULL,
                  startbounds=c(-2,2),
                  startmean=FALSE,
                  nstart=1,
                  prob=NULL,
                  delta=1e-2,
                  constrain=TRUE,
                  stopc=2,
                  lbound=1e-9,
                  trace=1,
                  stabil=TRUE,
                  gamma=1,
                  gamma2=1,
                  newton=10,
                  lambda=0 # Stabilizing factor (avoid singularities of I)
                  )

    
    if (!missing(control))
        optim[names(control)] <- control
    if ("iter.max"%in%names(optim)) optim$maxiter <- optim$iter.max

    sqem.idx <- match(c("K","method","square","step.min0","step.max0","mstep",
                            "objfn.inc","kr","tol","maxiter","trace"),
                          names(optim)) 
    sqem.control <- optim[na.omit(sqem.idx)]

        
    if (k==1) {
        if (is.list(x))
            res <- estimate(x[[1]],data,...)
        else
            res <- estimate(x,data,...)
        return(res)
    }
    if (class(x)[1]=="lvm") {
        index(x) <- reindex(x,zeroones=TRUE,deriv=TRUE)
        x <- rep(list(x),k)
    }  

    
    mg <- multigroup(x,rep(list(data),k),fix=FALSE)
    ## Bounds on variance parameters
    npar <- with(mg, npar+npar.mean)
    parpos <- modelPar(mg,1:npar)$p  
    lower <- rep(-Inf, mg$npar);
    offdiagpos <- c()
    varpos <- c()
    for (i in 1:k) {
        vpos <- sapply(mg$parlist[[i]][variances(mg$lvm[[i]])], function(y) as.numeric(substr(y,2,nchar(y))))
        offpos <- sapply(mg$parlist[[i]][offdiags(mg$lvm[[i]])], function(y) as.numeric(substr(y,2,nchar(y))))
        varpos <- c(varpos, vpos)
        offdiagpos <- c(offdiagpos,offpos)
        if (length(vpos)>0)
            lower[vpos] <- optim$lbound  ## Setup optimization constraints
    }
    lower <- c(rep(-Inf,mg$npar.mean), lower)
    constrained <- which(is.finite(lower))
    if (!any(constrained)) optim$constrain <- FALSE

    mymodel <- list(multigroup=mg,k=k,data=data); class(mymodel) <- "lvm.mixture"

    if (is.null(optim$start)) {
        constrLogLikS <- function(p) {      
            if (optim$constrain) {
                p[constrained] <- exp(p[constrained])
            }
            -logLik(mymodel,p=p,rep(1/k,k))
        }

        start <- runif(npar,optim$startbounds[1],optim$startbounds[2]);
        if (length(offdiagpos)>0)
            start[mg$npar.mean + offdiagpos] <- 0
        if (optim$nstart>1) {
            myll <- constrLogLikS(start)
            for (i in 1:optim$nstart) {
                newstart <- runif(npar,optim$startbounds[1],optim$startbounds[2]);
                newmyll <- constrLogLikS(newstart)
                if (newmyll<myll) {
                    start <- newstart
                }
            }
        }##    start <- optim(1:4,constrLogLikS,method="SANN",control=list(maxit=50))
        optim$start <- start
    }
    
    if (is.null(optim$prob))
        optim$prob <- rep(1/k,k-1)
    thetacur <- optim$start
    probcur <- with(optim, c(prob,1-sum(prob)))
    probs <- rbind(probcur);

    thetas <- rbind(thetacur)
    if (optim$constrain) {
        thetas[constrained] <- exp(thetas[constrained])
    }
    
    gamma <- t(rmultinom(nrow(data),1,probs))
    newgamma <- gamma
    ##  gammas <- list()
    curloglik <- logLik(mymodel,p=thetacur,prob=probcur,model=MODEL)
    vals <- c(curloglik)
    i <- count <- 0
    member <- rep(1,nrow(data))
    E <- dloglik <- Inf
    

    ## constrLogLik <- function(p,prob) {      
    ##   if (optim$constrain) {
    ##     p[constrained] <- exp(p[constrained])
    ##   }
    ##   logLik(mymodel,p=p,prob=prob)
    ## }
    ## EM algorithm:
    myObj <- function(p) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- modelPar(mg,p)$p
        ##    save(p,file="p.rda")
        ##      lf <- sapply(1:k, function(i) logLik(mg$lvm[[i]],p=myp[[i]],data=data,indiv=TRUE))
        ##    print(lf)
        ##      ff <- exp(lf)
        ff <- sapply(1:k, function(j) logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL))
        return(-sum(gamma*ff))
        ## Previous:
        ##      ff <- sapply(1:k, function(j) exp(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)))
        ##      return(-sum(log(rowSums(gamma*ff))))
    }
    myGrad <- function(p) {
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- modelPar(mg,p)$p
        D <- lapply(1:k, function(j) gamma[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL))
        D0 <- matrix(0,nrow(data),length(p))
        for (j in 1:k) D0[,parpos[[j]]] <- D0[,parpos[[j]]]+D[[j]]
        S <- -colSums(D0)
        if (optim$constrain) {
            S[constrained] <- S[constrained]*p[constrained]
        }
        return(S)
        ## Previous:
        ##      ff <- sapply(1:k, function(j) exp(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL)))
        ##      gammaff <- gamma*ff
        ##      f0 <- rowSums(gammaff)
        ##      D <- lapply(1:k, function(j) 1/f0*(gammaff)[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data,model=MODEL))
        ##      D0 <- matrix(0,nrow(data),length(p))
        ##      for (k in 1:k) D0[,parpos[[k]]] <- D0[,parpos[[k]]]+D[[k]]
        ##      -colSums(D0)
    }
    myInformation <- function(p) {
        p0 <- p
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- modelPar(mg,p)$p    
        I <- lapply(1:k, function(j) probcur[j]*information(mg$lvm[[j]],p=myp[[j]],n=nrow(data),data=data,model=MODEL))
        I0 <- matrix(0,length(p),length(p))
        for (j in 1:k) {
            I0[parpos[[j]],parpos[[j]]] <- I0[parpos[[j]],parpos[[j]]] + I[[j]]
        }
        if (optim$constrain) {
            I0[constrained,-constrained] <- apply(I0[constrained,-constrained,drop=FALSE],2,function(x) x*p[constrained]);
            I0[-constrained,constrained] <- t(I0[constrained,-constrained])
            D <- -myGrad(p0)  
            if (length(constrained)==1)
                I0[constrained,constrained] <- I0[constrained,constrained]*outer(p[constrained],p[constrained]) + D[constrained]
            else
                I0[constrained,constrained] <- I0[constrained,constrained]*outer(p[constrained],p[constrained]) + diag(D[constrained])
        }
        return(I0)
    }
    Scoring <- function(p,gamma) {
        p.orig <- p
        if (optim$constrain) {
            p[constrained] <- exp(p[constrained])
        }
        myp <- modelPar(mg,p)$p
        D <- lapply(1:k, function(j) gamma[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE,model=MODEL))
        D0 <- matrix(0,nrow(data),length(p))
        for (j in 1:k) D0[,parpos[[j]]] <- D0[,parpos[[j]]]+D[[j]]
        S <- colSums(D0)
        if (optim$constrain) {
            S[constrained] <- S[constrained]*p[constrained]
        }
        I <- lapply(1:k, function(j) probcur[j]*information(mg$lvm[[j]],p=myp[[j]],n=nrow(data),data=data,model=MODEL
                                                            ))
        I0 <- matrix(0,length(p),length(p))
        for (j in 1:k) {
            I0[parpos[[j]],parpos[[j]]] <- I0[parpos[[j]],parpos[[j]]] + I[[j]]
        }
        if (optim$constrain) {
            I0[constrained,-constrained] <- apply(I0[constrained,-constrained,drop=FALSE],2,function(x) x*p[constrained]);
            I0[-constrained,constrained] <- t(I0[constrained,-constrained])
            if (length(constrained)==1)
                I0[constrained,constrained] <- I0[constrained,constrained]*p[constrained]^2 + S[constrained]
            else
                I0[constrained,constrained] <- I0[constrained,constrained]*outer(p[constrained],p[constrained]) + diag(S[constrained])
        }
        ##    print(paste(S,collapse=","))
        if (optim$stabil) {
            ##      I0 <- I0+S%*%t(S)
            if (optim$lambda>0)
                sigma <- optim$lambda
            else
                sigma <- (t(S)%*%S)[1]^0.5
            I0 <- I0+optim$gamma2*(sigma)*diag(nrow(I0))
        }       
        p.orig + optim$gamma*Inverse(I0)%*%S
        ##   p.orig + Inverse(I0+optim$lambda*diag(nrow(I0)))%*%S
    }
    ## env <- new.env()
    ## assign("mg",mg,env)
    ## assign("k",k,env)
    ## assign("data",data,env)
    ## assign("MODEL",MODEL,env)
    ## assign("optim",optim,env)
    ## assign("parpos",parpos,env)
    ## assign("constrained",constrained,env)

    mytheta <- thetacur
    if (optim$constrain) {
        mytheta <- thetacur
        mytheta[constrained] <- exp(mytheta[constrained])
    }

    p <- c(thetacur,probcur)
    EMstep <- function(p,all=FALSE) {        
        thetacur <- p[seq_along(thetacur)]
        thetacur0 <- thetacur
        if (optim$constrain) {
            thetacur0[constrained] <- exp(thetacur[constrained])
        }
        probcur <- p[seq(length(thetacur)+1,length(p))]
        pp <- modelPar(mg,thetacur0)$p
        logff <- sapply(1:k, function(j) (logLik(mg$lvm[[j]],p=pp[[j]],data=data,indiv=TRUE,model=MODEL)))
        logplogff <- t(apply(logff,1, function(z) z+log(probcur)))
        ## Log-sum-exp (see e.g. NR)
        zmax <- apply(logplogff,1,max)
        logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax    
        gamma <- exp(apply(logplogff,2,function(y) y - logsumpff)) ## Posterior class probabilities
        ## M-step:
        ##cat("thetacur=",thetacur,"\n")
        ##cat("probcur=",probcur,"\n")
        probcur <- colMeans(gamma)
        count2 <- 0
        for (jj in 1:optim$newton) {
            count2 <- count2+1
            oldpar <- thetacur
            thetacur <- Scoring(thetacur,gamma)
            if (frobnorm(oldpar-thetacur)<optim$delta) break;
        }
        p <- c(thetacur,probcur)
        if (all) {
            res <- list(p=p,gamma=gamma,
                        theta=rbind(thetacur0),
                        prob=rbind(probcur))
            return(res)
        }
        return(p)
    }
    opt <- squarem(p,fixptfn=EMstep,control=sqem.control)
    val <- EMstep(opt$par,all=TRUE)
    val <- c(val, list(member=apply(val$gamma,1,which.max),
                       k=k,
                       data=data,
                       parpos=parpos,
                       multigroup=mg,
                       model=mg$lvm,
                       logLik=NA))
    class(val) <- "lvm.mixture"
    val$vcov <- Inverse(information.lvm.mixture(val))
    return(val)
}
 
###}}} mixture

##' @export
model.frame.lvm.mixture <- function(formula,...) {
    return(formula$data)
}

##' @export
iid.lvm.mixture <- function(x,...) {
    bread <- vcov(x)
    structure(t(bread%*%t(score(x,indiv=TRUE))),bread=bread)
}

##' @export
manifest.lvm.mixture <- function(x,...) {
    manifest(x$multigroup,...)
}

##' @export
predict.lvm.mixture <- function(object,x=NULL,p=coef(object,full=TRUE),...) {
    p0 <- coef(object)
    pp <- p[seq_along(p0)]
    pr <- p[length(p0)+seq(length(p)-length(p0))];
    if (length(pr)<object$k) pr <- c(pr,1-sum(pr))
    myp <- modelPar(object$multigroup,p=pp)$p
    
    logff <- sapply(seq(object$k), function(j) (logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE)))
    logplogff <- t(apply(logff,1, function(y) y+log(pr)))
    zmax <- apply(logplogff,1,max)
    logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
    aji <- apply(logplogff,2,function(x) exp(x-logsumpff))
    gamma <- exp(apply(logplogff,2,function(y) y - logsumpff)) ## Posterior class probabilities, conditional mean
    Vgamma <- gamma-gamma^2 ## conditional variance
    M <- 0; V <- 0
    for (i in seq(object$k)) {
        m <- Model(object$multigroup)[[i]]
        P <- predict(m,data=object$data,p=myp[[i]],x=x)
        M <- M+gamma[,i]*P
        V <- V+gamma[,i]^2*attributes(P)$cond.var
    }
    structure(M,cond.var=V)
}

###{{{ logLik

ll  <- function(object,p=coef(object),prob,model="normal") {
  myp <- modelPar(object$multigroup,p)$p
  ff <- sapply(1:object$k, function(j) exp(logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE,model=model)))
  if (missing(prob))
    prob <- coef(object,prob=TRUE)
  ##  gamma <- tail(object$gamma,1)[[1]]
  ##  sum(log(colSums((prob*t(ff*gamma)))))
  ##  sum(log(colSums((prob*t(ff)))))
  loglik <- sum(log(colSums((prob*t(ff)))))
  npar <- length(prob)-1 + length(p)
  nobs <- nrow(object$data)
  attr(loglik, "nall") <- nobs
  attr(loglik, "nobs") <- nobs-npar
  attr(loglik, "df") <- npar
  class(loglik) <- "logLik"  
  return(loglik)
}

##' @export
score.lvm.mixture <- function(x,theta=c(p,prob),p=coef(x),prob,indiv=FALSE,model="normal",...) {
  myp <- modelPar(x$multigroup,p)$p
  if (missing(prob))
    prob <- coef(x,prob=TRUE)
  if (length(prob)<x$k)
    prob <- c(prob,1-sum(prob))
  logff <- sapply(seq(x$k), function(j) (logLik(x$multigroup$lvm[[j]],p=myp[[j]],data=x$data,indiv=TRUE)))
  logplogff <- t(apply(logff,1, function(y) y+log(prob)))
  zmax <- apply(logplogff,1,max)
  logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
  aji <- apply(logplogff,2,function(x) exp(x-logsumpff))
  
  scoref <- lapply(score(x$multigroup,p=p,indiv=TRUE,model=model),
                   function(x) { x[which(is.na(x))] <- 0; x })

  Stheta <- matrix(0,ncol=ncol(scoref[[1]]),nrow=nrow(scoref[[1]]))
  Spi <- matrix(0,ncol=x$k-1,nrow=nrow(Stheta))
  for (j in 1:x$k) {
    Stheta <- Stheta + apply(scoref[[j]],2,function(x) x*aji[,j])
    if (j<x$k)
      Spi[,j] <- aji[,j]/prob[j] - aji[,x$k]/prob[x$k]
  }
  S <- cbind(Stheta,Spi)
  if (!indiv)
    return(colSums(S))
  return(S)
}

##' @export
information.lvm.mixture <- function(x,...) {
  S <- score.lvm.mixture(x,indiv=TRUE,...)
  res <- t(S)%*%S
  attributes(res)$grad <- colSums(S)
  return(res)
}

##' @export
logLik.lvm.mixture <- function(object,theta=c(p,prob),p=coef(object),prob,model="normal",...) {
  myp <- modelPar(object$multigroup,p)$p
  if (missing(prob))
    prob <- coef(object,prob=TRUE)
  if (length(prob)<object$k)
      prob <- c(prob,1-sum(prob))
  logff <- sapply(1:object$k, function(j) (logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE,model=model)))
  logplogff <- t(apply(logff,1, function(y) y+log(prob)))
  ## Log-sum-exp (see e.g. NR)
  zmax <- apply(logplogff,1,max)
  logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax    
  loglik <- sum(logsumpff)
  npar <- length(prob)-1 + length(p)
  nobs <- nrow(object$data)
  attr(loglik, "nall") <- nobs
  attr(loglik, "nobs") <- nobs-npar
  attr(loglik, "df") <- npar
  class(loglik) <- "logLik"  
  return(loglik)
}

###}}} logLik

###{{{ vcov

##' @export
vcov.lvm.mixture <- function(object,...) {
  return(object$vcov)
}

###}}}

###{{{ summary/print

##' @export
summary.lvm.mixture <- function(object,labels=0,...) {
  mm <- object$multigroup$lvm
  p <- coef(object,list=TRUE)
  p0 <- coef(object)
  myp <- modelPar(object$multigroup,1:length(p0))$p
  coefs <- list()
  ncluster <- c()
  for (i in 1:length(mm)) {
    ## cc <- coef(mm[[i]],p=p[[i]],vcov=vcov(object)[myp[[i]],myp[[i]]],data=NULL)
    ## nn <- coef(mm[[i]],mean=TRUE,labels=labels,symbol="<-")
    ## nm <- index(mm[[i]])$npar.mean
    ## if (nm>0) {
    ##   nn <- c(nn[-(1:nm)],nn[1:nm])
    ## }
    ## rownames(cc) <- nn
    ## attributes(cc)[c("type","var","from","latent")] <- NULL
    cc <- CoefMat(mm[[i]],p=p[[i]],vcov=vcov(object)[myp[[i]],myp[[i]]],data=NULL,labels=labels)
    coefs <- c(coefs, list(cc))
    ncluster <- c(ncluster,sum(object$member==i))
  }
  res <- list(coef=coefs,ncluster=ncluster,prob=tail(object$prob,1),
              AIC=AIC(object),s2=sum(score(object)^2))
  class(res) <- "summary.lvm.mixture"
  return(res)
}

##' @export
print.summary.lvm.mixture <- function(x,...) {
  space <- paste(rep(" ",12),collapse="")
  for (i in 1:length(x$coef)) {
    cat("Cluster ",i," (n=",x$ncluster[i],", Prior=", formatC(x$prob[i]),"):\n",sep="")
    cat(rep("-",50),"\n",sep="")
    print(x$coef[[i]], quote=FALSE)
    if (i<length(x$coef)) cat("\n")
  }
  cat(rep("-",50),"\n",sep="")
  cat("AIC=",x$AIC,"\n")
  cat("||score||^2=",x$s2,"\n")
  invisible(par)  
}

##' @export
print.lvm.mixture <- function(x,...) {
  space <- paste(rep(" ",12),collapse="")
  for (i in 1:x$k) {
    cat("Cluster ",i," (n=",sum(x$member==i),"):\n",sep="")
    cat(rep("-",50),"\n",sep="")
    print(coef(x)[x$parpos[[i]]], quote=FALSE)
    cat("\n")   
  }
  invisible(par)
}

###}}}

###{{{ plot

##' @export
plot.lvm.mixture <- function(x,type="l",...) {
  matplot(x$theta,type=type,...)
}

###}}} plot

###{{{ coef

##' @export
coef.lvm.mixture <- function(object,iter,list=FALSE,full=FALSE,prob=FALSE,class=FALSE,...) {
  N <- nrow(object$theta)
  res <- object$theta
  if (class) return(object$gammas)
  if (list) {
      res <- list()
      for (i in 1:object$k) 
          res <- c(res, list(coef(object)[object$parpos[[i]]]))
      return(res)
  }
  if (full) {
      res <- cbind(res,object$prob[,seq(ncol(object$prob)-1)])
  }   
  if (prob) {
      res <- object$prob
  }
  if (missing(iter))
      return(res[N,])
  else
      return(res[iter,])
}

###}}} coef
