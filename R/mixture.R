
###{{{ mixture
mixture <- function(x, data, k=length(x), control, FUN, type=c("standard","CEM","SEM"),...) {    

  optim <- list(start=NULL,
                startbounds=c(-2,2),
                startmean=FALSE,
                nstart=3,
                prob=NULL,
                iter.max=500,
                delta=1e-2,
                stabil=TRUE,
                gamma=1,
                gamma2=1,
                newton=20,
                tol=1e-6,
                method="scoring",
                constrain=TRUE,
                lbound=1e-9,
                trace=1,
                lambda=0, # Stabilizing factor (avoid singularities of I)
                tol=1e-12)
  
  type <- tolower(type[1])
  if (!missing(control))
    optim[names(control)] <- control
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
  gammas <- list()
  curloglik <- logLik(mymodel,p=thetacur,prob=probcur)
  vals <- c(curloglik)
  i <- count <- 0
  member <- rep(1,nrow(data))
  E <- Inf

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
    ff <- sapply(1:k, function(j) logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE))
    return(-sum(gamma*ff))
    ## Previous:
    ##      ff <- sapply(1:k, function(j) exp(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE)))
    ##      return(-sum(log(rowSums(gamma*ff))))
  }
  myGrad <- function(p) {
    if (optim$constrain) {
      p[constrained] <- exp(p[constrained])
    }
    myp <- modelPar(mg,p)$p
    D <- lapply(1:k, function(j) gamma[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE))
    D0 <- matrix(0,nrow(data),length(p))
    for (j in 1:k) D0[,parpos[[j]]] <- D0[,parpos[[j]]]+D[[j]]
    S <- -colSums(D0)
    if (optim$constrain) {
      S[constrained] <- S[constrained]*p[constrained]
    }
    return(S)
    ## Previous:
    ##      ff <- sapply(1:k, function(j) exp(logLik(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE)))
    ##      gammaff <- gamma*ff
    ##      f0 <- rowSums(gammaff)
    ##      D <- lapply(1:k, function(j) 1/f0*(gammaff)[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data))
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
    I <- lapply(1:k, function(j) probcur[j]*information(mg$lvm[[j]],p=myp[[j]],n=nrow(data),data=data))
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
  Scoring <- function(p) {
    p.orig <- p
    if (optim$constrain) {
      p[constrained] <- exp(p[constrained])
    }
    myp <- modelPar(mg,p)$p
    D <- lapply(1:k, function(j) gamma[,j]*score(mg$lvm[[j]],p=myp[[j]],data=data,indiv=TRUE))
    D0 <- matrix(0,nrow(data),length(p))
    for (j in 1:k) D0[,parpos[[j]]] <- D0[,parpos[[j]]]+D[[j]]
    S <- colSums(D0)
    if (optim$constrain) {
      S[constrained] <- S[constrained]*p[constrained]
    }
    I <- lapply(1:k, function(j) probcur[j]*information(mg$lvm[[j]],p=myp[[j]],n=nrow(data),data=data))
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
    p.orig + optim$gamma*solve(I0)%*%S
##   p.orig + solve(I0+optim$lambda*diag(nrow(I0)))%*%S
  }
  

  on.exit(
          {
            member <- apply(gamma,1,which.max)
            res <- list(prob=probs,theta=thetas, objective=vals, gamma=gammas, k=k, member=member, data=data, parpos=parpos, multigroup=mg, model=mg$lvm, logLik=vals);
            class(res) <- "lvm.mixture"
            return(res)
          }
          )


  mytheta <- thetacur
  if (optim$constrain) {
    mytheta <- thetacur
    mytheta[constrained] <- exp(mytheta[constrained])
  }
  
  while ((i<optim$iter.max) && (E>optim$tol)) {
    if (!missing(FUN)) {
##      if (!missing(FUN) & i>0) {
      member <- apply(gamma,1,which.max)
      res <- list(prob=probs,theta=thetas, objective=vals, gamma=gammas, k=k, member=member, data=data, parpos=parpos, multigroup=mg, model=mg$lvm);      
      class(res) <- "lvm.mixture"
      dummy <- FUN(res)
    }
    i <- i+1
    count <- count+1
    if (count==optim$trace) {
      cat("Err: ", E, "\n")
##      print(E>optim$tol)      
      print(as.vector(mytheta)); count <- 0
    }

    probs <- rbind(probs,probcur)
    pp <- modelPar(mg,mytheta)$p
    logff <- sapply(1:k, function(j) (logLik(mg$lvm[[j]],p=pp[[j]],data=data,indiv=TRUE)))
##    print(probcur)
##    print(dim(logff))
##    print(logff)
##    pff <- t(apply(exp(logff),1, function(y) y*probcur))

    
    logplogff <- t(apply(logff,1, function(z) z+log(probcur)))
    ## Log-sum-exp (see e.g. NR)
    zmax <- apply(logplogff,1,max)
    logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax    

    oldloglik <- curloglik
    curloglik <- sum(logsumpff)
    vals <- c(vals,curloglik)
    
##    sumpff <- rowSums(exp(logplogff));
##    sumpff[sumpff==0] <- 1e-6
    
    gamma <- exp(apply(logplogff,2,function(y) y - logsumpff)) ## Posterior class probabilities
##    print(gamma)
    
##    gamma <- apply(logplogff,2,function(y) y - log(sumpff)) ## Posterior class probabilities
    
##    gamma <- apply(pff,2,function(y) y/sumpff) ## Posterior class probabilities
    ##  opt <- nlminb(pcur,myObj,grad=myGrad,control=list(trace=1))
    mythetaold <- mytheta
##    opt <- nlminb(thetacur,myObj,grad=myGrad, lower=lower, control=list(...))
    ##    print(thetacur)    

    
    #### M-step    
    if (type%in%c("sem","cem")) {
      if (type=="sem")
        idx <- apply(gamma,1,function(rr) which(rmultinom(1,1,rr)==1))
      else
        idx <- apply(gamma,1,function(rr) which.max(rr))

      mydata <- list()
      for (ii in 1:k) {
        mydata <- c(mydata, list(data[idx==ii,,drop=FALSE]))
      }      
      mymg <- multigroup(x,mydata,fix=FALSE)
      e <- estimate(mymg,silent=TRUE,control=list(start=mytheta))
      mytheta <- pars(e)
      probcur <- sapply(1:k,function(j) sum(idx==j)/length(idx))
      
    } else {    
      if (optim$method=="BFGS") {
        opt <- nlminb(thetacur,myObj,gradient=myGrad,hessian=myInformation,lower=lower, control=list(iter.max=10))  
        thetacur <- opt$par
      } else {     
        probcur <- colMeans(gamma)
        oldpar <- thetacur
        count2 <- 0
        for (jj in 1:optim$newton) {
          count2 <- count2+1
          ##        aa <- list(thetacur=thetacur, gamma=gamma, mg=mg, optim=optim, constrained=constrained, parpos=parpos, probcur=probcur)
          ##        save(aa, file="aa.rda")
          ##        thetacur <- thetacur + solve(myInformation(thetacur))%*%(-myGrad(thetacur))
          oldpar <- thetacur
          ##        print(thetacur)
          thetacur <- Scoring(thetacur)
          if (frobnorm(oldpar-thetacur)<optim$delta) break;
        }
        cat(count2, " Scoring iterations.\n")
      }    
      if (optim$constrain) {
        mytheta <- thetacur
        mytheta[constrained] <- exp(mytheta[constrained])
      }
    }
      thetas <- rbind(thetas,as.vector(mytheta))
      gammas <- c(gammas, list(gamma))
      E <- sum((mytheta-mythetaold)^2)
  }
  

  member <- apply(gamma,1,which.max)
  res <- list(prob=probs,theta=thetas, objective=vals, gamma=gammas, k=k, member=member, data=data, parpos=parpos, multigroup=mg, model=mg$lvm, logLik=vals);
  class(res) <- "lvm.mixture"
  return(res)
}
###}}} mixture

###{{{ logLik

ll  <- function(object,p=coef(object),prob) {  
  myp <- modelPar(object$multigroup,p)$p
  ff <- sapply(1:object$k, function(j) exp(logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE)))
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

logLik.lvm.mixture <- function(object,p=coef(object),prob) {
  myp <- modelPar(object$multigroup,p)$p
  if (missing(prob))
    prob <- coef(object,prob=TRUE)
  if (length(prob)<object$k)
    prob <- c(prob,1-sum(prob))
  logff <- sapply(1:object$k, function(j) (logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE)))
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

vcov.lvm.mixture <- function(object,...) {
  return(NULL)
}

###}}}

###{{{ summary/print

summary.lvm.mixture <- function(object,...) {
  object
}

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

plot.lvm.mixture <- function(x,...) {
  matplot(x$theta,...)
}

###}}} plot

###{{{ coef

coef.lvm.mixture <- function(object,iter,list=FALSE,prob=FALSE,class=FALSE,...) {
  N <- nrow(object$theta)
  res <- object$theta
  if (list) {
    res <- list()
    for (i in 1:object$k) 
      res <- c(res, list(coef(object)[object$parpos[[i]]]))
    return(res)
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
