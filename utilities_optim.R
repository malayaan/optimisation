# set of utilities for simple optimization program
# 
# Rodolphe Le Riche, CNRS LIMOS, july 2021
# 

# calculate f and gradient of f by forward finite difference
f.gradf <- function(x,f,h=1.e-8){
  d<-length(x)
  res <- list()
  res$gradf <- rep(NA,d)
  res$fofx <- f(x)
  for (i in 1:d){
    xp <- x
    xp[i] <- x[i]+h
    res$gradf[i] <- (f(xp)-res$fofx)/h
  }
  return(res)
}

# record points online
updateRec <- function(arec,x,f,t){
  if (is.null(arec$X)) {arec$X<-x} else {arec$X <- rbind(arec$X,x)}
  if (is.null(arec$F)) {arec$F<-f} else {arec$F <- c(arec$F,f)}
  if (is.null(arec$Time)) {arec$Time<-t} else {arec$Time <- c(arec$Time,t)}
  return(arec)
}

# L2 norm
l2norm <- function(x){
  return(sqrt(sum(x^2)))
}

# plot contour of function when d==2
plot_contour <- function(LB,UB,f){
  no.grid <- 100
  x1 <- seq(LB[1], UB[1], length.out=no.grid)
  x2 <- seq(LB[2], UB[2], length.out=no.grid)
  x.grid <- expand.grid(x1, x2)
  z <- apply(x.grid, 1, f)
  z.grid <- matrix(z, no.grid)
  contour(x1, x2, z.grid, nlevels=20, xlab="x1", ylab="x2")
  # add global optimum
  glob_xstar <- f(xx = rep(0,2),get_glob_xstar=TRUE)
  points(x = glob_xstar[1],y=glob_xstar[2],pch=3,cex=1.5,col="green")
}

#### 3D plot
plot_2Dfun_in_3D <- function(LB,UB,f,do_rgl=FALSE){
  no.grid <- 100
  x1 <- seq(LB[1], UB[1], length.out=no.grid)
  x2 <- seq(LB[2], UB[2], length.out=no.grid)
  x.grid <- expand.grid(x1, x2)
  z <- apply(x.grid, 1, f)
  z.grid <- matrix(z, no.grid)
  mypersp <- persp(x = x1,y=x2,z=z.grid,zlab = "f")
  # add global optimum
  glob_xstar <- f(xx = rep(0,2),get_glob_xstar=TRUE)
  points(x = trans3d(glob_xstar[1],glob_xstar[2],f(glob_xstar),pmat=mypersp), pch=3,cex=1.5,col="green")
  ## Below is the nicer interactive 3D RGL version
  if (do_rgl) {
    library("rgl")
    open3d()
    surface3d(x1, x2, z.grid, col= "lightblue")
    points3d(glob_xstar[1], glob_xstar[2], f(glob_xstar), pch=19, col="green", size=10)
    title3d("a 2D function", col="blue", font=4)
    decorate3d()
    aspect3d(1, 1, 1)    
  }
}


# increasing sequence: useful to plot first points and then fewer and fewer
inc.geom.seq <- function(from=1,to=1000,coef=1.4)
{
  s <- c(round(from))
  x <- from
  i <- 1
  ieff <- 1
  done <- FALSE
  while (!done){
    x <- x*coef
    sp <- round(x)
    if (sp != s[ieff]){
      s <- c(s,sp)
      ieff <- ieff+1
      if (sp>to) done<-TRUE
    }
    i <- i+1
  }
  s<-s[-ieff]
  return(s)
}

# plot results of an optimization run
plot_res_optim <- function(PbFormulation,optAlgoParam,res,printlevel,nout=200)
{
  
  if (printlevel >= 4){resF <- res$rec$F}
  else {resF <- res$rBest$F}
  observedLowest <- min(resF)
  if (observedLowest >= 0){delta <- 1} 
  else {delta <- ceiling(1+abs(observedLowest))}
  ylim = c(min(log(delta+resF)),max(log(delta+resF)))
  
  if (length(res$rBest$Time)>=2){
    cleaned_rBest <- approx(x = res$rBest$Time,y=log(delta+res$rBest$F),
                            xout = seq(1,max(optAlgoParam$budget,tail(res$rBest$Time,n=1)),
                                       length.out=nout),rule = 2)
    
    plot(cleaned_rBest,type = "l",xlab = "nb. fct. eval",ylab=paste0("log(",delta,"+f)"),col="red",ylim = ylim)
    if (printlevel >= 4){
      lines(x = res$rec$Time,y = log(delta+res$rec$F),col="blue")
    }
  } else {cat("\n no convergence plot done, only 1 point to draw\n")}
  # x plots if in d=2
  if (d==2){
    plot_contour(LB=pbFormulation$LB,UB=pbFormulation$UB,f=pbFormulation$fun) 
    if (printlevel>=4){
      points(res$rec$X[,1], res$rec$X[,2], pch=20, col="blue", cex=0.5)  
      iplot <- inc.geom.seq(from=1,to=length(res$rec$Time),coef=1.4) # selection of points, otherwise unreadable
      text(x=res$rec$X[iplot,1], y=res$rec$X[iplot,2],labels=res$rec$Time[iplot], pos=3, cex=1.0)        
    } else if (printlevel==2){
      points(res$rBest$X[,1], res$rBest$X[,2], pch=20, col="red", cex=0.5)  
      iplot <- inc.geom.seq(from=1,to=length(res$rBest$Time),coef=1.4) # selection of points, otherwise unreadable
      text(x=res$rBest$X[iplot,1], y=res$rBest$X[iplot,2],labels=res$rBest$Time[iplot], pos=3, cex=1.0)        
    }
  }
}

# repeat runs for the optimiser
# doplots concerns only the repeatedRuns actions. printlevel is passed to the optimizer level
repeatedRuns <- function(optimizer,pbFormulation, algoParam,no_test=10,printlevel=1,silent=FALSE) {
  if (!silent) {cat("******* Start testing \n")}
  hist <- list()
  maxTime <- 0
  for(i in 1:no_test){
    if (!silent) {cat("   test no.",i,"\n")}
    res <- optimizer(pbFormulation=pbFormulation,algoParam=optAlgoParam,printlevel=printlevel)
    hist[[i]] <- res
    maxTime <- max(c(maxTime,optAlgoParam$budget,res$nbFun))
  }
  if (!silent) {cat("******* Done testing \n")}
  
  # Post-processing results
  ngrid <- 200
  tgrid <- seq(from=1,to=maxTime,length.out=ngrid)
  cleanedF <- matrix(data=NA,nrow = no_test, ncol = ngrid)
  allFbest <- rep(NA,no_test)
  for(i in 1:no_test){
    oneF <- approx(x=hist[[i]]$rBest$Time,y=hist[[i]]$rBest$F,xout = tgrid, rule = 2)
    cleanedF[i,]<-oneF$y
    allFbest[i]<-hist[[i]]$fbest
  }
  quantiles<-apply(X = cleanedF , MARGIN = 2 , FUN = quantile , c(0.25,0.5,0.75))
  # plot all curves
  if (!silent){
    observedLowest <- min(cleanedF)
    if (observedLowest >= 0){delta <- 1} 
    else {delta <- ceiling(1+abs(observedLowest))}
    matplot(x = tgrid,y = t(log(delta+cleanedF)),type="l",xlab="nb. fct. calls",ylab = paste0("log(",delta,"+f)"),lty=1)
    # plot mediam and quartiles
    matplot(x = tgrid,y = t(log(delta+quantiles)),type=rep("l",3),
            xlab="nb. fct. calls",ylab = paste0("log(",delta,"+quantiles(f))"),col=rep(2,3),lty = c(2,1,2))
    hist(allFbest,xlab="f",main="histogram of final fbest's")
  }
  return(list(Fquantiles=quantiles,Time=tgrid,Fbest=allFbest))
}

