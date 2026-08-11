# Internal implementation for the msrlos package.
# Source: rgev.fit.park_17Oct25.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

#--------------------------------------------------------------
#-------------------------------------------------------    
rgev.fit.park = function(xdat=NULL, r=NULL, num_inits=10, 
                         reltol=1e-6, show=FALSE, start.para=NULL){
  
  rfit=list()  
  ntry=num_inits
  if(is.null(r)) r=dim(xdat)[2]
  numr=r
  dtr=xdat
  
  if(numr==1){
    rfit[[3]]= gev.max.consT(xdat=as.vector(dtr[,1]), 
                             ntry=num_inits, lowb= -0.7, 
                             reltol=reltol, const=TRUE, 
                             start.para=start.para) 
    # gev.fit(dtr[,1])
    
  }else{
    rfit[[3]]= rgevmle.park(xdat=dtr, numr=numr, 
                            ntry=num_inits, lowb= -0.7, 
                            reltol=reltol, const=TRUE, qpro=c(.99),
                            start.para=start.para)
    
    rfit[[3]]$mle = rfit[[3]]$rmle.theta
  }
  invisible(rfit)
  return(rfit[[3]]) 
}
#------------------------------------------------------
#-------------------------------------------------------------  
rgevmle.park = function(xdat, numr=NULL, ntry=10, lowb= -1.0, 
                        reltol=1e-6, const=TRUE, qpro=NULL,
                        start.para=NULL){
  
  zz=list(); k=list(); z=list()
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat[,1],ntry)
  if(!is.null(start.para)) init[ntry,]=start.para
  nllh= rep(NA, ntry)
  
  tryCatch( 
    for(i in 1:nrow(init)){
      
      value= try( rlarg.fit.consT.stnry(xdat[,1:numr],r=numr, init=init[i,1:3], 
                                        reltol=reltol, show=FALSE, 
                                        lowb=lowb, const=const) 
                  , silent=TRUE)                               # coles style para
      
      if(is(value)[1]=="try-error"){
        k[[i]] <- list(value=10^6)
        #       cat("i try error= ", i,"\n" )
      }else{
        k[[i]] <- value
        nllh[i]= k[[i]]$nllh
      }
      
    } #for  
  ) #tryCatch
  
  selc_num = which.min(nllh)
  
  x  <-k[[selc_num]]
  
  z$conv <- x$conv
  z$nllh <- x$nllh
  z$mle <- x$mle             # coles style parameter
  z$mle[3] = - x$mle[3]       # Hosking style para
  
  
  if( z$mle[3] <= lowb & const==TRUE ) {
    
    nllh= rep(NA, ntry)
    tryCatch( 
      for (i in 1:nrow(init)) {
        
        value= try( rlarg.fit.consT.stnry(xdat[,1:numr],r=numr, 
                                          init=init[i,1:3], 
                                          reltol=reltol, show=FALSE, 
                                          lowb=lowb, const=const) 
                    , silent=TRUE)                       # coles style para
        
        if(is(value)[1]=="try-error"){
          k[[i]] <- list(value=10^6)
        }else{
          k[[i]] <- value
          nllh[i]= k[[i]]$nllh
        }
        
      } #for
    ) #tryCatch
    
    selc_num = which.min(nllh)
    x  <-k[[selc_num]]
    
    z$conv <- x$conv
    z$nllh <- x$nllh
    z$mle <- x$mle              # coles style parameter
    z$mle[3] = - x$mle[3]       # Hosking style para
    
  }
  
  zz$conv = z$conv
  zz$nllh = z$nllh
  zz$rmle.rl = lmomco::quagev(qpro, lmomco::vec2par(z$mle, type='gev'))
  zz$rmle.theta = z$mle                                # hosking style parameter
  
  return(zz)
}
#--------------------------------------------------------------
#  The following function is just a simple modification of 
#  'rlarg.fit' function in the 'ismev' package.
#--------------------------------------------------------------

rlarg.fit.consT.stnry = function (xdat, r = dim(xdat)[2], init=NULL, ydat = NULL, mul = NULL, sigl = NULL, 
                                  shl = NULL, mulink = identity, siglink = identity, shlink = identity, 
                                  muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
                                  method = "Nelder-Mead", maxit = 1000, lowb=-1, const=NULL, 
                                  reltol= 1e-6, ...) 
{
  

  z <- list()                             # coles style para
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat[, 1]))/pi
  in1 <- mean(xdat[, 1]) - 0.57722 * in2
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(muinit)) 
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(siginit)) 
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  xdatu <- xdat[, 1:r, drop = FALSE]
  
  #    init <- c(muinit, siginit, shinit)    # park modified these 2 lines
  init = init
  
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  u <- apply(xdatu, 1, min, na.rm = TRUE)
  
#---------------------------------------  
  rlarg.lik <- function(a) {
    
    mu <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
    sc <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    xi <- a[3] #shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
    
    if (any(sc <= 0)) 
      return(10^6)

    y <- 1 + xi * (xdatu - mu)/sc
    if (min(y, na.rm = TRUE) <= 0) 
      l <- 10^6
    else {
      y <- (1/xi + 1) * log(y) + log(sc)
      y <- rowSums(y, na.rm = TRUE)
      l <- sum((1 + xi * (u - mu)/sc)^(-1/xi) + y)
    }
    l
  }
#-----------------------------------------------------
  
  if(is.null(const) | const==FALSE){
    
      x <- optim(init, rlarg.lik, hessian = FALSE, method = method, 
               control = list(maxit = maxit, reltol=reltol,...))
  }else{
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

     x <- optim(init, rlarg.lik, hessian = FALSE, method = c("L-BFGS-B"), 
             lower= c(-Inf, 0, -1.0), upper=c(Inf, Inf, -lowb),     # coles style para
             control = list(maxit = maxit, reltol=reltol,...))
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
  mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
  sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
  xi <- shlink(drop(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)])))
  z$conv <- x$convergence
  z$nllh <- x$value
  #  z$data <- xdat
  if (z$trans) {
    for (i in 1:r) z$data[, i] <- -log((1 + (as.vector(xi) * 
                                               (xdat[, i] - as.vector(mu)))/as.vector(sc))^(-1/as.vector(xi)))
  }
  z$mle <- x$par
  # z$cov <- solve(x$hessian)
  # z$se <- sqrt(diag(z$cov))
  #  z$vals <- cbind(mu, sc, xi)
  z$r <- r
  if (show) {
    if (z$trans) 
      print(z[c(2, 3)])
    print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "rlarg.fit"
  invisible(z)
  
  return(z)
}
# ----------------------------------------------------------------
#------------------------------------------------------------
gev.max.consT=function (xdat, ntry=10, lowb= -1.0, 
                        reltol=1e-6, const=TRUE, start.para=NULL) 
{
  z <- list();  k =list()           # hosking style para
  n=ntry
  
  nsample=length(xdat)
  z$nsample=nsample
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat,ntry)
  if(!is.null(start.para)) init[ntry,]=start.para
  
  #--------------------------------------------------------- 
  # The following function is a simple modification of 
  # 'gev.lik' function in the 'ismev' package.
  # --------------------------------------------------------
  gev.lik.max <- function(a) {
    
    mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- a[3]      #shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    y <- (xdat - mu)/sc
    y <- 1 - xi * y        # park modify to negative, for xi in hosking
    
    for (i in 1:nsample){
      y[i] = max(0, y[i], na.rm=T) }
    
    if (any(y <= 0) || any(sc <= 0)) 
      return(10^6)
    
    if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi       # park modify to xi in hosking
    }  else  {ooxi=sign(xi)*10^5}
    
    zz=nsample*(log(sc)) + sum( exp(ooxi *log(y)) ) + sum(log(y) * (1-(ooxi)) ) 
    
    return(zz)
  }
  #-------------------------------------------------------------
  tryCatch(
    for(i in 1:nrow(init)){
      
      value <- try(Rsolnp::solnp(init[i,], fun=gev.lik.max, 
                         LB =c(-Inf,0,lowb),UB =c(Inf,Inf,1.0),     # hosking style para
                         control=list(trace=0, outer.iter=40,
                                      delta=1.e-6, inner.iter=100, 
                                      tol=reltol) ))
      
      if(is(value)[1]=="try-error"){
        k[[i]] <- list(value=10^6)
      }else{
        k[[i]] <- value
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:n,value=sapply(k, function(x) x$value[which.min(x$value)]))
  
  optim_table1 <-optim_value[order(optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <-k[[selc_num]]
  
  #  mu <- x$par[1];  sc <- x$par[2];  xi <- x$par[3] 
  
  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$mle <- x$par                            # hosking style parameter
  
  return(z)
}


#--------------------------------------------------------------
ginit.max <-function(data,ntry){
  
  n=ntry
  init <-matrix(rep(0,n*3),ncol=3)
  
  lmom_init = lmomco::lmoms(data,nmom=5)
  lmom_est <- lmomco::pargev(lmom_init)
  
  init[1,1]    <-lmom_est$para[1]
  init[1,2]    <-lmom_est$para[2]
  init[1,3]    <-lmom_est$para[3]
  
  maxm1=ntry; maxm2=maxm1-1
  init[2:maxm1,1] <- init[1,1]+ rnorm(n=maxm2,mean=0,sd = 5)
  init[2:maxm1,2] <- abs( init[1,2]+ rnorm(n=maxm2,mean=5,sd = 5)) +1
  init[2:maxm1,3] <- runif(n=maxm2,min= -0.5,max=0.5)
  #    init[2:maxm1,2] = max(0.1, init[2:maxm1,2])
  
  return(init)
}

#---------------------------------------------    
rgev.lik.cvnll.park <- function(a, r=NULL, xdat=NULL) {
  mu <- a[1] # mulink(drop(mumat %*% (a[1:npmu])))
  sc <- a[2] # siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
  
  xi <- -a[3] # hosking style xi         # shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
  
  if (any(sc <= 0)) 
    return(10^6)
  
  xdatu <- xdat[, 1:r, drop = FALSE]
  u <- apply(xdatu, 1, min, na.rm = TRUE)
  
  y <- 1 + xi * (xdatu - mu)/sc            # xi: coles style
  if (min(y, na.rm = TRUE) <= 0) {
    l <- 10^6
  }else {
    y <- (1/xi + 1) * log(y) + log(sc)
    y <- rowSums(y, na.rm = TRUE)
    l <- sum((1 + xi * (u - mu)/sc)^(-1/xi) + y)
  }
  l
}
#-------------------------------------------------------   
