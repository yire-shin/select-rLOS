# Internal implementation for the msrlos package.
# Source: rgev2.fit.park_22Aug25.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

#--------------------------------------------------------------
#-------------------------------------------------------    
rgev2.fit.park = function(xdat=NULL, r=NULL, k.fix=NULL, num_inits=20, 
                         reltol=1e-6, show=F){
  
  # rgev.fit for k.fix cases; 
  # here k.fix is only a scalar value, not a vector
  
  rfit=list()  
  if(is.null(k.fix)) stop("k.fix should not be null")
  dtr= as.matrix(xdat)
  ntry=num_inits
  if(is.null(r)) r=dim(xdat)[2]
  numr=r
  
  if(numr==1){
      rfit[[3]]= gev.fit.kfix(xdat=as.vector(dtr[,1]), 
                              k.fix=k.fix, ntry=num_inits, varcom=FALSE) 
    # gev.fit(dtr[,1])
    
  }else{
    rfit[[3]]= rgevmle.park.kfix(xdat=dtr, numr=numr, 
                                 k.fix=k.fix, ntry=num_inits, 
                                 reltol=reltol, const=FALSE,
                                 qpro=c(.99))
    rfit[[3]]$mle = rfit[[3]]$rmle.theta
  }
  invisible(rfit)
  return(rfit[[3]])
}
#------------------------------------------------------
# gev estimation under xi fixed --single program
#------------------------------------------------------------
gev.fit.kfix =function (xdat, k.fix= NULL, ntry=10, varcom=FALSE)   
{
  # rgev.fit for k.fix cases; 
  # here k.fix is only a scalar value, not a vector

  zx <- list();  kx =list(); value=list()
  if(is.null(k.fix)) stop("k.fix should not be null")

  xifix=k.fix
  
  init.xi= matrix(0, nrow=ntry, ncol=2)
  
  init.xi <- ginit.xifix(xdat, ntry)
  xi=xifix
  
  # upsig= EnvStats::iqr(xdat)*5
  # upmu=  abs(median(xdat))*3
  
  #  par.start=rep(1,2)
  
  tryCatch(
    for(itry in 1:nrow(init.xi)){
      
      #      par.start[1:2]=init.xifix[itry,1:2]
      
      value <- try(optim(init.xi[itry,], fn=gev.xilik, 
                         #lower =c(-upmu,0), upper =c(upmu,upsig), 
                         method = "Nelder-Mead",   #method="L-BFGS-B",
                         xifix=xifix, xdat=xdat) )
      
      # control=list(trace=1, outer.iter=40, 
      #              inner.iter=200, tol=1.e-5,
      #              delta=1.e-6) ))
      
      if(is(value)[1]=="try-error"){
        kx[[itry]] <- list(value=10^6)
      }else{
        kx[[itry]] <- value
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:ntry,value=sapply(kx, function(x) x$value[which.min(x$value)]))
  
  optim_table1 <-optim_value[order(optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <- kx[[selc_num]]
  
  zx$conv <- x$convergence
  zx$nllh <- x$value[which.min(x$value)]
  
  if(zx$conv != 0){
    zx$mle = NA
  }else{
    zx$mle <- c(x$par, xifix)    # hosking style
    # zx$mle[3]= xifix
    #     zx$grad <- numDeriv::grad(gev.xilik, x$par, xifix=xifix, xdat=xdat)
    
    # if(varcom==T){
    #   Hess = PrescottW(x$par, xifix=xifix, nsam=length(xdat))
    #   #numDeriv::hessian(gev.xilik, x$par, xifix=xifix, xdat=xdat)
    #   if(is.na(Hess[1,1])) {
    #     zx$cov =matrix(NA,2,2)
    #   }else if(det(Hess) <= 0) {
    #     zx$cov =matrix(NA,2,2)
    #   }else{
    #     zx$cov <- solve(Hess)
    #   }
    # } #end if varcom
    
  }
  class(zx) <- "gev.xifix"
  return(zx)
  invisible(zx)
}
#---------- end of gev.xifix program -----------------
ginit.xifix =function(data,ntry){
  
  initx <-matrix(0,nrow=ntry,ncol=2)
  
  lmom_init = lmomco::lmoms(data,nmom=5)
  lmom_est <- lmomco::pargum(lmom_init, checklmom=T )
  
  initx[1,1]    <- lmom_est$para[1]
  initx[1,2]    <- lmom_est$para[2]
  
  maxm1=ntry; maxm2=maxm1-1
  initx[2:maxm1,1] <- initx[1,1]+rnorm(n=maxm2,mean=0,sd = 7)
  initx[2:maxm1,2] <- initx[1,2]+rnorm(n=maxm2,mean=2,sd = 2)
  initx[2:maxm1,2] = max(0.1, initx[2:maxm1,2])
  
  return(initx)
}
#----------------------------------------------------------
gev.xilik <- function(a, xifix=NULL, xdat=NULL) {
  
  mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
  sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
  # xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
  
  xi=xifix
  if(sc <= 0) return(10^6)
  
  nsam= length(xdat)
  y= rep(0,nsam)
  y <- (xdat - mu)/sc
  y <- 1 - xi * y      # park modify to negative, for xi in hosking
  
  for (i in 1:nsam){
    y[i] = max(0, y[i], na.rm=T) }
  
  if (any(y <= 0)) return(10^6)
  if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi
  }  else  {ooxi=sign(xi)*10^5}
  
  b2= sum(log(y)) * (1-ooxi)
  b3= sum(exp(ooxi * log(y)) )
  
  nllh =  nsam*(log(sc)) + b2 + b3
  
  zz=nllh 
  return(zz)
}
#------------------------------------------------------------
#--------------------------------------------------------------
#  The following function is just a simple modification of 
#  'rlarg.fit' function in the 'ismev' package.
#--------------------------------------------------------------

rlarg.fit.kfix = function (xdat, r = dim(xdat)[2], 
                          k.fix=NULL, init=NULL, ydat = NULL, mul = NULL, 
                          sigl = NULL, 
                                  shl = NULL, mulink = identity, siglink = identity, shlink = identity, 
                                  muinit = NULL, siginit = NULL, 
                          shinit = NULL, show = TRUE, 
                                  method = "Nelder-Mead", maxit = 1000, 
                          lowb= -1.0, const=FALSE, 
                                  reltol= 1e-6, ...) 
{
  # rgev.fit for k.fix cases; 
  # here k.fix is only a scalar value, not a vector

  if(is.null(k.fix)) stop("k.fix should not be null")

  z <- list()                             # hosking style para
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
  # if (is.null(shl)) {
  #   shmat <- as.matrix(rep(1, dim(xdat)[1]))
  #   if (is.null(shinit)) 
  #     shinit <- 0.1
  # }
  # else {
  #   z$trans <- TRUE
  #   shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
  #   if (is.null(shinit)) 
  #     shinit <- c(0.1, rep(0, length(shl)))
  # }
  xdatu <- xdat[, 1:r, drop = FALSE]
  
  if(is.null(init)) init <- c(muinit, siginit)  #, shinit)    # park modified these 2 lines
  init=init
  
#  z$model <- list(mul, sigl) #, shl)
#  z$link <- deparse(substitute(c(mulink, siglink))  #, shlink)))
  u <- apply(xdatu, 1, min, na.rm = TRUE)

#-----------------------------  
  rlarg.lik.kfix <- function(a, k.fix=k.fix) {
    
    mu <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
    sc <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    xi <- -k.fix #shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
    
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
#--------------------------------------
#  if(is.null(const) | const==FALSE){
    
    x <- optim(init, rlarg.lik.kfix, hessian = FALSE, method = method, 
               control = list(maxit = 1000, reltol=reltol),
               k.fix=k.fix)
  
  # }else{
  #   # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  #   
  #   x <- optim(init, rlarg.lik.kfix, hessian = FALSE, method = c("L-BFGS-B"), 
  #              lower= c(-Inf, 0), upper=c(Inf, Inf),     # hosking style para
  #              control = list(maxit = 1000),
  #              k.fix=k.fix)
  # }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
  mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
  sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
#  xi <- shlink(drop(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)])))
  
  z$conv <- x$convergence
  z$nllh <- x$value
  #  z$data <- xdat
  # if (z$trans) {
  #   for (i in 1:r) z$data[, i] <- -log((1 + (as.vector(xi) * 
  #       (xdat[, i] - as.vector(mu)))/as.vector(sc))^(-1/as.vector(xi)))
  # }
  
#  x$par[3]= k.fix 
  z$mle <- c(x$par, k.fix)      # hosking style para
  
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
#-------------------------------------------------------------  
rgevmle.park.kfix = function(xdat, numr=NULL, 
                             k.fix=k.fix, ntry=20, lowb= -1.0, 
                             reltol=1e-6, const=FALSE, qpro=NULL){
  
  # rgev.fit for k.fix cases; 
  # here k.fix is only a scalar value, not a vector

  if(is.null(k.fix)) stop("k.fix should not be null")
  
  zz=list(); k=list(); z=list()
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat[,1],ntry)
  nllh= rep(NA, ntry)
  
  tryCatch( 
    for(i in 1:nrow(init)){
      
      value= try( rlarg.fit.kfix(xdat[,1:numr],r=numr, init=init[i,1:2], 
                            k.fix=k.fix, reltol=reltol, show=F, 
                            lowb=lowb, const=const) 
                  , silent=T)              # hosking style para
      
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
  z$mle <- x$mle           # hosking style parameter
  
  zz$conv = z$conv
  zz$nllh = z$nllh
  zz$rmle.rl = lmomco::quagev(qpro, lmomco::vec2par(z$mle, type='gev'))
  zz$rmle.theta = z$mle                  # hosking style parameter
  
  return(zz)
}
