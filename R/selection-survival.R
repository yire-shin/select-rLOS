# Internal implementation for the msrlos package.
# Source: sel_surv_batch_11Aug26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.


#----------Algoritm survival ----------------------------------
sel.surv = function(xdat, maxr=NULL, h.fix=NULL,
                         crit="AIC", num_inits=10,
                         show=TRUE, sig.ed=0.05,
                         qq=c(.98,.99,.995,.998),
                         true.para= NULL, dmin=0.03,
                         survtest=NULL){
  
  xdat=na.omit(xdat)
  z=list(); 
  nsam= dim(xdat)[1]
  R= maxr

  if(is.null(maxr)) R= maxr= dim(xdat)[2]
  if(maxr < 1) stop("maxr should be greater than or equal to 1")
  
  if(is.null(h.fix)){
    h.sel= h.select(xdat, sig.ed=sig.ed, dmin=dmin)
    h.fix= h.sel$h.fix
    z$optr_rk4d = h.sel$optr
  }
#  z$h.fix=h.fix
  
  nhfix=length(h.fix)
  nmod= 6+length(h.fix)
  
  k3d.name=NULL
  if(!is.null(h.fix)) k3d.name= rep('rk3d',nhfix)
  k3d.name2=NULL
  if(!is.null(h.fix)) {
    k3d.name2= paste(paste('rk3d','(h=',sep="_"),
                     round(h.fix,3),')',sep="") }

  name.m = c('rglo','rggd', 'rgev', 'rk4d', 'rld','rgd', 
             k3d.name)
  name.m2 = c(name.m[1:6], k3d.name2)
  
  z$model.names = name.m2
  
  npar=c(3,3,3,4,2,2,rep(3,nhfix))
#  work.fit=paste(name.m,'.fit.park',sep="")
  
  wid=rep(NA,maxr); worst= rep(NA,maxr); 
#  work.order =min(maxr, nmod-1)
  caic=matrix(NA, maxr,nmod)

  numr=1
  while(numr <= maxr){

     if(numr==1) surv= c(1:nmod)       
    # rmle=list()
     aic=rep(NA,nmod)
    
    for (mid in surv){
     if(numr==1){
      aic[mid]= as.numeric(survtest$redtest[[mid]][R,11])
     }else if(numr >=2) {
      aic[mid]= as.numeric(survtest$redtest[[mid]][numr-1,11])
     }
    } # end for
     
     wid[numr] = which.max(aic)
     if(max(aic, na.rm=TRUE)-min(aic, na.rm=TRUE) > 4) {
       worst[numr] = name.m2[wid[numr]] }
    
    recov= which.min(aic)
    surv= unique(sort(c(surv,recov)))

    if(!is.na(wid[numr])){
      surv = surv[-which(surv== wid[numr])]}
    
    if(show==TRUE){
      cat("numr=", numr," elim=",wid[numr]," best=",recov,
          " / surv= ", surv,"\n")
    }

     caic[numr,]=aic
     
     if(length(surv)==1){
       bestm= surv.model =name.m2[surv]
       bid=surv
       break
     }

    if(numr < maxr) {numr=numr+1
    }else if(numr==maxr){break}

  } #end while numr in 1:maxr ------------------------
  
  surv.model = name.m2[surv]

  if(length(surv) >= 1){
    bid= which.min(aic)
    bestm= name.m2[bid]
  }
    delta= aic-min(aic, na.rm=TRUE)     # smooth AIC
    bid.good= which(delta <= 2)         # choose good models
    bestm.good= name.m2[bid.good] 
  
#  z$surv.model= surv.model
  z$algorithm= "surv"
  z$mstar =bid
  # z$aic= t(caic[1:(numr),])
  # z$elim.model=worst[1:(numr-1)]; 

  # find the optimal r for the best model

  opt_r=rep(NA,nmod); theta=matrix(NA, nmod, 4)
  opt_r= survtest$opt_r
    
    for(mid in bid.good){
      if(opt_r[mid]==1) {rcol=R
      }else if(opt_r[mid] >= 2) {rcol=opt_r[mid]-1}
      mat = as.matrix(survtest$redtest[[mid]])
      theta[mid,1:4]= mat[rcol,c(4:6,10)]
    }
  
#  z$opt_r = opt_r[bid.good]
  z$rstar= opt_r[bid]
  z$best.model= paste(bestm,opt_r[bid],sep=" ")
  z$theta = theta[bid,]
  z
}
#----------------------------------------------------------------  
h.select = function(xdat, numh=8,
                    sig.ed= 0.05, dmin=0.03, h.min= -1.5,
                    low.xi= -0.7, uph= 2.20){
  
  maxr= ncol(xdat)  # park
  
  somek4d= multi.rEdtest.park(xdat, model="rk4d", 
                              method="ed",par = NULL)
  
  opt_r4= max(one_optr(somek4d, mid=4, sigL=sig.ed, 
                       h.fix=NULL),2)
  
  rdata= xdat[,1:opt_r4]
  rmle = c(somek4d$est.loc[opt_r4-1],
           somek4d$est.scale[opt_r4-1],
           somek4d$est.shape[opt_r4-1],
           somek4d$est.shape2[opt_r4-1])
  
  # low.xi= -0.7
  # uph= 2.20
  klmom= lmomco::lmoms(xdat[, 1])
  kappar <- lmomco::parkap(klmom, snap.tau4=TRUE, 
                           nudge.tau4=1.e-3)$para
  if(is.na(kappar[1])){
    kappar <- lmomco::pargev(klmom)$para
    kappar[4]= -0.01
  }
  kappa2=kappar[2]
  ropt= opt_r4
  
  hess= optimHess(rmle, rk4d.lik.park.fisher,
                  uph=uph, low.xi=low.xi,
                  rdata=rdata, kappa2=kappa2, ropt=ropt)
  
  determ = det(hess)
  if(determ > 1.e-5) {
    cov <- solve(hess)
    se <- sqrt(diag(cov))
  }else{
    se=rep(0.05,4)
  }
  
  # set h values
  L <- qnorm(0.005)   # -2.575829
  U <- qnorm(0.995)   #  2.575829
  pp <- pnorm(L) + (0:(numh-1)) * (pnorm(U) - pnorm(L))/(numh-1)
  c_j <- qnorm(pp)

  h.fix= rmle[4]+ c_j*max(se[4],0.05)

  # filter out h values close to h hat at opt_r or h=-1 or 0
  if(maxr >= 2){
   did = which(h.fix >= 1/(maxr-1))
   if(length(did) >0) h.fix= h.fix[ -did ]
  }
  
  if(ropt >= 2){  #---------------------- added
    did = which(h.fix >= 1/(ropt-1))
    if(length(did) >0) h.fix= h.fix[ -did ]
  }  # ----------------------------------
  
  did2= integer(0)
  if(rmle[3] < 0){
    did2= which(h.fix >= -1/rmle[3] & h.fix < 0)
  }else if(rmle[3] > 0){
    did2= which(h.fix <= -1/rmle[3] & h.fix < 0)
  }
  if(length(did2) > 0) h.fix= h.fix[ -did2 ]
  
  disth = abs(h.fix-somek4d$est.shape2[opt_r4-1] )
  if(length(which(disth < dmin)) > 0){
    h.fix = h.fix[-which(disth < dmin)]}
  disth = abs(h.fix-(-1.0) )
  if(length(which(disth < dmin)) > 0){
    h.fix = h.fix[-which(disth < dmin)]}
  disth = abs(h.fix-(0.0) )
  if(length(which(disth < dmin)) > 0){
    h.fix = h.fix[-which(disth < dmin)]}
  
  if( is.null(h.fix) | all(is.na(h.fix)) ) {
     h.fix= c(-0.4, -0.3, -0.2, -0.1)
     
  }else if( !is.null(h.fix) & !all(is.na(h.fix)) ){
     if(min(h.fix) >= 0) h.fix= unique(c(-0.3,-0.2,-0.1, h.fix))
  }
  
  did2= integer(0)
  if(rmle[3] < 0){
    did2= which(h.fix >= -1/rmle[3] & h.fix < 0)
  }else if(rmle[3] > 0){
    did2= which(h.fix <= -1/rmle[3] & h.fix < 0)
  }
  if(length(did2) > 0) h.fix= h.fix[ -did2 ]
  
  did3 = which(h.fix <= h.min)  # park
  if(length(did3) > 0) h.fix= h.fix[ -did3 ]

  z=list()
  z$se=se
  z$h.fix= h.fix
  z$optr= opt_r4
  z$rmle= rmle
  z
}
#----------------------------------------------------
#------------------------------------------------------
rk4d.lik.park.fisher <- function(a, uph=uph, low.xi=low.xi,
                                 rdata=rdata, kappa2=kappa2,
                                 ropt=ropt) {
  
  # if( z$trans == FALSE){
  mu=a[1]
  sc=a[2]
  xi=a[3]
  h=a[4]
  penalty=0
  xdat=rdata
  penh=penk=NULL
  r=ropt
  
  # }else{
  #   mu <- drop(mumat %*% (a[1:npmu]))
  #   sc <- drop(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
  #   xi <- drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
  #   h <- drop(hmat %*% (a[seq(npmu + npsc + nph + 1, length = nph)]))
  # }
  
  if( any(is.na(a)) | any(is.null(a))) return(10^6)  # park
  if( xi <=  low.xi | xi > 2.999) return(10^6) # park
  if( abs(xi) <= 0.001) xi = sign(xi)*0.001
  if( abs(h) <= 0.001) h = sign(h)*0.001
  if( sc <= 0.001) sc = 0.001
  #    if( h <  -.999 | h > .999) return(10^6) # park
  if( h <=  -4.5 | h > uph) return(10^6) # park
  if(sc > 5*kappa2) return(10^6) # park
  if( h < 0 & xi*h <= -1) return(10^6) # park
  
  ri <- (r - seq(1:(r)))
  cr <- (1 - ri * h[1])
  if (any(sc <= 0) | any(cr <= 0)) 
    return(10^6)
  y <- 1 - xi * (xdat - mu)/sc
  
  if ( min(y, na.rm=T) <= 0) return(10^6)
  zr <- as.matrix(xdat[, r], ncol = 1)
  
  f <- 1 - h * (1 - xi * (zr - mu)/sc)^(1/xi)
  
  po = 1 - xi * (zr - mu)/sc
  if( min(po, na.rm=T) <= 0) {
    pou=NA
  }else{
    pou = po^(1/xi)
  }
  
  if(is.null(pou)) pou=1e-11
  pou[is.na(pou)]= 1e-11
  
  if( any(abs(pou) <= 1e-10) ) {
    id= which(abs(pou) <= 1e-10) 
    pou[id]= sign(pou[id])*1e-10
  }
  pool = 1/pou
  
  if (max(h) > min(pool, na.rm = TRUE) ) 
    #      if (max(h) > min(1/(1 - xi * (zr - mu)/sc)^(1/xi), na.rm = TRUE)) 
    return(10^6)
  if (r >= 2 && min(h) > (1/(r - 1))) 
    return(10^6)
  if (min(y, na.rm = TRUE) <= 0 || min(f, na.rm = TRUE) <= 
      0 || max(f^(1/h), na.rm = TRUE) > 1) 
    return(10^6)
  y <- log(sc) + (1 - 1/xi) * log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)
  
  
  if (is.null(penk) == F) {
    if (penk == "CD") {
      if (xi[1] >= 0) {
        p_k <- 1
      }
      else if (xi[1] > -1 && xi[1] < 0) {
        p_k <- exp(-((1/(1 + xi[1])) - 1))
      }
      else if (xi[1] <= -1) {
        p_k <- 0
      }
    }
    else if (penk == "MS") {
      if (xi[1] >= 0.5 || xi[1] <= -0.5) 
        return(10^6)
      p_k <- ((0.5 + xi[1])^(6 - 1)) * ((0.5 - xi[1])^(9 - 
                                                         1))/beta(6, 9)
    }
  }
  if (is.null(penh) == F) {
    if (penh == "MS") {
      if (h[1] <= -0.5 || h[1] >= 0.5) 
        return(10^6)
      Bef <- function(x) {
        ((0.5 + x)^(6 - 1)) * ((0.5 - x)^(9 - 1))
      }
      Be <- stats::integrate(Bef, lower = -0.5, upper = 0.5)[1]$value
      p_h <- ((0.5 + h[1])^(6 - 1)) * ((0.5 - h[1])^(9 - 
                                                       1))/Be
    }
    else if (penh == "MSa") {
      if (h[1] >= 1.2 || h[1] <= -1.2) 
        return(10^6)
      Bef <- function(x) {
        ((1.2 + x)^(6 - 1)) * ((1.2 - x)^(9 - 1))
      }
      Be <- stats::integrate(Bef, lower = -1.2, upper = 1.2)[1]$value
      p_h <- ((1.2 + h[1])^(6 - 1)) * ((1.2 - h[1])^(9 - 
                                                       1))/Be
    }
  }
  if (is.null(penk) == T & is.null(penh) == T) {
    penalty <- 0
  }
  else if (is.null(penk) == F & is.null(penh) == F) {
    penalty <- r * log(p_k * p_h)
  }
  else if (is.null(penk) == F) {
    penalty <- r * log(p_k)
  }
  else if (is.null(penh) == F) {
    penalty <- r * log(p_h)
  }
  
  if( any(abs(y) > 10^6) ) return(10^6)
  
  sum((r * h - 1)/h * log(f) + y, na.rm = TRUE) - penalty
}
#----------------------------------------------------  
