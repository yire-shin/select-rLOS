# Internal implementation for the msrlos package.
# Source: redtest_all_10Aug26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

#----------------------------------------------------------  
redtest.all = function(xdat, h.fix=NULL, numh=8,
                       sig.ed=0.05, dmin=0.03){
  
  z=list()
  if(is.null(h.fix)){
    h.sel= h.select(xdat, numh=numh, sig.ed=sig.ed, dmin=dmin)
    z$h.fix = h.fix= h.sel$h.fix
    # z$optr_rk4d = h.sel$optr
  }
  nhfix=length(h.fix)
  nmod= 6+ nhfix
  R = ncol(xdat)
  
  k3d.name=NULL
  if(!is.null(h.fix)) k3d.name= rep('rk3d',nhfix)

  name.m = c('rglo','rggd', 'rgev', 'rk4d', 'rld','rgd', 
             k3d.name)

  npar=c(3,3,3,4,2,2,rep(3,nhfix))
  opt_r=rep(NA,nmod)
  
  U = theta.req1(xdat, surv=c(1:nmod), num_inits=5,
                 h.fix=h.fix)

  xtest=list()
  for (mid in 1:nmod){
    
    h.mid <- if (mid >= 7L) h.fix[mid - 6L] else NULL
    xtest[[mid]] = multi.rEdtest.park(xdat, model=name.m[mid],
                           h.fix=h.mid,
                           method="ed", num_inits=5)

    xtest[[mid]][R,1]=1
    xtest[[mid]][R,c(4:6,10)] = U[mid,]

    opt_r[mid]= one_optr(xtest[[mid]],mid=mid,
                          sigL=sig.ed, h.fix=h.fix)

    theta= rep(NA,4)
    for (r in 1:R){
      
      if(r==1){ theta = as.matrix(xtest[[mid]])[R,c(4:6,10)]
      
      }else if(r >= 2){
        theta = as.matrix(xtest[[mid]])[r-1,c(4:6,10)] 
      }
      
      if(any(is.na(theta)) | theta[4] > 1/(r-1)) {
       # cat("mid,r,theta=",mid,r,theta,"\n")
        aic= 10^10
        
      }else{ 
        aic= com.aic.fast(xdat, r=r, mid=mid, 
                         theta=theta, npar, nhfix) 
      }
      
      if(r==1) {xtest[[mid]][R,11] = aic
      }else if(r >=2){xtest[[mid]][r-1,11] = aic}

    } # end for r
    
  } # end for mid

  z$opt_r = opt_r
  z$redtest = xtest
  z
}
#-----------------------------------------------------