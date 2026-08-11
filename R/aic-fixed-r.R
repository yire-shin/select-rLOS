# Internal implementation for the msrlos package.
# Source: AIC1_model_rfix.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.


# for given rfix, compute AIC-2 (for 2nd los) for a fixed rlos model
#-----------------------------------------------------------------
AIC1.model.rfix = function(xdat, h.fix=NULL, r=NULL, 
                           mid=NULL, name.m,
                           edtest1=NULL){
  
  nhfix=length(h.fix)
  nmod= 6+length(h.fix)
  opt_r=r
  theta=matrix(NA, nmod,4)
 
    if(r==1) {rcol=ncol(xdat)
    }else if(r >= 2) {rcol=r-1}
    mat = as.matrix(edtest1$redtest[[mid]])
    theta[mid,1:4]= mat[rcol,c(4:6,10)]

  npar=c(3,3,3,4,2,2,rep(3,nhfix))
  
#  for(mid in 1:nmod){
    aic= com.aic(xdat, r=1, mid=mid, theta=theta, npar, nhfix) 
    #  cat("mid,aic=",mid,aic[mid],"\n")
#  }
  return(aic)
}