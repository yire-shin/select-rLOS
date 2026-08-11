# Internal implementation for the msrlos package.
# Source: sel_rmed_batch_11Aug26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.


# Pohang.rainfall=readRDS("D:Pohang_rainfall")

#sel.rmed(xdat, sig.ed=0.05, h.fix=batch$h.fix,
#              medtest=batch)

#----------Algoritm r_Median ----------------------------------
sel.rmed = function(xdat, maxr=NULL, choose="median",
                         h.fix=NULL, num_inits=10,
                         sig.ed=0.05, qq=c(.98,.99,.995,.998),
                         true.para= NULL, dmin=0.03,
                         medtest=NULL){
  z=list()
  xdat= na.omit(xdat)
  nsam= dim(xdat)[1]
  R = maxr  
  
  if(is.null(maxr)) R = maxr= dim(xdat)[2]
  if(maxr < 2) stop("marx should be greater than or equal to 2")
  
  if(is.null(h.fix)){
    h.sel= h.select(xdat, sig.ed=sig.ed, dmin=dmin)
    h.fix= h.sel$h.fix
#    z$optr_rk4d = h.sel$optr
  }
#  z$h.fix = h.fix

  nhfix=length(h.fix)
  nmod= 6+length(h.fix)

  k3d.name=NULL
  if(!is.null(h.fix)) k3d.name= rep('rk3d',nhfix)
  k3d.name2=NULL
  if(!is.null(h.fix)) { 
    k3d.name2= paste("rk3d","(h=",round(h.fix,3),")",sep="") }
  
  name0= c('rglo','rggd', 'rgev', 'rk4d', 'rld','rgd')
  name.m = c(name0, k3d.name)
  name.m2 = c(name0, k3d.name2)
  
  opt_r = medtest$opt_r

     otr  = round(quantile(opt_r, probs=0.5, type=2,
                                na.rm=T)+1e-5, digits=0)
     otrv= min(otr,maxr)

  aic=rep(NA,nmod); 

  for(mid in 1:nmod){
    if(otr==1) {rcol=R
    }else if(otr >= 2) {rcol=otr-1}
    aic[mid]= medtest$redtest[[mid]][rcol,11]
  }

    last= which.min(aic[1:nmod])
    delta= aic-min(aic, na.rm=TRUE)
    bid.good= which(delta < 2)
    bestm.good= name.m2[bid.good]

  # z$maxr = maxr;# z$theta = theta
  # z$model.names = name.m2; z$optr = opt_r;
    
  z$algorithm="rmed"
  z$rhats =opt_r
  z$rstar= otr 
  # z$aic= aic; 
  z$mstar= last; 
  z$best.model= paste(name.m2[last],otr,sep=" ")
  
  # z$good.models= bestm.good
  # 
  theta=matrix(NA,nmod,4)
  for(mid in bid.good){
    if(otr==1) {rcol=R
    }else if(otr >= 2) {rcol=otr-1}
    mat = as.matrix(medtest$redtest[[mid]])
    theta[mid,1:4]= mat[rcol,c(4:6,10)]
  }
#  z$theta.good.model = theta[bid.good,]
  z$theta = theta[last,]
  z
}
