# Internal implementation for the msrlos package.
# Source: optr_find_24July26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.


#-----------------------------------------------------
find_optr=function(redtest,bigR,nmod, sigL=0.05, h.fix=NULL){
  
  opt_r= rep(bigR, nmod)
  nhfix=length(h.fix)
  
  for(mid in 1:nmod){
    
    R <- bigR
    if(mid >= 7 &  mid <= 6+nhfix){
      if(h.fix[mid-6] >= 1/(bigR-1) ){
        R = min(floor(1+ 1/h.fix[mid-6] -0.0001), bigR) }
    }
    opt_r[mid] = min( which(redtest$p.values <= sigL),
                      R )
    if(opt_r[mid]==R){
      opt_r[mid] = min( which(redtest$StrongStop <= sigL),
                        R )
      if(opt_r[mid]==R){
        opt_r[mid] = min( which(redtest$ForwardStop <= sigL),
                          R ) 
      }
    }
  } # end for mid
  opt_r
}
#------------------------------------------------------
one_optr=function(Edtest, mid, sigL=0.05, h.fix=NULL,
                        method='ed'){
  
  nhfix=length(h.fix)
  # if(method=='ed') {bigR= length(Edtest$r)
  # }else if(method=='pbscore'){bigR= length(Edtest$r)}
  R <- bigR <- length(Edtest$r)
  
  if(mid >= 7 &  mid <= 6+nhfix){
    if(h.fix[mid-6] >= 1/(bigR-1) ){
      R = min(floor(1+ 1/h.fix[mid-6] -0.0001), bigR) }
  }
  opt_r = min( which(Edtest$p.values <= sigL),
               R )
  if(opt_r==R){
    opt_r = min( which(Edtest$StrongStop <= sigL),
                 R )
    if(opt_r==R){
      opt_r = min( which(Edtest$ForwardStop <= sigL),
                   R ) 
    }
  }
  opt_r
}
#------------------------------------------------------   
#------------------------------------------------------
 # one_optr_batch(Edtest=batch, mid=9,
 #          h.fix=c(-0.3020, -0.2022, -0.0643))
#------------------------------------------------------
one_optr_batch=function(Edtest, mid, sigL=0.05, h.fix=NULL,
                  method='ed'){
  
  nhfix=length(h.fix)
  # if(method=='ed') {bigR= length(Edtest$r)
  # }else if(method=='pbscore'){bigR= length(Edtest$r)}
  R <- bigR <- length(Edtest$redtest[[mid]]$r)
  
  if(mid >= 7 &  mid <= 6+nhfix){
    if(h.fix[mid-6] >= 1/(bigR-1) ){
      R = min(floor(1+ 1/h.fix[mid-6] -0.0001), bigR) }
  }
  opt_r = min( which(Edtest$redtest[[mid]]$p.values <= sigL),
               R )
  if(opt_r==R){
    opt_r = min( which(Edtest$redtest[[mid]]$StrongStop <= sigL),
                 R )
    if(opt_r==R){
      opt_r = min( which(Edtest$redtest[[mid]]$ForwardStop <= sigL),
                   R ) 
    }
  }
  opt_r
}
#------------------------------------------------------   
#------------------------------------------------------
one_optr2=function(Edtest, mid, sigL=0.05, 
                   h.fix=NULL, method='ed'){
  
  sigL1 = sigL
  sigL2= sigL+0.05
  
  nhfix=length(h.fix)
  # if(method=='ed') {bigR= length(Edtest$r)
  # }else if(method=='pbscore'){bigR= length(Edtest$r)}
  R <- bigR <- length(Edtest$r)
  
  if(mid >= 7 &  mid <= 6+nhfix){
    if(h.fix[mid-6] >= 1/(bigR-1) ){
      R = min(floor(1+ 1/h.fix[mid-6] -0.0001), bigR) }
  }
  opt_r = min( which(Edtest$p.values <= sigL1),
               R )
  if(opt_r==R){
    opt_r = min( which(Edtest$p.values <= sigL2),
                 R )}
    if(opt_r==R){
      opt_r = min( which(Edtest$StrongStop <= sigL1),
                 R )}
      if(opt_r==R){
        opt_r = min( which(Edtest$StrongStop <= sigL2),
                     R )}
        if(opt_r==R){
           opt_r = min( which(Edtest$ForwardStop <= sigL1),
                   R ) }
           if(opt_r==R){
             opt_r = min( which(Edtest$ForwardStop <= sigL2),
                          R ) }
  opt_r
}
#------------------------------------ 
