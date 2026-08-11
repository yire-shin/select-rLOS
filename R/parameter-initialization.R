# Internal implementation for the msrlos package.
# Source: theta.req1_24July26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.


#----------------------------------------------------------------------
theta.req1 = function(xdat, surv=NULL, h.fix=NULL, 
                      num_inits=5){
  
  nhfix=length(h.fix)
  nmod= 6 + nhfix
  newdt= as.matrix(xdat[,1])
  
  name.m = c('rglo','rggd', 'rgev', 'rk4d', 'rld','rgd')
 
  work.fit=paste(name.m,'.fit.park',sep="")

  retheta=matrix(-0.0001,nmod,4)

  for (mid in surv){
    
    theta=rep(-0.0001,4)
    
    if(mid >= 7 &  mid <= 6+nhfix){
      theta = rk3d.fit.park(newdt, r=1, h.fix=h.fix[mid-6],
                            num_inits=num_inits, reltol=1e-5)$mle
      
      if(theta[3] < -0.5) {
        theta= rk3d.fit.park(newdt,r=1,show=F,
                             h.fix=h.fix[mid-6], penk="CD",
                             num_inits=num_inits, reltol=1e-4)$mle }
      
    }else if(mid <= 6){
      theta= match.fun(work.fit[mid])(newdt, r=1, 
                                      num_inits=num_inits, show=F,
                                      reltol=1e-5)$mle
      
      if(mid==4){
        if(theta[4] > 1.0 | theta[3] < -0.5 |
           theta[3] > 0.5){
          theta= match.fun(work.fit[mid])(newdt,r=1,show=F,
                                          penk="CD", penh="MSa", 
                                          num_inits=num_inits, 
                                          reltol=1e-4)$mle }
      } #end if mid=4
      
      if(mid==1 & theta[3] < -0.5){

        theta= rk3d.fit.park(newdt,r=1,show=F,
                             penk="CD", h.fix= -1.0,
                             num_inits=num_inits,
                             reltol=1e-4)$mle }
      
    } #end if mid>7
    
    if(mid==6) {theta[3]= 0.0001; theta[4]= -1e-4}
    if(mid==5) {theta[3]= 0.0001; theta[4]= -.9999}
    if(mid >= 7 &  mid <= 6+nhfix) theta[4]= h.fix[mid-6]
    if(mid > 6+nhfix) theta[4] = -.0001 
    if(mid ==3) theta[4]= -1e-4
    if(mid ==1) theta[4] = -.9999
    if(mid ==2) {theta[4]= theta[3]; theta[3]= -1e-4}
    
    retheta[mid,1:4]= theta[1:4]

  } # end for mid 
  
 retheta
}
