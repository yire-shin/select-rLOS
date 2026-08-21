# Internal implementation for the msrlos package.
# Source: Alteria_new_batch_20Aug26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.



# alternating iterative algorithm
#------------------------------------------------------------------
Alteria.new = function(xdat, sig.ed=0.1, 
                       start=c("surv","rmed","rand_r","rand_m",
                               "fix_r","fix_m"),
                       mid.best=NULL, rhat.best=NULL,
                       dmin=0.03, h.fix=NULL,
                       altest=NULL, altering=TRUE){
  
  z=list()
  iter=1
  maxiter=15
  model= rhat= rep(NA,maxiter)
  
  z$algorithm = paste("alt ",start,sep="")
  R= maxr = ncol(xdat)
  z$ifail="converge"
  
  h.fix= altest$h.fix
  z$h.fix=h.fix
  nhfix=length(h.fix)
  nmod= 6+length(h.fix)
  
  if(start=="surv"){
    # x1=list()
    # x1 = sel.rmod.surv(xdat, maxr=NULL, h.fix=h.fix,
    #                    sig.ed=sig.ed, dmin=dmin,
    #                    survtest=altest, show=FALSE)
    
    model[iter]= mid.best  #x1$bid
    rhat[iter]=  rhat.best #x1$rhat
    
  }else if(start=="rmed"){
    # x2=list()
    # x2 = sel.rmod.rmed(xdat, maxr=NULL, h.fix=h.fix,
    #                    sig.ed=sig.ed, dmin=dmin,
    #                    medtest=altest)
    
    model[iter]= mid.best  #x2$bid
    rhat[iter]=  rhat.best #x2$rhat
    
  }else if(start=="rand_r"){
    
    rhat[iter] = sample.int(n=ncol(xdat), size=1)
    model[1]= 0
    
  }else if(start=="rand_m"){
    
    model[iter] =sample.int(n=nmod, size=1)
    rhat[1]= 0
    
  }else if(start=="fix_r"){
    rhat[1]= round(maxr/2+.01)
    model[1]= 0
    
  }else if(start=="fix_m"){
    model[1]= 4
    rhat[1]= 0
  }
  
  k3d.name=NULL; k3d.name= rep('rk3d',nhfix)
  name.m = c('rglo','rggd', 'rgev', 'rk4d', 'rld','rgd', 
             k3d.name)
  k3d.name2=NULL
  k3d.name2= paste(paste('rk3d','(h=',sep="_"),
                   round(h.fix,3),')',sep="") 
  name.m2 =c('rglo','rggd', 'rgev', 'rk4d', 'rld','rgd', 
             k3d.name2)
  
  # cat("iter, rhat[1], model[1]=",
  #     iter, rhat[iter], model[iter],"\n")
  aic=rep(NA,nmod)
  
  #----------------------- alternating --------------------
  if(altering==FALSE){
    z$rstar = rhat[iter]
    z$mstar = model[iter]
    z$model.names =  name.m2[z$mstar]
    z$ifail="converge"
    return(z)
    
  }else if(altering==TRUE){  #-------
    
    if(start=="surv" | start=="rand_r" | start=="fix_r"){
      
      while(iter < maxiter){
        
        iter=iter+1
        
        otr= rhat[iter-1]
        for(mid in 1:nmod){
          if(otr==1) {rcol=R
          }else if(otr >= 2) {rcol=otr-1}
          aic[mid]= altest$redtest[[mid]][rcol,11]
        }
        model[iter] = which.min(aic)
        
        # cat("iter, rhat[iter-1]--> model[iter]=",
        #    iter, rhat[iter-1], model[iter],"\n")
        # 
        if(model[iter]==model[iter-1]){
          rhat[iter]= z$rstar = rhat[iter-1]
          z$mstar = model[iter-1]
          z$model.names =  name.m2[z$mstar]
          #        z$model= model[1:iter]; z$rhat= rhat[1:iter]
          z$ifail="converge"
          return(z)
          
        }else if(model[iter] != model[iter-1]){
          
          rhat[iter]= altest$opt_r[ model[iter] ]
          
          # cat("iter, model[iter]--> rhat[iter]=",
          #     iter,  model[iter],rhat[iter], "\n")
          
          if(rhat[iter]==rhat[iter-1]){
            z$rstar = rhat[iter]
            z$mstar = model[iter]
            z$model.names =  name.m2[z$mstar]
            #          z$model= model[1:iter]; z$rhat= rhat[1:iter]
            z$ifail="converge"
            return(z)
          } #end if rhat[iter]
        } # end if model[iter]
        
        if(iter >= 3){
          prev= seq(1:(iter-2))
          if(any(rhat[iter] == rhat[prev]) 
             | any(model[iter] == model[prev])){
            z$ifail="cycle"
            iter.cy = iter
            # cat("iter, rhat, model, fail=",
            #     iter, rhat[iter],model[iter],z$ifail,"\n")
            break
          }
        } #end if iter
        
      } # end while
      
    }else if(start=="rmed" | start=="rand_m" | start=="fix_m"){
      
      while(iter < maxiter){
        
        iter=iter+1
        # edtest = multi.rEdtest.park(xdat, model=name.m[model[iter-1]],
        #                 h.fix= as.numeric(h.fix[model[iter-1]-6]), 
        #                 method="ed")
        # rhat[iter]= one_optr(edtest,mid=model[iter-1], sigL=sig.ed, 
        #                      h.fix=as.numeric(h.fix[model[iter-1]-6]))
        
        rhat[iter]= altest$opt_r[ model[iter-1] ]
        
        # cat("iter, model[iter-1]--> rhat[iter]=",
        #     iter,  model[iter-1],rhat[iter], "\n")
        
        if(rhat[iter]==rhat[iter-1]){
          z$rstar = rhat[iter]
          model[iter]= z$mstar = model[iter-1] 
          z$model.names =  name.m2[z$mstar]
          #          z$model= model[1:iter]; z$rhat= rhat[1:iter]
          z$ifail="converge"
          return(z)
          
        }else if(rhat[iter] != rhat[iter-1]){
          
          # aic= AIC.all.rfix(xdat, h.fix=h.fix, r=rhat[iter], name.m) 
          
          otr= rhat[iter]
          for(mid in 1:nmod){
            if(otr==1) {rcol=R
            }else if(otr >= 2) {rcol=otr-1}
            aic[mid]= altest$redtest[[mid]][rcol,11]
          }
          model[iter] = which.min(aic)
          
          # cat("iter, rhat[iter]--> model[iter]=",
          #     iter, rhat[iter], model[iter],"\n")
          
          #         model[iter] = which.min(x$aic[,rhat[iter]])  # need modif
          
          if(model[iter]==model[iter-1]){
            z$rstar = rhat[iter]
            z$mstar = model[iter]
            z$model.names =  name.m2[z$mstar]
            #            z$model= model[1:iter]; z$rhat= rhat[1:iter]
            z$ifail="converge"
            return(z)
          } #end if model
        } # end if rhat
        
        if(iter >= 3){
          prev= seq(1:(iter-2))
          if(any(rhat[iter] == rhat[prev]) 
             | any(model[iter] == model[prev])){
            z$ifail="cycle"
            iter.cy = iter
            # cat("iter, rhat, model, fail=",
            #     iter, rhat[iter],model[iter],z$ifail,"\n")
            break
          }
        } #end if iter
      } # end while
      
    } # end if start
    
    if(z$ifail=="cycle"){
      
      cycle = solve.cycle(xdat, start, h.fix=h.fix, sig.ed, 
                          dmin, name.m, rhat=rhat, 
                          model=model, iter.cy=iter.cy,
                          predt=altest)
      
      z$rstar = cycle$rstar
      z$mstar = cycle$mstar
      z$model.names =  name.m2[cycle$mstar]
      z$ifail = "cycle_solved"
      
    }else{
      
      # z$model= model[1:iter]
      # z$rhat= rhat[1:iter]
      z$ifail="No sol: reach max iteration"
    }
    z
  } #end if alter
}
#-----------------------------------------------------------
#------------------------------------------------
solve.cycle= function(xdat, start, h.fix=NULL, sig.ed, 
                      dmin, name.m, rhat=NULL, 
                      model=NULL, iter.cy=NULL,
                      predt=NULL){
  
  cy=list()
  max.ic= 2*iter.cy-1  
  rcy=mcy=rep(0, max.ic)
  comb= matrix(NA,max.ic,2)
  prev= c(1:(max.ic-1))
  
  if( start=="surv" ){
    for(it in 1:max.ic){   
      
      rcy[it] = rhat[round(it/2+.01)]
      mcy[it] = model[round(it/2+0.51)]
      #       comb[it,1:2] = c(rcy[it],mcy[it])
    } # end for
    
  }else if(start=="rand_m" | start=="rmed" | start=="fix_m"){
    
    for(it in 1:max.ic){  
      mcy[it] = model[round(it/2+.01)]
      rcy[it] = rhat[round(it/2+0.51)]
      #        comb[it,1:2] = c(rcy[it],mcy[it])
    }
    
  }else if( start=="rand_r" | start =="fix_r"){
    
    for(it in 1:max.ic){  
      mcy[it] = model[round(it/2+.51)]
      rcy[it] = rhat[round(it/2+0.01)]
      comb[it,1:2] = c(rcy[it],mcy[it])
    }
    if( any(comb[max.ic,1:2]==comb[prev,1:2]) ) {
      rcy[max.ic] =0; mcy[max.ic]=0
    }
    
  }# end if start
  
  aic1=rep(1e10,max.ic) 
  for (it in 1:max.ic){
    
    if( (rcy[it] !=0 ) & (mcy[it] != 0) ){
      
      aic1[it]= AIC1.model.rfix(xdat, h.fix=h.fix, 
                                r=rcy[it], mid=mcy[it], name.m,
                                edtest1=predt)
      
      # cat("it, rcy, mcy, aic12=", it,rcy[it],mcy[it],
      #     aic1[it],"\n")
    }
    
  } # end for
  
  star= which.min(aic1)
  
  cy$mstar=mcy[star]
  cy$rstar=rcy[star]
  cy$ifail="cycle_solved"
  cy
}
#----------------------------------------
# #----------------------------------------------------------
# decide.MR = function(xdat, name.m, h.fix=NULL, sig.ed,x1,x2){
#   
#   cat("x1$bid, x2$bid, x1$rhat, x2$rhat=",
#       x1$bid, x2$bid, x1$rhat, x2$rhat,"\n")
#   
#   if(x1$bid==x2$bid & x1$rhat==x2$rhat){
#     
#     cat("cycle, model.cy, rhat.cy=",
#         x1$bid, x1$rhat, "\n")
#     
#     z$rstar = x1$rhat
#     z$mstar = x1$bid
#     z$model.names =  name.m2[z$mstar]
#     z$ifail="cycle_but_opt_result"
#     return(z)
#   }
#   
#   if(x1$bid==x2$bid){
#     edtest = multi.rEdtest.park(xdat, model=name.m[x1$bid],
#                                 h.fix= as.numeric(h.fix[x1$bid-6]), 
#                                 method="ed")
#     rhat.cy= one_optr(edtest,mid=x1$bid, sigL=sig.ed, 
#                       h.fix=as.numeric(h.fix[x1$bid-6]))
#     
#     cat("cycle, model.cy, rhat.cy=",
#         x1$bid, rhat.cy, "\n")
#     
#     z$rstar = rhat.cy
#     z$mstar = x1$bid
#     z$model.names =  name.m2[z$mstar]
#     z$ifail="cycle_but_opt_result"
#     return(z)
#     
#   }else if(x1$rhat==x2$rhat){
#     
#     model.cy = which.min(x1$aic[,x1$rhat])  
#     
#     cat("cycle, model.cy, rhat.cy=",
#         model.cy, x1$rhat, "\n")
#     
#     z$rstar = x1$rhat
#     z$mstar = model.cy
#     z$model.names =  name.m2[z$mstar]
#     z$ifail="cycle_but_opt_result"
#     return(z)
#     
#   }else{
#     z$ifail="still_cycle"
#     return(z)
#   }
# }
# #---------------------------------------------------------
