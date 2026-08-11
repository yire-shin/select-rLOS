# Internal implementation for the msrlos package.
# Source: com.aic_10Aug26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.


#----------------------------------------------------------------
com.aic = function(data,r=NULL,mid=NULL,theta=NULL,
                   npar=NULL,nhfix){

  if(mid==1 | mid==5) {
    temp = -2* sum(rgloLh(as.matrix(data[,1:r]),
                          theta[mid,1:3]), na.rm=T) # -2*llh

  }else if(mid==2 | mid==6){
    temp = -2* sum(rggdLh(as.matrix(data[,1:r]),
                          theta[mid,c(1,2,4)]), na.rm=T) # -2*llh

  }else if(mid==3 | mid > 6+nhfix){
    temp = -2* sum(eva::dgevr(as.matrix(data[,1:r]),
                         loc= theta[mid,1],
                         scale= theta[mid,2],
                         shape= -theta[mid,3],
                         log.d=TRUE), na.rm=T)
    # coles style

  }else if(mid==4 | (mid >= 7 & mid <= 6+nhfix)){
    temp =  -2* sum(rk4d.lh(as.matrix(data[,1:r]), 
                            theta[mid,]), na.rm=T) # -2*nllh
  }
  
  if(is.na(temp) | temp==-Inf | temp == Inf) temp = 10^6
  temp =  temp + 2* npar[mid]
  return(temp)
}

#----------------------------------------------------------------
com.aic.fast = function(data,r=NULL,mid=NULL,theta=NULL,
                   npar=NULL,nhfix){
  
  workd= as.matrix(data[,1:r])
  
  if(mid==1 | mid==5) {
    temp = -2* sum(rgloLh(workd,
                          theta[1:3]), na.rm=T) # -2*llh
    
  }else if(mid==2 | mid==6){
    temp = -2* sum(rggdLh(workd,
                          theta[c(1,2,4)]), na.rm=T) # -2*llh
    
  }else if(mid==3 | mid > 6+nhfix){
    temp = -2* sum(eva::dgevr(workd,
                         loc= theta[1],
                         scale= theta[2],
                         shape= -theta[3],
                         log.d=TRUE), na.rm=T)
    # coles style
    
  }else if(mid==4 | (mid >= 7 & mid <= 6+nhfix)){
    temp =  -2* sum(rk4d.lh(workd, 
                            theta[1:4]), na.rm=T) # -2*nllh
  }
  
  if(is.na(temp) | temp==-Inf | temp == Inf) temp = 10^6
  temp =  temp + 2* npar[mid]
  return(temp)
}
