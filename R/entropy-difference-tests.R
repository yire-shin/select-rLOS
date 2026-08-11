# Internal implementation for the msrlos package.
# Source: EDtest_subs_16July26.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

# multi.rEdtest.park(data, model="rggd")
# multi.rEdtest.park(data, model="rglo")
# multi.rEdtest.park(data, model="rgev")
# multi.rEdtest.park(data, model="rk4d")
# multi.rEdtest.park(data, model="rgd")
# multi.rEdtest.park(data, model="rld")
# rlarg.fit

#----------------------------------------
rgev2Ed.park = function (data, k.fix=NULL, theta = NULL,
                         num_inits=5)
{
  # rgev.fit for k.fix cases; 
  # here k.fix is only a scalar value, not a vector
  
  if(is.null(k.fix) | is.na(k.fix)) stop("k.fix should not be null nor NA")
  
  data <- as.matrix(data)
  R <- ncol(data)
  if (R == 1) 
    stop("R must be at least two")
  n <- nrow(data)
  if (is.null(theta)) {
    
#    cat("k.fix in rgev2Ed.park =", k.fix,"\n")
    
    y <- tryCatch(rgev2.fit.park(data, k.fix=k.fix,
                                 num_inits=num_inits,
                                 reltol=1e-5)
                  # y = eva::gevrFit(data, method = "mle")
    )
    if (is.null(y)) {
      y <- tryCatch(rgev2.fit.park(data, num_inits=num_inits,
                                   k.fix=k.fix, reltol=1e-4)
                    #eva::gevrFit(data, method = "mle")
      )}
    
    if (is.null(y)) {
      warning("Maximum likelihood failed to converge at initial step")
      y$mle =lmomco::pargev(lmomco::lmoms(data[,1]))$para
      theta = y$mle
      theta[3]= -y$mle[3]   # coles style
    }else{ 
      theta <- y$mle        # hosking style     #y$par.ests
      theta[3]= -y$mle[3]   # coles style
    }
  }
  Diff <- eva::dgevr(data[, 1:R], loc = theta[1], scale = theta[2], 
                shape = theta[3], log.d = TRUE) - eva::dgevr(data[, 1:(R - 
                                                                    1)], loc = theta[1], scale = theta[2], shape = theta[3], 
                                                        log.d = TRUE)
  
  EstVar <- sum((Diff - mean(Diff))^2)/(n - 1)
  FirstMom <- -log(theta[2]) - 1 + (1 + theta[3]) * digamma(R)
  Diff <- sum(Diff)/n
  Diff <- sqrt(n) * (Diff - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - pnorm(abs(Diff)))
  
  theta1= theta
  theta1[3]= -theta[3]   # hosking style
  
  out <- list(statistics = as.numeric(Diff), p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(FirstMom))
  out
}
#-------------------------------------------------------   
#-------------------------------------------------
rk3dEd.park = function (data, h.fix=NULL, par=NULL,
                        num_inits=5) 
{
  if(is.null(h.fix) | is.na(h.fix)){
    stop("h.fix should not be null nor NA")
  }
  R <- ncol(data)
  if(h.fix >= 1/(R-1)) stop("h.fix should be less 1/(R-1)")
  
  y <- suppressWarnings(rk3d.fit.park(data, h.fix=h.fix, show = F, 
                                num_inits=num_inits, reltol=1.e-5))
  theta1 <- y$mle
  
  if(theta1[3] < -0.5 | y$conv !=0 ) {
    y <- suppressWarnings(rk3d.fit.park(data, h.fix=h.fix, 
                              show = F, num_inits=num_inits, 
                              penk="CD", low.xi = -0.6, reltol=1.e-4))
    theta1 <- y$mle
  }
  
  #  cat("theta 1=", theta1, "\n")
  if(is.null(theta1)){
    # y <- suppressWarnings(rk3d.fit.park(data, h.fix=h.fix, 
    #                               show = F, num_inits=num_inits, 
    #                               penk="CD",reltol=1.e-4))
    # theta1 <- y$mle
    # #    cat("theta 2=", theta1, "\n")
    # if(is.null(theta1)){
      y <- suppressWarnings(rgev.fit.park(data, show = F, 
                                  num_inits=num_inits,reltol=1.e-4))
      theta1 <- y$mle
      theta1[4] = h.fix
   # }
  }
  
  mu <- theta1[1]
  sc <- theta1[2]
  xi <- theta1[3]
  h <- h.fix
  theta1[4] = h.fix

  nr <- nrow(data)
  Diff <- ( rk4d.lh(data[, 1:R], theta1) 
          - rk4d.lh(as.matrix(data[,1:(R - 1)], ncol =R-1),theta1) )
  
  EstVar <- sum((Diff - mean(Diff))^2)/(nr - 1)
  if (h > 0) {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma(ar) - digamma(ar + 
                                                        R))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma(ar1) - digamma(ar1 + 
                                                               R - 1))
    term5 <- (1 - xi) * ((digamma(R) - digamma(ar + R)) - 
                           log(h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  else {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma((1/-h) + R) - digamma(1/-h))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma((1/-h) + R - 
                                                1) - digamma(1/-h))
    term5 <- (1 - xi) * (digamma(R) - digamma(1/-h) - log(-h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
 
  Diff1 <- sum(Diff)/nr
  Stat <- sqrt(nr) * (Diff1 - eta)/sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Stat)))
  out <- list(statistics = as.numeric(Stat), 
              p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(Diff1))
  out
}
#------------------------------------------------------------
# h.fix=hw.fix =0.15; model='rgev2'; data=xdat
#------------------------------------------------------------
multi.rEdtest.park =function(data=NULL, model=NULL, method="ed",
                             h.fix=NULL, k.fix=NULL, par = NULL,
                             num_inits=10) 
{
  # here k.fix and h.fix are only scalar values, not a vector

  modelEd=paste(model,"Ed.park", sep="")
  hw.fix= as.numeric(h.fix)
  
  data <- as.matrix(na.omit(data))
  R <- ncol(data)
  if(model=="rk3d"){ 
    if(hw.fix >= 1/(ncol(data)-1) ){
      R = min(floor(1+ 1/hw.fix -0.0001), ncol(data)) }
  }

  if(R==1){ 
     result <- matrix(0, 1, 10) 
     colnames(result) <- c("r", "p.values", "statistic", "est.loc", 
                           "est.scale", "est.shape", "ybar", "ForwardStop",
                           "StrongStop", "est.shape2")
     return(as.data.frame(result))
  }
  result <- matrix(0, R - 1, 10)
 
  if(model=="rgev"){
    #result= gevrSeqTests.park(data, simulation=TRUE) # for simulation study
    result= gevrSeqTests.park(data, bootnum=500,
                              method=method,  
                              num_inits=num_inits)  # for real data analysis
    # result= gevrSeqTests(data)

  }else{
    for (i in 2:R) {
      result[i - 1, 1] <- i
      if(model=="rk3d"){
        fit = rk3dEd.park(data[, 1:i], h.fix=hw.fix,
                          num_inits=num_inits)
      }else if(model=="rgev2"){
        fit = rgev2Ed.park(data[, 1:i], k.fix=k.fix,
                           num_inits=num_inits)
      }else{
        fit <- match.fun(modelEd)(data[, 1:i], par,
                                  num_inits=num_inits)
      }
      
      if(is.na(fit$p.value)) fit$p.value=0
      
      result[i - 1, 2] <- fit$p.value
      result[i - 1, 3] <- fit$statistics
      result[i - 1, 7] <- fit$ybar
      
      if(model=="rk4d" | model=="rk3d"){
        result[i - 1, 4:6] <- fit$theta[1:3]  # hosking style
        result[i - 1, 10] = fit$theta[4] 
        
      }else if(model=="rglo" | model=="rld"){
        result[i - 1, 4:6] <- fit$theta   # hosking style 
        result[i - 1, 10] = -0.9999
        
      }else if(model=="rgd" | model=="rggd"){
        result[i - 1, 4:5] <- fit$theta[1:2]   
        result[i - 1, 6] <-  1e-4
        result[i - 1, 10] = fit$theta[3]   # hosking style 
      }

    result[, 8] <- rev(eva::pSeqStop(rev(result[, 2]))$ForwardStop)
    result[, 9] <- rev(eva::pSeqStop(rev(result[, 2]))$StrongStop)
    
   } # end for
  } # end if model
  
  coll = c("r", "p.values", "statistic", "est.loc", 
           "est.scale", "est.shape", "ybar", "ForwardStop",
           "StrongStop", "est.shape2")
  colnames(result) <- c(coll)
  
  as.data.frame(result)
}
#--------------------------------------------------------
gevrEd.park1 =function (data, theta = NULL,num_inits=5) 
{
  data <- as.matrix(data)
  R <- ncol(data)
  if (R == 1) 
    stop("R must be at least two")
  n <- nrow(data)
  ifail=0
  
  if (is.null(theta)) {
    y <- tryCatch(eva::gevrFit(data, method = "mle"), error = function(w) {
      return(NULL)
    }, warning = function(w) {
      return(NULL)
    })
    if (is.null(y)) {
      warning("Maximum likelihood failed to converge at initial step in gevrFit")
      out=list(); out$ifail= 1
      
      y= rgev.fit.park(data,reltol=1e-5,
                       num_inits=num_inits)  # park
      theta = y$mle           # park  y$mle is hosking style
      theta[3]= -theta[3]
      
      # return(out)           # park
      
    }else{
      theta <- y$par.ests  # coles style
#      theta[3]= -theta[3]  # hosking style
    }
  }
  
  if(theta[2] < 1e-4) theta[2]= 1e-4
  Diff1 <- eva::dgevr(data[, 1:R], loc = theta[1], scale = theta[2], 
                shape = theta[3], log.d = TRUE)
  
  Diff2= eva::dgevr(data[, 1:(R-1)], loc = theta[1], scale = theta[2], 
                shape = theta[3], log.d = TRUE)
  
  Diff= Diff1 - Diff2
    
  EstVar <- sum((Diff - mean(Diff))^2)/(n - 1)
  FirstMom <- -log(theta[2]) - 1 + (1 +theta[3]) * digamma(R)
  Diff <- sum(Diff)/n
  Diff <- sqrt(n) * (Diff - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - pnorm(abs(Diff)))

  out <- list(as.numeric(Diff), as.numeric(p.value), theta, ifail)
  names(out) <- c("statistic", "p.value", "theta", "ifail")
  return(out)
}
#----------------------------------------
gevrEd.park2 = function (data, theta = NULL, num_inits=5)
{
  data <- as.matrix(data)
  R <- ncol(data)
  if (R == 1) 
    stop("R must be at least two")
  n <- nrow(data)

  if (is.null(theta)) {
    y <- tryCatch(rk3d.fit.park(data, h.fix= -0.0001,
                                num_inits=num_inits,
                                penk="CD",
                                reltol= 1e-5)
                  # y = eva::gevrFit(data, method = "mle")
    )
    if (is.null(y)) {
      y <- tryCatch(rk3d.fit.park(data, num_inits=num_inits,
                                  h.fix= -0.0001, penk="CD",
                                  reltol= 1e-4)
    )}
    
    if (is.null(y)) {
      warning("Maximum likelihood failed to converge at initial step")
      y$mle =lmomco::pargev(lmomco::lmoms(data[,1]))$para
      theta = y$mle
    }else{ 
      theta <- y$mle[1:3]        # hosking style     #y$par.ests
      theta[3]= -y$mle[3]        # coles style
    }
  }
  Diff <- eva::dgevr(data[, 1:R], loc = theta[1], scale = theta[2], 
                shape = theta[3], log.d = TRUE) - eva::dgevr(data[, 1:(R - 
                                                                    1)], loc = theta[1], scale = theta[2], shape = theta[3], 
                                                        log.d = TRUE)

  EstVar <- sum((Diff - mean(Diff))^2)/(n - 1)
  FirstMom <- -log(theta[2]) - 1 + (1 + theta[3]) * digamma(R)
  Diff <- sum(Diff)/n
  Diff <- sqrt(n) * (Diff - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - pnorm(abs(Diff)))
  out <- list(as.numeric(Diff), as.numeric(p.value), theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
#-------------------------------------------------------   
gevrSeqTests.park = function (data, bootnum = NULL, 
                              method = c("ed", "pbscore", "multscore"), 
                              information = c("expected", "observed"), 
                              allowParallel = FALSE, 
                              numCores = 1, num_inits=5,
                              low.xi= -0.7) 
{
  data <- as.matrix(data)
  R <- ncol(data)
  method <- match.arg(method)
  
  if (method != "ed") {
    if (is.null(bootnum)) 
      stop("Must enter the number of bootstrap replicates!")
    information <- match.arg(information)
    result <- matrix(0, R, 10)
    for (i in 1:R) {
      result[i, 1] <- i
      if (method == "multscore") 
        fit <- eva::gevrMultScore(data[, 1:i], bootnum, information)
      if (method == "pbscore") 
        fit <- eva::gevrPbScore(data[, 1:i], bootnum, information, 
                           allowParallel, numCores)
      result[i, 2] <- fit$p.value
      result[i, 3] <- fit$statistic
     # result[i, 6:8] <- fit$theta
      result[i, 4:5] <- fit$theta[1:2]
      result[i, 6] <-  -fit$theta[3]   #----- hosking style
      result[i,10] = 1e-4
    }
  }
  else {
    if (R == 1)
      stop("R must be at least two")
    result <- matrix(0, R - 1, 10)
    for (i in 2:R) {
      result[i - 1, 1] <- i
      fit <- gevrEd.park1(data[, 1:i])

      if(fit$ifail==1){
        cat("gevrEd.park1: ifail, theta=", fit$ifail, fit$theta,"\n")
        mle <- rk3d.fit.park(data[, 1:i], h.fix=-0.001, 
                             penk="CD",
                             num_inits=num_inits)$mle

        mle[3]= -mle[3]
        fit = gevrEd.park1(data[, 1:i], theta=mle[1:3],
                           num_inits=num_inits)
      }
      result[i - 1, 2] <- fit$p.value
      result[i - 1, 3] <- fit$statistic
      result[i - 1, 4:5] <- fit$theta[1:2]
      result[i - 1, 6] <-  -fit$theta[3]   #----- hosking style
      result[i-1,10] = 1e-4
    }
  }
  result[,7] <- 0.0
  result[, 8] <- rev(eva::pSeqStop(rev(result[, 2]))$ForwardStop)
  result[, 9] <- rev(eva::pSeqStop(rev(result[, 2]))$StrongStop)
  colnames(result) <- c("r", "p.values",
                        "statistic", "est.loc", "est.scale",
                        "est.shape", "ybar", 
                        "ForwardStop", "StrongStop", 
                        "est.shape2")
  as.data.frame(result)
}

#-------------------------------------------------
rk4dEd.park = function (data, par=NULL,
                        num_inits=5) 
{
  y <- suppressWarnings(rk4d.fit.park(data, show = F, 
                                      num_inits=num_inits,
                                      reltol=1.e-5))
  theta1 <- y$mle
  
  if(theta1[3] < -0.5 | theta1[4] > 1.0 | y$conv !=0){
    y= rk4d.fit.park(data, show = F, num_inits=num_inits, 
                     penk="CD",penh="MSa", reltol=1.e-4,
                     low.xi=-0.6)
    theta1 <- y$mle
  }

  if(is.null(theta1)){
    y <- suppressWarnings(rk4d.fit.park(data, show = F, 
                                        num_inits=num_inits, 
                          penk="CD",penh="MSa", reltol=1.e-4,
                          low.xi=-0.5))
    theta1 <- y$mle

    if(is.null(theta1)){
      y <- suppressWarnings(rgev.fit.park(data, show = F, 
                                          num_inits=num_inits))
      theta1 <- y$mle
      theta1[4] = 0.001
    }
  }
  
  mu <- theta1[1]
  sc <- theta1[2]
  xi <- theta1[3]
  h <- theta1[4]
  R <- ncol(data)
  nr <- nrow(data)
  Diff <- ( rk4d.lh(data[, 1:R], theta1) 
          - rk4d.lh(as.matrix(data[,1:(R-1)],ncol=R-1),theta1) )
 
  EstVar <- sum((Diff - mean(Diff))^2)/(nr - 1)
  if (h > 0) {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma(ar) - digamma(ar + 
                                                        R))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma(ar1) - digamma(ar1 + 
                                                               R - 1))
    term5 <- (1 - xi) * ((digamma(R) - digamma(ar + R)) - 
                           log(h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  else {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma((1/-h) + R) - digamma(1/-h))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma((1/-h) + R - 
                                                1) - digamma(1/-h))
    term5 <- (1 - xi) * (digamma(R) - digamma(1/-h) - log(-h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  Diff1 <- sum(Diff)/nr
  Stat <- sqrt(nr) * (Diff1 - eta)/sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Stat)))
  out <- list(statistics = as.numeric(Stat), p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(Diff1))
  out
}
#------------------------------------------------------------
#--------------------------------------------------
rggdEd.park= function (data, par=NULL, num_inits=5) 
{
  y <- rggd.fit.park(data, show = F, num_inits=num_inits,
                     reltol=1.e-5)
  theta1 <- y$mle
  #cat("theta 1=", theta1, "\n")
  if(is.null(theta1)){
    y <- rggd.fit.park(data, show = F, num_inits=num_inits, 
                       maxit=100, reltol=1.e-4)
    theta1 <- y$mle
    #cat("theta 4=", theta1, "\n")
    if(is.null(theta1)){
      y <- rgd.fit.park(data, show = F, num_inits=num_inits,
                        reltol=1.e-4)
      theta1 <- y$mle
      theta1[3] = -0.001
    }
  }

  mu <- theta1[1]
  sc <- theta1[2]
  h <- theta1[3]
  R <- ncol(data)
  nr <- nrow(data)
  Diff <- rggdLh(data[, 1:R], theta1) - rggdLh(as.matrix(data[, 
                                                              1:(R - 1)], ncol = R - 1), theta1)
  EstVar <- sum((Diff - mean(Diff))^2)/(nr - 1)
  ri <- (R - seq(1:(R)))
  cr <- prod((1 - ri * h))
  ri1 <- (R - 1 - seq(1:(R - 1)))
  cr1 <- prod((1 - ri1 * h))
  if (h > 0) {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma(ar) - digamma(ar + 
                                                        R))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma(ar1) - digamma(ar1 + 
                                                               R - 1))
    term5 <- ((digamma(R) - digamma(ar + R)) - log(h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  else {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma((1/-h) + R) - digamma(1/-h))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma((1/-h) + R - 
                                                1) - digamma(1/-h))
    term5 <- (digamma(R) - digamma(1/-h) - log(-h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  Diff1 <- sum(Diff)/nr
  Diff <- sqrt(nr) * (Diff1 - eta)/sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Diff)))
  out <- list(statistics = as.numeric(Diff), p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(Diff1))
  out
}

#----------------------------------------------------
rgloEd.park = function (data, par = NULL,num_inits=5) 
{
  if (is.null(par)) {
    y <- rglo.fit.park(data, show = F, num_inits=num_inits)
    theta1 <- y$mle
  }else{
    theta1 <- par
  }
  R <- ncol(data)
  nr <- nrow(data)
  Diff <- rgloLh(data[, 1:R], theta1) - rgloLh(as.matrix(data[, 
                                                              1:(R - 1)], ncol = R - 1), theta1)
  
  FirstMom <- -log(theta1[2]) + log(R) - ((R + 1)/(R)) + theta1[3] * 
    (digamma(1) - digamma(R))
  EstVar <- sum((Diff - FirstMom)^2)/(nr)
  
  Diff1 <- sum(Diff)/nr
  Diff <- sqrt(nr) * (Diff1 - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Diff)))
  out <- list(statistics = as.numeric(Diff), p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(Diff1))
  out
}

#-----------------------------------------------
rldEd.park = function (data, par = NULL, num_inits=5) 
{

    y <- rld.fit.park(data, show = F, num_inits=num_inits)
    theta1 <- c(y$mle, 0.0001)

  R <- ncol(data)
  nr <- nrow(data)
  Diff <- rgloLh(data[, 1:R], theta1) - rgloLh(as.matrix(data[, 
                                                              1:(R - 1)], ncol = R - 1), theta1)
  
  FirstMom <- -log(theta1[2]) + log(R) - ((R + 1)/(R)) + theta1[3] * 
    (digamma(1) - digamma(R))
  EstVar <- sum((Diff - FirstMom)^2)/(nr)
  
  Diff1 <- sum(Diff)/nr
  Diff <- sqrt(nr) * (Diff1 - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Diff)))
  out <- list(statistics = as.numeric(Diff), p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(Diff1))
  out
}
#----------------------------------------
rgdEd.park =function (data, par=NULL,num_inits=5) 
{

  R <- ncol(data)
    y <- rgd.fit.park(data, r=R, show = F, num_inits=num_inits)
    theta1 <- c(y$mle, -0.0001)

  mu <- theta1[1]
  sc <- theta1[2]
  h <- theta1[3]
  nr <- nrow(data)
  Diff <- rggdLh(data[, 1:R], theta1) - rggdLh(as.matrix(data[, 
                                                              1:(R - 1)], ncol = R - 1), theta1)
  EstVar <- sum((Diff - mean(Diff))^2)/(nr - 1)
  ri <- (R - seq(1:(R)))
  cr <- prod((1 - ri * h))
  ri1 <- (R - 1 - seq(1:(R - 1)))
  cr1 <- prod((1 - ri1 * h))
  if (h > 0) {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma(ar) - digamma(ar + 
                                                        R))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma(ar1) - digamma(ar1 + 
                                                               R - 1))
    term5 <- ((digamma(R) - digamma(ar + R)) - log(h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  else {
    ar <- (1 - (R - 1) * h)/h
    ar1 <- (1 - (R - 2) * h)/h
    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h)/h) * (digamma((1/-h) + R) - digamma(1/-h))
    term4 <- ((1 - (R - 1) * h)/h) * (digamma((1/-h) + R - 
                                                1) - digamma(1/-h))
    term5 <- (digamma(R) - digamma(1/-h) - log(-h))
    eta <- -term1 + term2 + term3 - term4 + term5
  }
  Diff1 <- sum(Diff)/nr
  Diff <- sqrt(nr) * (Diff1 - eta)/sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Diff)))
  out <- list(statistics = as.numeric(Diff), p.value = as.numeric(p.value), 
              theta = theta1, ybar = as.numeric(Diff1))
  out
}

#---------------------------------------------------------
rk4d.lh <- function(data, par) {
  
  
  R <-ncol(data)
  nr<-nrow(data)
  
  mu <-par[1]
  sc <-par[2]
  xi <-par[3]
  h  <-par[4]
  
  ri  <- (R - seq(1:(R)))
  cr  <- (1 - ri * h)
  
  y <- 1 - xi * (data - mu) / sc
  f <- 1 - h * (1 - xi * (data[,R] - mu) / sc)^(1/xi)
  
  y <- log(sc) + (1 - 1/xi) * log(y) - log(cr)
  y <- rowSums(y, na.rm=TRUE)
  
  log.den3 <- (R * h - 1) / h * log(f) + y
  
  return(-log.den3)
  
}
#----------------------------------------------------
rggdLh <- function(data,par) {
  
  R <-ncol(data)
  nr<-nrow(data)
  
  mu <-par[1]
  sc <-par[2]
  h  <-par[3]
  
  ri  <- (R-seq(1:(R))) # r-i
  cr  <- (1-ri*h)     # c_r
  
  y <- exp(-(data - mu)/sc)
  f <- 1 - h * exp(-(data[,R] - mu)/sc)
  
  
  y <- log(sc) - log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)
  
  log.den3 = (R*h - 1)/h * log(f) + y
  
  return(-log.den3)
  
}
#-------------------------------------------------------
rgloLh <- function(data,par) {
  
  R <-ncol(data)
  nr<-nrow(data)
  
  mu <-par[1]
  sc <-par[2]
  xi <-par[3]
  
  ri  <- (R-seq(1:(R))) # r-i
  cr  <- (1+ri)    # c_r
  
  log.den3 = -R*log(sc) + sum(log(cr)) + (1+R) * log(1 / (1+((1-xi*(data[,R]-mu)/sc)^(1/xi))))+
    rowSums(((1/xi) -1)*log(1-xi*((data-mu)/sc)))
  
  return(log.den3)
  
}
# #-------------------------------------------------------------------