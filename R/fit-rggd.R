# Internal implementation for the msrlos package.
# Source: rggd.fit.park_15Sep25.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

#-----------------------------------------------------------------------
rggd.fit.park =function (xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, 
                         hl = NULL, mulink = identity, 
                         siglink = identity, hlink = identity, 
                         num_inits = 20, muinit = NULL, 
                         siginit = NULL, hinit = NULL, 
                         show = TRUE, method = "Nelder-Mead", maxit = 1000, 
                         reltol=1e-6, ...) 
{
  z <- list()
  if (is.null(r)) {
    if (is.vector(xdat)) {
      xdat <- matrix(xdat, ncol = 1)
      r <- 1
    }
    else {
      r <- dim(xdat)[2]
    }
  }
  else {
    if (r == 1) {
      if (is.vector(xdat)) {
        xdat <- matrix(xdat, ncol = 1)
      }
      else {
        xdat <- matrix(xdat[, 1:r], ncol = 1)
      }
    }
    else {
      xdat <- as.matrix(xdat[, 1:r], ncol = r)
    }
  }

  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  nph <- length(hl) + 1
  z$trans <- FALSE
  mu_names <- if (is.null(mul)) 
    "mu"
  else c("mu", paste0("mu0", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) 
    "sigma"
  else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  h_names <- if (is.null(hl)) 
    "h"
  else c("h", paste0("h0", seq_len(nph - 1)))
  
  xpar <- lmomco::parkap(lmomco::lmoms(xdat[, 1]), snap.tau4=TRUE,
                         nudge.tau4=1.e-4)$para
  ggdpar = c(xpar[1:2],xpar[4])   # park
  if(is.na(ggdpar[1])) {
    work= lmomco::pargev(lmomco::lmoms(xdat[, 1]))$para
    ggdpar =c(work[1:2], -0.1)
  }
  
  init_list <- list()
  init_list[[1]] = ggdpar   # park

  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(muinit)) 
      muinit <- ggdpar[1]
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(ggdpar[1], rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(siginit)) 
      siginit <- ggdpar[2]
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(ggdpar[2], rep(0, length(sigl)))
  }
  if (is.null(hl)) {
    hmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(hinit)) 
      hinit <- ggdpar[3]
  }
  else {
    z$trans <- TRUE
    hmat <- cbind(rep(1, dim(xdat)[1]), ydat[, hl])
    if (is.null(hinit)) 
      hinit <- c(ggdpar[3], rep(0, length(hl)))
  }
  z$model <- list(mul, sigl, hl)
  z$link <- deparse(substitute(c(mulink, siglink, hlink)))
  z1 <- as.matrix(xdat[, 1], ncol = 1)
  zr <- as.matrix(xdat[, r], ncol = 1)
  init <- c(muinit, siginit, hinit)
  names(init) <- c(mu_names, sigma_names, h_names)
  
  work= lmomco::pargev(lmomco::lmoms(xdat[, 1]))$para
  init_list[[2]] = c(work[1:2], init[3])
  
  for (i in 3:(num_inits)) {
    new_init <- init + c(stats::rnorm(npmu, mean = 0, sd = 1), 
                         abs(stats::rnorm(npsc, mean = 0, sd = 1)), 
                         stats::rnorm(nph, mean = 0, sd = 0.5))
    init_list[[i]] <- new_init
  }
  
#--------------------------------------------------------------  
  ggd.lik <- function(a) {
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    h <- hlink(hmat %*% (a[seq(npmu + npsc + 1, length = nph)]))
    y <- exp(-(xdat - mu)/sc)
    
    if(any(h <= -1.5)) return(10^6)  # park
    
    if (any(h > 0)) {
      f <- 1 - h * exp(-(xdat - mu)/sc)
      if (any(f < 0, na.rm = T)) 
        return(10^6)
      if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
          any(f^(1/h) < 0, na.rm = T) || any(f^(1/h) > 
                                             1, na.rm = T)) 
        return(10^6)
      sum(log(sc)) - sum(log(y)) - sum(((1 - h)/h) * log(f))
    }
    else {
      f <- 1 - h * exp(-(xdat - mu)/sc)
      if (max(f, na.rm = T) < 0) 
        return(10^6)
      if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
          any(f^(1/h) < 0, na.rm = T) || any(f^(1/h) > 
                                             1, na.rm = T)) 
        return(10^6)
      sum(log(sc)) - sum(log(y)) - sum(((1 - h)/h) * log(f))
    }
  }
#--------------------------------------------------------  
  rggd.lik <- function(a) {
    
    mu <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
    sc <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    h <-  a[3] #hlink(drop(hmat %*% (a[seq(npmu + npsc + 1, length = nph)])))
    if (r >= 2) {
      if (min(h) > (1/(r - 1))) 
        return(10^6)
      if( h[1] > 1/(r-1) ) return(10^6)
    }
    
    if(any(h <= -1.5)) return(10^6)  # park
    
    ri <- (r - seq(1:(r)))
    cr <- (1 - ri * h[1])
    if (any(sc <= 0) || any(cr < 0)) 
      return(10^6)
    if (any(h > 0)) {
      y <- exp(-(xdat - mu)/sc)
      f <- 1 - h[1] * exp(-(zr - mu)/sc)
      if (any(f < 0, na.rm = T)) 
        return(10^6)
      if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
          any(f^(1/h) < 0, na.rm = T) || any(f^(1/h) > 
                                             1, na.rm = T)) 
        return(10^6)
      if (any(min(h) > 1/exp(-(zr - mu)/sc), na.rm = TRUE)) 
        return(10^6)
      y <- log(sc) - log(y) - log(cr)
      y <- rowSums(y, na.rm = TRUE)
      sum((r * h - 1)/h * log(f) + y, na.rm = T)
    }
    else {
      y <- exp(-(xdat - mu)/sc)
      f <- 1 - h * exp(-(zr - mu)/sc)
      if (any(f < 0, na.rm = T)) 
        return(10^6)
      if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
          any(f^(1/h) < 0, na.rm = T) || any(f^(1/h) > 
                                             1, na.rm = T)) 
        return(10^6)
      if (any(min(h) > 1/exp(-(zr - mu)/sc), na.rm = TRUE)) 
        return(10^6)
      y <- log(sc) - log(y) - log(cr)
      y <- rowSums(y, na.rm = TRUE)
      sum((r * h - 1)/h * log(f) + y, na.rm = T)
    }
  }
#----------------------------------------
  
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      if (z$trans == F) {
        stats::optim(init, ggd.lik, hessian = F, method = method, 
                     control = list(maxit = maxit, trace = 0,
                                    reltol=reltol))
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, ggd.lik, 
                                       control = list(trace = 0)))
      }
    }
    else {
      if (z$trans == F) {
        stats::optim(init, rggd.lik, hessian = F, 
                     method = method, control = list(maxit = maxit, 
                                      reltol=reltol, trace = 0))
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, rggd.lik, 
                                       control = list(trace = 0)))
      }
    }
  })
  optim_value <- data.frame(num = 1:length(optim_results), 
                            nllh = sapply(optim_results, function(res) {
                              if (z$trans) 
                                min(res$values)
                              else res$value
                            }) #, grad = sapply(optim_results, function(res) {
                            #   sum(abs(if (r == 1) numDeriv::grad(ggd.lik, res$par) else numDeriv::grad(rggd.lik, 
                            #  res$par)))
                            # })
  )
  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  optim_value <- optim_value[order(optim_value$nllh), ]
  
  best_result <- optim_results[[optim_value$num[1]]]
  
  #-------------------------------------------------------
  # mu <- drop(mumat %*% (best_result$par[1:npmu]))
  # sc <- drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
  # h <- drop(hmat %*% (best_result$par[seq(npmu + npsc + 1, 
  #                                         length = nph)]))
  #  park made the above three lines as comments
  #-------------------------------------------------------
  
  z$r <- r
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle <- best_result$par
  
  if( is.null(best_result) ) {
    z = rgd.fit.park(xdat, r=r, num_inits=10, show=F)
    
    z$mle[3] = init[3]
    if( z$mle[3] > 1/(r-1)) {
        z$mle[3] = 1/(r-1) - 0.1
    }
    
    if(r==1){
      z$nllh= ggd.lik(z$mle)
    }else{ z$nllh= rggd.lik(z$mle) }
    
   # z$nllh = rggd.lik.cvnll.park(z$mle, r=r, xdat)
  }
  
  #--------------------------------------------------
  # z$cov <- solve(best_result$hessian)
  # z$se <- sqrt(diag(z$cov))
  z$vals = z$mle
  mu= z$mle[1]; sc=z$mle[2]; h = z$mle[3]; k=0
  
  if( h > 1/(r-1)) {
    cat("==== warning ====: h > 1/(r-1) at rggd.fit ","\n")
    z$mle[3] = h = 1/(r-1) - 0.1
  }
  # z$vals <- cbind(mu, sc, h)   # park made it comments
  #----------------------------------------------------
  
  if (show) {
    if (z$trans) 
      print(z[c(2, 3)])
    if (!z$conv) 
      print(z[c(4, 6, 8, 10)])
  }
  class(z) <- "rggd.fit"
     invisible(z)
  return(z)
}

#----------------------------------------------   
rgd.fit.park =function (xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, 
                        mulink = identity, siglink = identity, 
                        num_inits = 20, muinit = NULL, 
                        siginit = NULL, show = TRUE, 
                        method = "Nelder-Mead", maxit = 1000, 
                        ...) 
{
  z <- list()
  if (is.null(r)) {
    if (is.vector(xdat)) {
      xdat <- matrix(xdat, ncol = 1)
      r <- 1
    }
    else {
      r <- dim(xdat)[2]
    }
  }
  else {
    if (r == 1) {
      if (is.vector(xdat)) {
        xdat <- matrix(xdat, ncol = 1)
      }
      else {
        xdat <- matrix(xdat[, 1:r], ncol = 1)
      }
    }
    else {
      xdat <- as.matrix(xdat[, 1:r], ncol = r)
    }
  }
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  z$trans <- FALSE
  mu_names <- if (is.null(mul)) 
    "mu"
  else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) 
    "sigma"
  else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  gdpar <- lmomco::pargum(lmomco::lmoms(xdat[, 1]))$para
  init_list <- list(gdpar)
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(muinit)) 
      muinit <- gdpar[1]
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(gdpar[1], rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(siginit)) 
      siginit <- gdpar[2]
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(gdpar[2], rep(0, length(sigl)))
  }
  z$model <- list(mul, sigl)
  z$link <- deparse(substitute(c(mulink, siglink)))
  z1 <- as.matrix(xdat[, 1], ncol = 1)
  zr <- as.matrix(xdat[, r], ncol = 1)
  init <- c(muinit, siginit)
  names(init) <- c(mu_names, sigma_names)
  init_list <- list(init)
  for (i in 2:num_inits) {
    new_init <- init + c(stats::rnorm(npmu, mean = 0, sd = 1), 
                         abs(stats::rnorm(npsc, mean = 0, sd = 1)))
    init_list[[i]] <- new_init
  }
  #-----------------------------------------------
  rgd.lik.park <- function(a) {
    
    if( z$trans == FALSE){
      mu=a[1]
      sc=a[2]
    }else{
      mu <- mulink(mumat %*% (a[1:npmu]))
      sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    }
    
    if (any(sc <= 0)) 
      return(10^6)
    y <- (xdat - mu)/sc
    y <- y + log(sc)
    y <- rowSums(y, na.rm = TRUE)
    sum(exp(-(zr - mu)/sc) + y, na.rm = T)
  }
  #--------------------------------------------------
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      if (z$trans == F) {
        stats::optim(init, rgd.lik.park, hessian = F, method = method, 
                     control = list(maxit = maxit, trace = 0))
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, rgd.lik.park, 
                                       control = list(trace = 0)))
      }
    }
    else {
      if (z$trans == F) {
        stats::optim(init, rgd.lik.park, hessian = F, method = method, 
                     control = list(maxit = maxit, trace = 0))
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, rgd.lik.park, 
                                       control = list(trace = 0)))
      }
    }
  })
  optim_value <- data.frame(num = 1:length(optim_results), 
                            nllh = sapply(optim_results, function(res) {
                              if (z$trans) 
                                min(res$values)
                              else res$value
                            }) #, grad = sapply(optim_results, function(res) {
                            #   sum(abs(if (r == 1) numDeriv::grad(rgd.lik.park, res$par) 
                            #           else        numDeriv::grad(rgd.lik.park, res$par)))
                            # })
  )
  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  
  # optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]
  optim_value <- optim_value[order(optim_value$nllh), ]  # park
  
  best_result <- optim_results[[optim_value$num[1]]]
  mu <- drop(mumat %*% (best_result$par[1:npmu]))
  sc <- drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
  z$r <- r
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle <- best_result$par
  # z$cov <- solve(best_result$hessian)
  # z$se <- sqrt(diag(z$cov))
  #  z$vals <- cbind(mu, sc)
  if (show) {
    if (z$trans) 
      print(z[c(2, 3)])
    if (!z$conv) 
      print(z[c(4, 6, 8, 10)])
  }
  class(z) <- "rgd.fit"
  invisible(z)
}

#---------------------------------------------------------
#-----------------------------------------------------      
rggd.lik.cvnll.park <- function(a, r=NULL, xdat=NULL) {
  
  mu <- a[1]  #mulink(drop(mumat %*% (a[1:npmu])))
  sc <- a[2]  #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
  h <- a[3]  #hlink(drop(hmat %*% (a[seq(npmu + npsc + 1, length = nph)])))
  
  if (r >= 2) {
    if (min(h) > (1/(r - 1))) 
      return(10^6)
  }
  ri <- (r - seq(1:(r)))
  cr <- (1 - ri * h[1])
  if (any(sc <= 0) || any(cr < 0)) 
    return(10^6)
  
  zr <- as.matrix(xdat[, r], ncol = 1)  # park
  
  if (any(h > 0)) {
    y <- exp(-(xdat - mu)/sc)
    f <- 1 - h[1] * exp(-(zr - mu)/sc)
    if (any(f < 0, na.rm = T)) 
      return(10^6)
    if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
        any(f^(1/h) < 0, na.rm = T) || any(f^(1/h) > 
                                           1, na.rm = T)) 
      return(10^6)
    if (any(min(h) > 1/exp(-(zr - mu)/sc), na.rm = TRUE)) 
      return(10^6)
    y <- log(sc) - log(y) - log(cr)
    y <- rowSums(y, na.rm = TRUE)
    sum((r * h - 1)/h * log(f) + y, na.rm = T)
  }
  else {
    y <- exp(-(xdat - mu)/sc)
    f <- 1 - h * exp(-(zr - mu)/sc)
    if (any(f < 0, na.rm = T)) 
      return(10^6)
    if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
        any(f^(1/h) < 0, na.rm = T) || any(f^(1/h) > 
                                           1, na.rm = T)) 
      return(10^6)
    if (any(min(h) > 1/exp(-(zr - mu)/sc), na.rm = TRUE)) 
      return(10^6)
    y <- log(sc) - log(y) - log(cr)
    y <- rowSums(y, na.rm = TRUE)
    sum((r * h - 1)/h * log(f) + y, na.rm = T)
  }
}
#-----------------------------------------------------
#---------------------------------------------------
rgd.lik.cvnll.park <- function(a, r=NULL, xdat=NULL) {
  mu <- a[1] #mulink(mumat %*% (a[1:npmu]))
  sc <- a[2]# siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
  if (any(sc <= 0)) 
    return(10^6)
  y <- (xdat - mu)/sc
  y <- y + log(sc)
  y <- rowSums(y, na.rm = TRUE)
  
  zr <- as.matrix(xdat[, r], ncol = 1)
  
  sum(exp(-(zr - mu)/sc) + y, na.rm = T)
}
#-----------------------------------------------
# mle= lmomco::vec2par(c(100,10,-0.2),'ggd')
# quaggd(prob, mle)
# -------------------------------------------
quaggd = function(prob, mle, paracheck= TRUE){
  
  prob[prob > 0.99999] =0.99999
  prob[prob < 0.00001] =0.00001
  year = 1/(1-prob)
  zggd = NULL
  zggd$mle = mle$para
  zggd$cov = matrix(1,3,3)
  y = rggd.rl.park(zggd, year, just=TRUE)$rl
  return(y)
}

#-----------------------------------------------    
rggd.rl.park = function (z, year = c(20, 50, 100, 200), 
                         just=FALSE, show = F) 
{
  del <- matrix(ncol = length(year), nrow = 3)
  mu = z$mle[1]
  sig = z$mle[2]
  h = z$mle[3]
  f = 1 - (1/year)
  f_inv = 1/f
  del[1, ] <- 1
  del[2, ] <- -log((1 - exp(h * log(f)))/h)
  del[3, ] <- sig * (1/h + ((f^h) * log(f))/(1 - (f^h)))
  del.t <- t(del)
  if (h > 0) {
    z$rl = mu - sig * log((1 - exp(h * log(f)))/h)
  }
  else {
    z$rl = mu + sig * log((-h)/(expm1(h * log(f))))
  }
  #      z$rlse <- diag(sqrt((del.t %*% z$cov %*% del)))
  if(just !=TRUE){
    names(z$rl) <- paste0(as.character(year), "y")
  }
  #      names(z$rlse) <- paste0(as.character(year), "y")
  if (show) 
    print(z[c(12, 13)])
  invisible(z)
}
