# Internal implementation for the msrlos package.
# Source: rk4d.fit.park.15Sep25.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

#  rk4d.fit.park(xdat)

#------------------------------------------------------------------
rk4d.fit.park = function (xdat, r = NULL, penk = NULL, penh = NULL, ydat = NULL, 
                          mul = NULL, sigl = NULL, 
                          shl = NULL, hl = NULL, mulink = identity, 
                          siglink = identity, shlink = identity, hlink = identity, 
                          num_inits = 20, muinit = NULL, 
                          siginit = NULL, shinit = NULL, 
                          hinit = NULL, show = TRUE, 
                          method = "Nelder-Mead", maxit = 1000, 
                          uph= 1.50, reltol=1.e-6, 
                          low.xi= -0.7,
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
  npsh <- length(shl) + 1
  nph <- length(hl) + 1
  z$trans <- FALSE
  
  # mu_names <- if (is.null(mul)) 
  #   "mu"
  # else c("mu0", paste0("mu", seq_len(npmu - 1)))
  # sigma_names <- if (is.null(sigl)) 
  #   "sigma"
  # else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  # xi_names <- if (is.null(shl)) 
  #   "xi"
  # else c("xi0", paste0("xi", seq_len(npsh - 1)))
  # h_names <- if (is.null(hl)) 
  #   "h"
  # else c("h0", paste0("h", seq_len(nph - 1)))
  # 
  # kappar <- lmomco::parkap(lmomco::lmoms(xdat[, 1]), 
  #                          snap.tau4=T)$para           #park
  
  #---------------------------------------------------------
  klmom= lmomco::lmoms(xdat[, 1])
  kappar <- lmomco::parkap(klmom, snap.tau4=TRUE, nudge.tau4=1.e-3)$para
  if(is.na(kappar[1])){
    kappar <- lmomco::pargev(klmom)$para
    kappar[4]= -0.01
  }

  #-------------------- park modified the above-----------------
  
  init_list <- list(kappar)
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(muinit)) 
      muinit <- kappar[1]
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(kappar[1], rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(siginit)) 
      siginit <- kappar[2]
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(kappar[2], rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(shinit)) 
      shinit <- kappar[3]
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(kappar[3], rep(0, length(shl)))
  }
  if (is.null(hl)) {
    hmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(hinit)) 
      hinit <- kappar[4]
  }
  else {
    z$trans <- TRUE
    hmat <- cbind(rep(1, dim(xdat)[1]), ydat[, hl])
    if (is.null(hinit)) 
      hinit <- c(kappar[4], rep(0, length(hl)))
  }
  z$model <- list(mul, sigl, shl, hl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink, hlink)))
  z1 <- as.matrix(xdat[, 1], ncol = 1)
  zr <- as.matrix(xdat[, r], ncol = 1)
  init <- c(muinit, siginit, shinit, hinit)
  
  #  names(init) <- c(mu_names, sigma_names, xi_names, h_names)
  init_list <- list(init)
  init_list[[2]]= c(lmomco::pargev(klmom)$para, -0.1) 
  
  for (i in 3:num_inits) {
    new_init <- init + c(stats::rnorm(npmu, mean = 0, sd = 1), 
                         abs(stats::rnorm(npsc, mean = 0, sd = 1)),
                         stats::rnorm(npsh, mean = 0, sd = 0.1), 
                         stats::rnorm(nph, mean = 0,  sd = 0.1))
    init_list[[i]] <- new_init
  }
 #----------------------------------------------------  
  k4d.lik.park <- function(a, uph=uph, low.xi=low.xi) {
  
    # if (z$trans == FALSE) { 
      mu=a[1]
      sc=a[2]
      xi=a[3]
      h=a[4]
      penalty=0
      
    # }else{
    #   mu <- drop(mumat %*% (a[1:npmu]))
    #   sc <- drop(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    #   xi <- drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    #   h <- drop(hmat %*% (a[seq(npmu + npsc + nph + 1, length = nph)]))
    # }
    
    #    cat("a=, k4d.lik=",a,"\n")
    
    if( any(is.na(a)) | any(is.null(a))) return(10^6)  # park
    if( xi <= low.xi | xi > 0.9999) return(10^6) # park
    if( abs(xi) <= 0.001) xi = sign(xi)*0.001
    if( abs(h) <= 0.001) h = sign(h)*0.001
        if( h <=  -1.5 | h > uph) return(10^6) # park
    if( sc > 5*kappar[2]) return(10^6) # park
    if( h < 0 & xi*h <= -1) return(10^6) # park
    
    y <- (xdat - mu)/sc
    y <- 1 - xi * y
    
    if ( min(y, na.rm=T) <= 0) return(10^6)
#    if (any(y <= 0) ) return(10^6)
    
    f <- 1 - h * y^(1/xi)
    if (any(y <= 0) || any(sc <= 0) || any(f <= 0) || any(f^(1/h) > 
                                                          1)) 
      return(10^6)
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
    sum(log(sc)) + sum(log(1 - h * y^(1/xi)) * ((h - 1)/h)) + 
      sum(log(y) * (1 - 1/xi)) - penalty
  }
  
  # a=c(  29.26484795,  8.97460686, -0.43531501, -1.10000000 )
  #------------------------------------------------------
  rk4d.lik.park <- function(a, uph=uph, low.xi=low.xi) {

   # if( z$trans == FALSE){
      mu=a[1]
      sc=a[2]
      xi=a[3]
      h=a[4]
      penalty=0
      
      if(r >=2 & h >= 1/(r-1)) return(10^6) # park
      
    # }else{
    #   mu <- drop(mumat %*% (a[1:npmu]))
    #   sc <- drop(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    #   xi <- drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    #   h <- drop(hmat %*% (a[seq(npmu + npsc + nph + 1, length = nph)]))
    # }
    
    if( any(is.na(a)) | any(is.null(a))) return(10^6)  # park
    if( xi <=  low.xi | xi > 0.999) return(10^6) # park
    if( abs(xi) <= 0.001) xi = sign(xi)*0.001
    if( abs(h) <= 0.001) h = sign(h)*0.001
    if( sc <= 0.001) sc = 0.001
    #    if( h <  -.999 | h > .999) return(10^6) # park
    if( h <=  -1.5 | h > uph) return(10^6) # park
    if( sc > 5*kappar[2]) return(10^6) # park
    if( h < 0 & xi*h <= -1) return(10^6) # park
    
    ri <- (r - seq(1:(r)))
    cr <- (1 - ri * h[1])
    if (any(sc <= 0) | any(cr <= 0)) 
      return(10^6)
    y <- 1 - xi * (xdat - mu)/sc
    
     if ( min(y, na.rm=T) <= 0) return(10^6)
    
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
 #----------------------------------------------------------------------------- 
  
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      if (z$trans == F) {
        upsc = kappar[2]*2.5
        tryCatch( try(stats::optim(init, k4d.lik.park, hessian = F
                                  ,method=method,
                                   control = list(maxit = maxit, reltol=reltol),
                                   uph=uph, low.xi=low.xi
                                  ) 
        )
        )
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, k4d.lik.park, 
                                       control = list(trace = 0),
                                       uph=uph, low.xi=low.xi ))
      }
    }
    else {
      if (z$trans == F) {
        uph = 1/(r-1) -0.01
        upsc = kappar[2]*2.5
        tryCatch( try( stats::optim(init, rk4d.lik.park, hessian = F,
                                    method=method,
                                    control = list(maxit = maxit, reltol=reltol),
                                    uph=uph, low.xi=low.xi
                                    )
        ))
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, rk4d.lik.park, 
                                       control = list(trace = 0),
                                       uph=uph, low.xi=low.xi))
      }
    }
  })
  
  optim_value <- data.frame(num = 1:length(optim_results), 
                            nllh = sapply(optim_results, function(res) {
                              if (z$trans) 
                                min(res$values)
                              else res$value
                            })
                            # , grad = sapply(optim_results, function(res) {
                            #   sum(abs(if (r == 1) numDeriv::grad(k4d.lik.park, res$par)
                            #           else numDeriv::grad(rk4d.lik.park,
                            #                               res$par)))
                            #, grad=0
                            )
  
  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  
  # optim_value <- optim_value[order(optim_value$grad, 
  #                                  optim_value$nllh), ]
  
  optim_value <- optim_value[order(optim_value$nllh), ]
  
  # minid = which.min(optim_value$nllh)
  # best_result <- optim_results[[minid]]
  
  best_result <- optim_results[[optim_value$num[1]]]
  
  if( z$trans == FALSE){
    mu=best_result$par[1]
    sc=best_result$par[2]
    xi=best_result$par[3]
    h=best_result$par[4]
  }else{
    mu <- drop(mumat %*% (best_result$par[1:npmu]))
    sc <- drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
    xi <- drop(shmat %*% (best_result$par[seq(npmu + npsc + 1, 
                                              length = npsh)]))
    h <- drop(hmat %*% (best_result$par[seq(npmu + npsc + nph + 
                                              1, length = nph)]))
  }
  
  z$r <- r
  z$conv <- best_result$convergence
  z$data <- xdat
  z$mle <- best_result$par
  
  if(is.null(penk) & is.null(penh)){
    z$nllh <- best_result$value
  }else{
    z$nllh.pen <- best_result$value
    penk=NULL; penh=NULL
    if(z$r >=2) z$nllh = rk4d.lik.park(z$mle, uph=uph, low.xi=low.xi)
    if(z$r ==1) z$nllh =  k4d.lik.park(z$mle, uph=uph, low.xi=low.xi)
  }

  #  z$cov <- solve(best_result$hessian)
  #  z$se <- sqrt(diag(z$cov))
  #  z$vals <- cbind(mu, sc, xi, h)
  # if (show) {
  #   if (z$trans) 
  #     print(z[c(2, 3)])
  #   if (z$conv != 0) 
  #     print(z[c(4, 6, 8, 10)])
  # }
  class(z) <- "rk4d.fit"

  if( is.null(best_result) ) {
    z = rgev.fit.park(xdat, r=r, num_inits=30, show=F)
    z$mle[4] = -0.0001
    z$conv = 10
    penk=NULL; penh=NULL
    if(r >=2) z$nllh = rk4d.lik.park(z$mle, uph=uph, low.xi=low.xi)
    if(r ==1) z$nllh =  k4d.lik.park(z$mle, uph=uph, low.xi=low.xi)
#    z$nllh = rk4d.lik.cvnll.park(z$mle, r=r, xdat)
  }
  invisible(z)
  return(z)
}
#----------------------------------------
#----------------------------------------------------------------
rk4d.lik.cvnll.park <- function(a, r=NULL, xdat=NULL, penk=NULL, 
                                penh=NULL, penalty=0) {
  mu <- a[1] #drop(mumat %*% (a[1:npmu]))
  sc <- a[2] #drop(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
  xi <- a[3] #drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
  h <- a[4] #drop(hmat %*% (a[seq(npmu + npsc + nph + 1, length = nph)]))
  
  ri <- (r - seq(1:(r)))
  cr <- (1 - ri * h[1])
  if (any(sc <= 0) | any(cr < 0)) 
    return(10^6)
  
  if( xi <  -0.4999 | xi > .999) return(10^6) # park
  if( abs(xi) <= 0.001) xi = sign(xi)*0.001
  if( abs(h) <= 0.001) h = sign(h)*0.001
  if( h <  -1.2 | h > 1.2) return(10^6) # park
  
  zr <- as.matrix(xdat[, r], ncol = 1)  # park
  
  y <- 1 - xi * (xdat - mu)/sc
  f <- 1 - h * (1 - xi * (zr - mu)/sc)^(1/xi)
  if (max(h) > min(1/(1 - xi * (zr - mu)/sc)^(1/xi), na.rm = TRUE)) 
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
    if (penh == "CDa") {
      if (h[1] >= 1.2 || h[1] <= -1.2) 
        return(10^6)
      if (h[1] >= 0) {
        p_h <- 1
      }
      else if (h[1] > -1.2 && h[1] < 0) {
        p_h <- exp(-((1/(1.5 + h[1])) - 0.67))
      }
      else if (h[1] <= -1.2 || h[1] >= 1.2) {
        p_h <- 0
      }
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
    else if (penh == "MSo") {
      if (h[1] <= -0.5 || h[1] >= 0.5) 
        return(10^6)
      Bef <- function(x) {
        ((0.5 + x)^(6 - 1)) * ((0.5 - x)^(9 - 1))
      }
      Be <- stats::integrate(Bef, lower = -0.5, upper = 0.5)[1]$value
      p_h <- ((0.5 + h[1])^(6 - 1)) * ((0.5 - h[1])^(9 - 
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
  sum((r * h - 1)/h * log(f) + y, na.rm = TRUE) - penalty
}
