# Internal implementation for the msrlos package.
# Source: rglo.fit.park.15Sep25.R
# Scientific calculations are retained from the supplied research code;
# package-level cleanup is limited to namespace handling and side-effect removal.

#' @name rglo.fit.park
#' @aliases rglo.fit.park
#' @title Generalized Logistic Distribution for $r$ Largest Order Statistics
#' @description
#' Maximum-likelihood fitting for the order statistic model,including generalized linear modelling of each parameter.
#' @param xdat A numeric matrix of data to be fitted. Each row should be a vector of decreasing order, containing the largest order statistics for each year (or time period). The first column therefore contains annual (or period) maxima. Only the first \code{r} columns are used for the fitted model. By default, all columns are used.If one year (or time period) contains fewer order statistics than another, missing values can be appended to the end of the corresponding row.#'
#' @param r The largest \code{r} order statistics are used for the fitted model.
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or \code{NULL} (the default) for stationary fitting). The number of rows should be the same as the number of rows of \code{xdat}.
#' @param mul,sigl,shl Numeric vectors of integers, giving the columns of \code{ydat} that contain covariates for generalized linear modelling of the location, scale and shape parameters repectively (or \code{NULL} (the default) if the corresponding parameter is stationary).
#' @param mulink,siglink,shlink Inverse link functions for generalized linear modelling of the location, scale and shape parameters repectively.
#' @param num_inits Specifies the number of initial parameter sets to be generated for the optimization process.
#' @param muinit,siginit,shinit numeric of length equal to total number of parameters used to model the location, scale or shape parameter(s), resp.  See Details section for default (NULL) initial values.
#' @param show Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{\link{optim}} for details).
#' @param maxit The maximum number of iterations.
#' @param \dots Other control parameters for the optimization. These are passed to components of the \code{control} argument of \code{optim}.
#'
#' @return {
#'  A list containing the following components. A subset of these components are printed after the fit. If \code{show} is \code{TRUE}, then assuming that successful convergence is indicated, the components \code{nllh}, \code{mle} and \code{se} are always printed.
#'
#'  \item{trans}{An logical indicator for a non-stationary fit.}\item{model}{A list with components \code{mul}, \code{sigl} and \code{shl}.} \item{link}{A character vector giving inverse link functions.} \item{conv}{The convergence code, taken from the list returned by \code{\link{optim}}. A zero indicates successful convergence.} \item{nllh}{The negative logarithm of the likelihood evaluated at the maximum likelihood estimates.}
#'  \item{model}{A list with components \code{mul}, \code{sigl}, \code{shl} and \code{hl}.}
#'  \item{link}{A character vector giving inverse link functions.}
#'  \item{conv}{The convergence code, taken from the list returned by \code{\link{optim}}. A zero indicates successful convergence.}
#'  \item{nllh}{The negative logarithm of the likelihood evaluated at the maximum likelihood estimates.}
#'  \item{data}{The data that has been fitted. For non-stationary models, the data is standardized.}
#'  \item{mle}{A vector containing the maximum likelihood estimates.}
#'  \item{cov}{The covariance matrix.}
#'  \item{se}{A vector containing the standard errors.}
#'  \item{vals}{A matrix with three columns containing the maximum likelihood estimates of the location, scale and shape parameters at each data point.}
#'  \item{r}{The number of order statistics used.}
#' }
#' @seealso \code{\link{optim}}
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rglo.fit.park(bangkok, num_inits=20, show=T)
#' }
#' 
#' #----------------------------------------------------------------
rglo.fit.park <- function(xdat, r = dim(xdat)[2], ydat = NULL, 
                          mul = NULL, sigl = NULL, shl = NULL,
                     mulink = identity, siglink = identity, 
                     shlink = identity, num_inits = 20,
                     muinit = NULL, siginit = NULL, 
                     shinit = NULL, show = FALSE,
                     maxit = 1000, lower=c(-Inf, 0, -1.0), upper=c(Inf, Inf, 1.0),
                     method = "Nelder-Mead",
                     low.xi= -0.7, ...){
  
  # method = "L-BFGS-B",
  # park changed the method to "L-BFGS-B", and new lower and upper
  z <- list()

  # Determine the number of parameters for each component (mu, sigma, xi, h)
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1

  z$trans <- FALSE

  # # Generate parameter names based on the length of each list
  # mu_names <- if (is.null(mul)) "mu" else c("mu", paste0("mu", seq_len(npmu - 1)))
  # sigma_names <- if (is.null(sigl)) "sigma" else c("sigma", paste0("sigma", seq_len(npsc - 1)))
  # xi_names <- if (is.null(shl)) "xi" else c("xi", paste0("xi", seq_len(npsh - 1)))

  # Set initial values based on L-moments of the data and user-specified predictors
  glopar <- lmomco::parglo(lmomco::lmoms(xdat[, 1]))$para

  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(glopar)

  # is.glo(para)
  # Set up the mu matrix and initial values if each component (mu, sigma, xi, h) is provided

  if(is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( muinit)) muinit <- glopar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if( is.null( muinit)) muinit <- c(glopar[1], rep(0, length(mul)))
  }
  if(is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( siginit)) siginit <- glopar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if( is.null( siginit)) siginit <- c(glopar[2], rep(0, length(sigl)))
  }
  if(is.null(shl)) {
    shmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( shinit)) shinit <- glopar[3]
  }  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
    if( is.null( shinit)) shinit <- c(glopar[3], rep(0, length(shl)))
  }

  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))

  z1 <- as.matrix(xdat[, 1],ncol=1)
  zr <- as.matrix(xdat[, r],ncol=1)

  init <- c(muinit, siginit, shinit)
#  names(init) <- c(mu_names, sigma_names, xi_names)

  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(init)

  for (i in 2:num_inits) {
    # Apply random noise to each component of init to create diversity in starting points
    new_init <- init + c(
      stats::rnorm(npmu, mean = 0, sd = 1),  # Random noise for mu parameters
      abs(stats::rnorm(npsc, mean = 0, sd = 1)),  # Random noise for sigma parameters
      stats::rnorm(npsh, mean = 0, sd = 0.5) # Random noise for xi parameters
    )
    init_list[[i]] <- new_init
  }

  # Define the log-likelihood function for the case when r = 1

#------------------------------------------------------------------------------------
  glo.lik.park <- function(a, low.xi=low.xi) {

    mu <- a[1] #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2] #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- a[3] #shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    nsam =nrow(xdat)
    if( xi <= low.xi | xi >= 1) return(10^6) # park

    y <- (xdat[,1] - mu)/sc  #park changed xdat
    y <- 1 - xi * y

    f <- 1 + y^(1/xi)
    
    F_ <- 1-(1-1/f)

    if (any(y <= 0,na.rm=T) || any(sc <= 0,na.rm=T) 
        || any(f <=0,na.rm=T) || any(F_ >1,na.rm=T))
      return(10^6)
    
    nllh = nsam*(log(sc)) + sum(log(f) * 2) + sum(log(y) * (1 - 1/xi))
    nllh[1] 

  }
#----------------------------------------------------------------------------------------
  # Define the log-likelihood function for the case when r > 1
  rglo.lik.park <- function(a, low.xi=low.xi) {

    mu <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
    sc <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    xi <- a[3] #shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))

    if( xi <= low.xi | xi >= 1) return(10^6) #park
    
    ri  <- (r-seq(1:(r))) # r-i
    cr  <- (1+ri)    # c_r
    nsam =nrow(xdat)

    # constraints 1 #

    if (any(sc <= 0,na.rm=T) | any(cr < 0,na.rm=T) ) return(10^6)

    y <- 1 - xi * (xdat[,1:r] - mu)/sc   #park changed xdat
    f <- 1 + (1 - xi * (zr - mu)/sc)^(1/xi)

    # constraints 2,3,4 #

    if (any(y<= 0, na.rm = TRUE)  || any(f<= 0, na.rm = TRUE) ) return(10^6)
    
    f2 = (1 - 1/xi) * log(y)

    w <-  f2 - log(cr)
    
    if(r>1){ wrs <- rowSums(w, na.rm = TRUE) 
    
    }else{
      wrs=w
    }

    nllh = sum((r + 1) * log(f), na.rm=T) + sum(wrs) + nsam*(r*log(sc))
    nllh[1]

  }

  # Apply optimization on each set of initial values and retain results
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      stats::optim(init, glo.lik.park, hessian=F, method=method, 
                   control=list(maxit=maxit),
                   low.xi=low.xi) # Park changed to Hessian= F
    } else {
      stats::optim(init, rglo.lik.park, hessian=F, method=method, 
                   control=list(maxit=maxit),
                   low.xi=low.xi) # Park changed to Hessian= F
    }
  })

  # Collect optimization results and filter out invalid results
  optim_value <- data.frame(
    num = 1:length(optim_results),
    nllh = sapply(optim_results, function(res) res$value)
    #,
     # grad = sapply(optim_results, function(res) {
     #   sum(abs(if (r == 1) numDeriv::grad(glo.lik, res$par) else numDeriv::grad(rglo.lik, res$par)))
     # })
  )


  optim_value <- optim_value[optim_value$nllh != 10^6, ]
#  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]
  
  optim_value <- optim_value[order(optim_value$nllh), ]
  
  best_result <- optim_results[[optim_value$num[1]]]

  # Extract and store the best-fit parameter values
  mu <- best_result$par[1] #drop(mumat %*% (best_result$par[1:npmu]))
  sc <- best_result$par[2] #drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
  xi <- best_result$par[3] #drop(shmat %*% (best_result$par[seq(npmu + npsc + 1, length = npsh)]))

  # Store results in the output list with dynamic parameter names
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle  <- best_result$par
#  z$cov  <- solve(best_result$hessian)  # Park changed to #
#  z$se   <- sqrt(diag(z$cov))    # Park changed to #
  z$vals <- cbind(mu, sc, xi)
  z$r    <- r
  if(show) {
    if(z$trans)
      print(z[c(2, 3, 4)])
    #else print(z[4])
    if(!z$conv)
      print(z[c(4, 5, 7, 9)])
  }

  class(z) <- "glo"
  invisible(z)
}
#---------------------------------------------------------------------
# rld.fit.park(xdat)
#----------------------------------------------------
rld.fit.park = function (xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, 
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
  ldpar <- lmomco::parglo(lmomco::lmoms(xdat[, 1]))$para[1:2]
  init_list <- list(ldpar)
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(muinit)) 
      muinit <- ldpar[1]
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(ldpar[1], rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if (is.null(siginit)) 
      siginit <- ldpar[2]
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(ldpar[2], rep(0, length(sigl)))
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
  ld.lik.park <- function(a) {
    mu <- a[1] #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2] #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    y <- exp(-(xdat - mu)/sc)
    f <- 1 + exp(-(xdat - mu)/sc)
    if (max(f, na.rm = T) < 0) 
      return(10^6)
    if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
        any(f < 0, na.rm = T)) 
      return(10^6)
    sum(log(sc)) - sum(log(y)) - sum(-2 * log(f))
  }
  rld.lik.park <- function(a) {
    mu <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
    sc <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    ri <- (r - seq(1:(r)))
    cr <- (1 + ri)
    if (any(sc <= 0) || any(cr < 0)) 
      return(10^6)
    y <- exp(-(xdat - mu)/sc)
    f <- 1 + exp(-(zr - mu)/sc)
    if (any(f < 0, na.rm = T)) 
      return(10^6)
    if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
        any(f < 0, na.rm = T)) 
      return(10^6)
    y <- log(sc) - log(y) - log(cr)
    y <- rowSums(y, na.rm = TRUE)
    sum((r + 1) * log(f) + y, na.rm=T)  # park
  }
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      if (z$trans == F) {
        suppressWarnings(stats::optim(init, rld.lik.park, hessian = F, method = method, 
                                      control = list(maxit = maxit, trace = 0)) )
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, rld.lik.park, 
                                       control = list(trace = 0)))
      }
    }
    else {
      if (z$trans == F) {
        suppressWarnings(stats::optim(init, rld.lik.park, hessian = F, method = method, 
                                      control = list(maxit = maxit, trace = 0)) )
      }
      else {
        suppressWarnings(Rsolnp::solnp(init, rld.lik.park, 
                                       control = list(trace = 0)))
      }
    }
  })
  optim_value <- data.frame(num = 1:length(optim_results), 
                            nllh = sapply(optim_results, function(res) {
                              if (z$trans) 
                                min(res$values)
                              else res$value
                            })) #, grad = sapply(optim_results, function(res) {
                              #sum(abs(if (r == 1) numDeriv::grad(ld.lik.park, res$par)
                              #        else numDeriv::grad(rld.lik.park, 
                              #                            res$par)))
  #})
  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  # optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]
  
  optim_value <- optim_value[order(optim_value$nllh), ]
  
  best_result <- optim_results[[optim_value$num[1]]]
  mu <- best_result$par[1] #drop(mumat %*% (best_result$par[1:npmu]))
  sc <- best_result$par[2] #drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
  z$r <- r
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle <- best_result$par
  # z$cov <- solve(best_result$hessian)
  # z$se <- sqrt(diag(z$cov))
  # z$vals <- cbind(mu, sc)
  if (show) {
    if (z$trans) 
      print(z[c(2, 3)])
    if (!z$conv) 
      print(z[c(4, 6, 8, 10)])
  }
  class(z) <- "rld.fit"
  invisible(z)
}
#<simpleError in grad.default(ld.lik.park, res$par): function returns NA at 0.0003325237652156431e-04 distance from x.>
#---------------------------------------------------------------    
rld.lik.cvnll.park <- function(a, r=NULL, xdat=NULL) {
  
  mu <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
  sc <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
  ri <- (r - seq(1:(r)))
  cr <- (1 + ri)
  if (any(sc <= 0) || any(cr < 0)) 
    return(10^6)
  y <- exp(-(xdat - mu)/sc)
  
  zr <- as.matrix(xdat[, r], ncol = 1)  # park
  
  f <- 1 + exp(-(zr - mu)/sc)
  if (any(f < 0, na.rm = T)) 
    return(10^6)
  if (any(y <= 0, na.rm = T) || any(sc <= 0, na.rm = T) || 
      any(f < 0, na.rm = T)) 
    return(10^6)
  y <- log(sc) - log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)
  sum((r + 1) * log(f) + y, na.rm=T)  # park
}
#-----------------------------------------------------------
#-------------------------------------------------------------------    
rglo.lik.cvnll.park <- function(a, r=NULL, xdat=NULL) {
  
  mu <- a[1]# mulink(drop(mumat %*% (a[1:npmu])))
  sc <- a[2]#siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
  xi <- a[3]#shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
  
  ri <- (r - seq(1:(r)))
  cr <- (1 + ri)
  if (any(sc <= 0, na.rm = T) | any(cr < 0, na.rm = T)) 
    return(10^6)
  y <- 1 - xi * (xdat - mu)/sc
  
  zr <- as.matrix(xdat[, r], ncol = 1)  # park
  
  f <- 1 + (1 - xi * (zr - mu)/sc)^(1/xi)
  if (any(y <= 0, na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) 
    return(10^6)
  y <- log(sc) + (1 - 1/xi) * log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)
  sum((r + 1) * log(f) + y, na.rm = T)
}
#-------------------------------------------------------------------  

