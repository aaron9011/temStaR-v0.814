library(functional)
library(nloptr)
library(pracma)
library(spatstat)
library(Matrix)

chf_subTS <- function(u, param){
  a <- param[1]
  th <- param[2]
  if ((length(param))>=3){
    dt <- param[3]
  } else {
    dt <- 1
  }  
  y <- exp( dt*(-2*th^(1-a/2)/a*((th-1i*u)^(a/2)-th^(a/2) ) ) );  
  return( y )
}

dsubTS <- function(x, subtsparam){
  alpha <- subtsparam[1]
  theta <- subtsparam[2]
  if ((length(subtsparam))>=3){
    t <- subtsparam[3]
  } else {
    t <- 1
  }
  param <- c(alpha, t*theta)
  pdf <- (1/t)*pdf_FFT( x/t, param,  functional::Curry(chf_subTS))
  return( pdf )
}

psubTS <- function(x, subtsparam){
  alpha <- subtsparam[1]
  theta <- subtsparam[2]
  if ((length(subtsparam))>=3){
    t <- subtsparam[3]
  } else {
    t <- 1
  }
  param <- c(alpha, t*theta)
  cdf_subts <- cdf_FFT_GilPelaez( x/t, param,  functional::Curry(chf_subTS))
  return( cdf_subts )
}

#' @export
qsubTS <- function(u, subtsparam){
  ipcts(u, subtsparam)
}

#' @export
rsubTS <- function(n, subtsparam){
  u <- rand(1,n)
  r <- ipsubTS(u, subtsparam)
  return( c(r) )
}

ipsubTS <-function( u, subtsparam, maxt = 10, du = 0.01){
  arg <- seq(from = du, to = maxt, by = du)
  c <- psubTS( arg, subtsparam )
  cc <- cleanupCdf(c, arg)
  x <- cc$x[cc$x>0]
  y <- cc$y[cc$y>0]
  x = c(0,x)
  y = c(0,y)
  
  ic <- pracma::pchip(y, x, u)
  return( ic )
}