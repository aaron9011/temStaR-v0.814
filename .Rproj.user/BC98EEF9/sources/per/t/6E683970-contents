#' @export
#' @title rmnts
#' @description \code{rmnts} generates random vector following the n dimensional NTS distribution.
#'
#' \eqn{r = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{strPMNTS$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0.00011, 0.00048 ),
#'                  sigma = c( 0.0162, 0.0231 ),
#'                  alpha = 1.23,
#'                  theta = 3.607,
#'                  beta =  c( -0.1209,  0.0905 ),
#'                  Rho = matrix( data = c(1.0, 0.55, 0.55, 1.0), nrow = 2, ncol = 2)
#' )
#' gensim <- rmnts( strPMNTS, 100 )
#' plot(gensim)
#'
rmnts <- function( strPMNTS, numofsample, rW = NaN, rTau = NaN )
  # r = \mu + \sigma X
  # X = \beta-1+\sqrt(\tau) \gamma^T \epsilon
  # \mu = [\mu_1, \mu_2, \cdots, \mu_N]^T
  # \sigma = [\sigma_1,\sigma_2, \cdots, \sigma_N]^T
  # \beta = [\beta_1, \beta_2, \cdots, \beta_N]^T
  # \gamma = [gamma_1, gamma_2, \cdots, \gamma_N]^T
  # \gamma_n = \sqrt{1-\beta_n^2(2-\alpha)/(2\theta)}
  # \epsilon \sim N(0,\Rho)
  # \Rho Correlation metrix of N-dim standard normal distribution
{
    if( is.nan(rW) ){
    rW <- randn(numofsample, strPMNTS$ndim)
  }

  if( is.nan(rTau) ){
    rTau <- rsubTS(numofsample, c(strPMNTS$alpha, strPMNTS$theta) )
  }
  #print(rW)
  N <- strPMNTS$ndim
  beta <- matrix(data = strPMNTS$beta, nrow = 1, ncol = N)
  gamma <- sqrt(1-beta^2*(2-strPMNTS$alpha)/(2*strPMNTS$theta))
  gamma <- matrix(data = gamma, nrow = 1, ncol = N)
  reps <- t(t(chol(strPMNTS$Rho))%*%t(rW))
  rstdmnts <- (rTau-1)%*%beta+reps%*%diag(as.numeric(gamma), N, N)
  #print(rstdmnts%*%diag(as.numeric(strPMNTS$sigma), N, N))
  res_rmnts <- rstdmnts%*%diag(as.numeric(strPMNTS$sigma), N, N) + t(matrix(data=as.numeric(strPMNTS$mu), nrow = N, ncol = numofsample))

  #print(res_rmnts)
  return( res_rmnts )
}

#' @export
#' @title fitmnts
#' @description \code{fitmnts} fit parameters
#' of the n-dimensional NTS distribution.
#'
#' \eqn{r = \mu + diag(\sigma) X}
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param returndata Raw data to fit the parameters.
#' The data must be given as a matrix form.
#' Each column of the matrix contains a sequence of asset returns.
#' The number of row of the matrix is the number of assets.
#'
#'@usage
#' \code{res <- fitmnts(returndata, n)}
#' \code{res <- fitmnts(returndata, n, alphaNtheta = c(alpha, theta))}
#' \code{res <- fitmnts(returndata, n, stdflag = TRUE ) }
#' \code{res <- fitmnts(returndata, n, alphaNtheta = c(alpha, theta), stdflag = TRUE)}
#'
#' @param n Dimension of the data. That is the number of assets.
#' @param alphaNtheta If \eqn{\alpha} and \eqn{\theta} are given,
#' then put those numbers in this parameter.
#' The function fixes those parameters and fits other remaining parameters.
#' If you set \code{alphaNtheta = NULL},
#' then the function fits all parameters including \eqn{\alpha} and \eqn{\theta}.
#'
#' @param stdflag If you want only standard NTS parameter fit, set this value be TRUE.
#'
#' @return Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{res$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{res$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{res$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{res$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{res$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix of the std NTS distribution (X),
#'                     which is correlation matrix of epsilon.
#'
#' \code{res$CovMtx} : Covariance matrix of return data \eqn{r}.
#'
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library(quantmod)
#'
#' getSymbols("^GSPC", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(GSPC$GSPC.Adjusted)
#' getSymbols("^DJI", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(DJI$DJI.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#' res <- fitmnts( returndata = returndata, n=2 )
#'
#'
#' #Fix alpha and theta.
#' #Estimate alpha dna theta from DJIA and use those parameter for IBM, INTL parameter fit.
#' getSymbols("^DJI", src="yahoo", from = "2020-8-25", to = "2020-08-31")
#' prDJ <- as.numeric(DJI$DJI.Adjusted)
#' ret <- diff(log(prDJ))
#' ntsparam <-  fitnts(ret)
#' getSymbols("IBM", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(IBM$IBM.Adjusted)
#' getSymbols("INTL", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(INTL$INTL.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#' res <- fitmnts( returndata = returndata,
#'                 n = 2,
#'                 alphaNtheta = c(ntsparam["alpha"], ntsparam["theta"])  )
#'
fitmnts <- function( returndata, n, alphaNtheta = NULL, stdflag = FALSE ){
  strPMNTS <- list(
    ndim = n,
    mu = matrix(data = 0, nrow = n, ncol = 1),
    sigma = matrix(data = 1, nrow = n, ncol = 1),
    alpha = 1,
    theta = 1,
    beta = matrix(data = 0, nrow = n, ncol = 1),
    Rho = matrix( nrow = n, ncol = n),
    CovMtx = matrix( nrow = n, ncol = n)
  )

  if( stdflag ){
    stdRetData <- returndata
    strPMNTS$CovMtx <- cor(returndata)
  }
  else{
    strPMNTS$CovMtx <- cov(returndata)
    strPMNTS$mu <- matrix(data = colMeans(returndata), nrow = n, ncol = 1)
    strPMNTS$sigma <- matrix(data = sqrt(diag(strPMNTS$CovMtx)), nrow = n, ncol = 1)
    muMtx <- t(matrix(data = strPMNTS$mu, nrow = n, ncol = length(returndata[,1])))
    sigmaMtx <- t(matrix(data = strPMNTS$sigma, nrow = n, ncol = length(returndata[,1])))
    stdRetData <- (returndata-muMtx)/sigmaMtx
  }

  athb <- matrix(nrow = n, ncol = 3)
  if (is.null(alphaNtheta)){
    for( k in 1:n ){
      stdntsparam <- fitstdnts(stdRetData[,k])
      athb[k,1] <- stdntsparam[1]
      athb[k,2] <- stdntsparam[2]
      athb[k,3] <- stdntsparam[3]
    }
    alphaNtheta = c(mean(athb[,1]), mean(athb[,2]))
  }
  strPMNTS$alpha <- alphaNtheta[1]
  strPMNTS$theta <- alphaNtheta[2]

  betaVec <- matrix(nrow = n, ncol = 1)
  gammaVec <- matrix(nrow = n, ncol = 1)
  for( k in 1:n ){
    betaVec[k,1] <- fitstdntsFixAlphaThata(stdRetData[,k], alpha = strPMNTS$alpha, theta = strPMNTS$theta)
    gammaVec[k,1] <- sqrt(1-betaVec[k,1]^2*(2-strPMNTS$alpha)/(2*strPMNTS$theta))
  }
  strPMNTS$beta <- betaVec
  strPMNTS$Rho <- (cov(stdRetData)-(2-strPMNTS$alpha)/(2*strPMNTS$theta)*(betaVec%*%t(betaVec)))
  #print(strPMNTS$Rho)
  igam <- diag(as.numeric(1/gammaVec))
  strPMNTS$Rho <- igam%*%strPMNTS$Rho%*%igam
  #print(strPMNTS$Rho)
  rho <- nearPD(strPMNTS$Rho, corr=TRUE)
  strPMNTS$Rho <- matrix(data = as.numeric(rho$mat), ncol = n, nrow = n)

  return(strPMNTS)
}

fitstdntsFixAlphaThata <- function( rawdat, alpha, theta, initialparam = NaN, maxeval = 100, ksdensityflag = 1){
  if (is.nan(sum(initialparam))){
    init <- 0
  } else {
    init = initialparam
  }

  if(ksdensityflag == 0){
    Femp = ecdf(rawdat)
    x = seq(from=min(rawdat), to = max(rawdat), by = 1000)
    y = Femp(x)
  } else{
    ks <- density(rawdat)
    cdfks <- CDF(ks)
    x <- ks$x
    y <- cdfks(x)
  }

  ntsp <- nloptr::bobyqa(init,
                         functional::Curry(
                           llhfntsFixAlphaTheta,
                           alpha = alpha,
                           theta = theta,
                           x = x,
                           cemp = y,
                           dispF = 0),
                         lower = -1,
                         upper = 1,
                         control = nloptr::nl.opts(list(maxeval = maxeval)))

  beta  <-  ntsp$par*sz(alpha, theta)
  return( beta )
}

llhfntsFixAlphaTheta <- function(betaparam, alpha, theta, x, cemp, dispF = 0){
  betaparam <- betaparam*sz(alpha, theta)
  Fnts <- pnts(x, c(alpha, theta, betaparam) )
  MSE <- mean((cemp-Fnts)^2)
  #SSE = sqrt(sum((cemp-Fcts)^2))/length(x)

  if (dispF == 1){
    cat(ntsparam, sep = ',')
    cat(";")
    cat(MSE)
    cat("\n")
  }
  return( MSE )
}


#' @title setPortfolioParam
#' @description Please use \code{getPortNTSParam} instead of \code{setPortfolioParam}.
#'
#' Portfolio return with capital allocation weight is \eqn{R_p=<w,r>},
#' which is a weighted sum of of elements in the N-dimensional NTS random vector.
#' \eqn{R_p} becomes an 1-dimensional NTS random variable.
#' \code{setPortfolioParam} find the parameters of \eqn{R_p}.
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{strPMNTS$Rho} : \eqn{\Sigma} matrix of the std NTS distribution (X).
#'
#' @param w Capital allocation weight vector.
#'
#' @usage
#' \code{res <- setPortfolioParam(strPMNTS,w)}
#'
#' @return The weighted sum follows 1-dimensional NTS.
#'
#' \eqn{R_p = <w, r> = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_1(\alpha, \theta, \beta, 1)}.
#'
#' Hence we obtain
#'
#' \code{res$mu} : \eqn{\mu} mean of \eqn{R_p}.
#'
#' \code{res$sigma} : \eqn{\sigma} standard deviation of \eqn{R_p}.
#'
#' \code{res$alpha} : \eqn{\alpha} of \eqn{X}.
#'
#' \code{res$theta} : \eqn{\theta} of \eqn{X}.
#'
#' \code{res$beta} : \eqn{\beta} \eqn{X}.
#'
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 9.876552e-05, 4.747343e-04 ),
#'                  sigma = c( 0.01620588, 0.02309643 ),
#'                  alpha = 0.1888129 ,
#'                  theta = 0.523042,
#'                  beta =  c( -0.04632938,  0.04063555 ),
#'                  Rho = matrix( data = c(1.0, 0.469883,
#'                                        0.469883, 1.0),
#'                                nrow = 2, ncol = 2)
#'                  CovMtx = matrix( data = c(0.0002626304, 0.0001740779,
#'                                          0.0001740779, 0.0005334452),
#'                                  nrow = 2, ncol = 2)
#' )
#' w <- c(0.3, 0.7)
#' res <- setPortfolioParam(strPMNTS,w)
#'
setPortfolioParam <- function(strPMNTS, w){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  w <- matrix(data = w, nrow = strPMNTS$ndim, ncol = 1)
  m <- sum(w*strPMNTS$mu)
  b <- sum(w*strPMNTS$sigma*strPMNTS$beta)
  gammaMtx <- diag(as.numeric(sqrt(1-strPMNTS$beta^2*(2-strPMNTS$alpha)/(2*strPMNTS$theta))))
  sigMtx <- diag(as.numeric(strPMNTS$sigma))
  g <- as.numeric(sqrt(t(w)%*%sigMtx%*%gammaMtx%*%strPMNTS$Rho%*%gammaMtx%*%sigMtx%*%w))
  a <- as.numeric(strPMNTS$alpha)
  th <- as.numeric(strPMNTS$theta)
  param <- change_ntsparam2stdntsparam(c(a, th, b, g, m))
  return(param)
}

#' @export
#' @title getPortNTSParam
#' @description Portfolio return with capital allocation weight is \eqn{R_p=<w,r>},
#' which is a weighted sum of of elements in the N-dimensional NTS random vector.
#' \eqn{R_p} becomes an 1-dimensional NTS random variable.
#' \code{getPortNTSParam} find the parameters of \eqn{R_p}.
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param stdform If \code{stdform} is \code{FALSE}, then the return parameter has the following representation
#'
#' \eqn{R_p = <w, r> = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_1(\alpha, \theta, \beta, 1)}.
#'
#' If \code{stdform} is \code{TRUE}, then the return parameter has the following representation
#'
#' \eqn{R_p = <w, r>} follows \eqn{stdNTS_1(\alpha, \theta, \beta, \gamma, \mu)}
#'
#' @usage
#' \code{res <- setPortfolioParam(strPMNTS,w)}
#' \code{res <- setPortfolioParam(strPMNTS,w, FALSE)}
#'
#' @return The weighted sum follows 1-dimensional NTS.
#'
#' \eqn{R_p = <w, r> = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_1(\alpha, \theta, \beta, 1)}.
#'
#' Hence we obtain
#'
#' \code{res$mu} : \eqn{\mu} mean of \eqn{R_p}.
#'
#' \code{res$sigma} : \eqn{\sigma} standard deviation of \eqn{R_p}.
#'
#' \code{res$alpha} : \eqn{\alpha} of \eqn{X}.
#'
#' \code{res$theta} : \eqn{\theta} of \eqn{X}.
#'
#' \code{res$beta} : \eqn{\beta} of \eqn{X}.
#'
#' @references
#' Proposition 2.1 of
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 9.876552e-05, 4.747343e-04 ),
#'                  sigma = c( 0.01620588, 0.02309643 ),
#'                  alpha = 0.1888129 ,
#'                  theta = 0.523042,
#'                  beta =  c( -0.04632938,  0.04063555 ),
#'                  Rho = matrix( data = c(1.0, 0.469883,
#'                                        0.469883, 1.0),
#'                                nrow = 2, ncol = 2)
#'                  CovMtx = matrix( data = c(0.0002626304, 0.0001740779,
#'                                          0.0001740779, 0.0005334452),
#'                                  nrow = 2, ncol = 2)
#'                  )
#' w <- c(0.3, 0.7)
#' res <- getPortNTSParam(strPMNTS,w)
#'
getPortNTSParam <- function(strPMNTS, w, stdform = TRUE){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  w <- matrix(data = w, nrow = strPMNTS$ndim, ncol = 1)
  sigmaBar <- as.numeric(sqrt(t(w)%*%strPMNTS$CovMtx%*%w))
  mbar <- sum(w*strPMNTS$mu)
  betbar <- sum(w*strPMNTS$sigma*strPMNTS$beta)/sigmaBar
  stdntsparam <- c(as.numeric(strPMNTS$alpha), as.numeric(strPMNTS$theta), betbar)
  names(stdntsparam) <- c("alpha", "theta", "beta")
  param <- list(stdparam = stdntsparam, mu = mbar, sig = sigmaBar)
  if (stdform == FALSE){
    param <- change_stdntsparam2ntsparam(param$stdparam, mbar, sigmaBar)
    param <- param[1:5]
  }
  return(param)
}

#' @export
#' @title portfolioVaRmnts
#' @description
#' Calculate portfolio value at risk on the NTS market model
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param eta significanlt level
#' @return portfolio value at risk on the NTS market model
#'
portfolioVaRmnts <- function(strPMNTS, w, eta){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  ntsparam <- getPortNTSParam(strPMNTS, w)
  VaR <- (- ntsparam$mu - ntsparam$sig * qnts(eta, ntsparam$stdparam))
  return(VaR)
}

#' @export
#' @title portfolioCVaRmnts
#' @description
#' Calculate portfolio conditional value at risk (expected shortfall) on the NTS market model
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param eta significanlt level
#' @return portfolio value at risk on the NTS market model
#'
portfolioCVaRmnts <- function(strPMNTS, w, eta){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  ntsparam <- getPortNTSParam(strPMNTS, w)
  CVaR <- (- ntsparam$mu + ntsparam$sig * cvarnts(eta, ntsparam$stdparam))
  return(CVaR)
}


#' @export
#' @title pmnts
#' @description \code{pmnts} calculates the cdf values of the multivariate NTS distribution:
#' \eqn{F(x_1, \cdots, x_n)=P(x_n<R_1, \cdots, x_n<R_n)}.
#' The multivariate NTS random vector \eqn{R = (R_1, \cdots, R_n)} is defined
#'
#' \eqn{R = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param x array of the \eqn{(x_1, \cdots, x_n)}
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{st$ndim} : dimension
#'
#' \code{st$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{st$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{st$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{st$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{st$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{st$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library(mvtnorm)
#' strPMNTS <- list(ndim = 2,
#'               mu = c( 0.5, -1.5 ),
#'               sigma = c( 2, 3 ),
#'               alpha = 0.1,
#'               theta = 3,
#'               beta =  c( 0.1, -0.3 ),
#'               Rho = matrix( data = c(1.0, 0.75, 0.75, 1.0),
#'                             nrow = 2, ncol = 2)
#' )
#' pmnts(c(0.6, -1.0), st = strPMNTS)
#'
#'
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0, 0, 0 ),
#'                  sigma = c( 1, 1, 1 ),
#'                  alpha = 0.1,
#'                  theta = 3,
#'                  beta =  c( 0.1, -0.3, 0 ),
#'                  Rho = matrix(
#'                      data = c(1.0, 0.75, 0.1, 0.75, 1.0, 0.2, 0.1, 0.2, 1.0),
#'                      nrow = 3, ncol = 3)
#' )
#' pmnts(c(0,0,0), st = strPMNTS)
#' dmnts(c(0,0,0), st = strPMNTS)
#'
pmnts <- function( x, st, subTS = NULL ){
  xadj <- (x-st$beta)/st$sigma
  return(pMultiStdNTS(xadj, st, subTS))
}

#' @export
#' @title dmnts
#' @description \code{dmnts} calculates the density of the multivariate NTS distribution:
#' \eqn{f(x_1, \cdots, x_n)=\frac{d^n}{dx_1\cdots dx_n}P(x_n<R_1, \cdots, x_n<R_n)}.
#' The multivariate NTS random vector \eqn{R = (R_1, \cdots, R_n)} is defined
#'
#' \eqn{R = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param x array of the \eqn{(x_1, \cdots, x_n)}
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{st$ndim} : dimension
#'
#' \code{st$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{st$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{st$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{st$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{st$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{st$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library(mvtnorm)
#' strPMNTS <- list(ndim = 2,
#'               mu = c( 0.5, -1.5 ),
#'               sigma = c( 2, 3 ),
#'               alpha = 0.1,
#'               theta = 3,
#'               beta =  c( 0.1, -0.3 ),
#'               Rho = matrix( data = c(1.0, 0.75, 0.75, 1.0),
#'                             nrow = 2, ncol = 2)
#' )
#' dmnts(c(0.6, -1.0), st = strPMNTS)
#'
#'
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0, 0, 0 ),
#'                  sigma = c( 1, 1, 1 ),
#'                  alpha = 0.1,
#'                  theta = 3,
#'                  beta =  c( 0.1, -0.3, 0 ),
#'                  Rho = matrix(
#'                      data = c(1.0, 0.75, 0.1, 0.75, 1.0, 0.2, 0.1, 0.2, 1.0),
#'                      nrow = 3, ncol = 3)
#' )
#' pmnts(c(0,0,0), st = strPMNTS)
#' dmnts(c(0,0,0), st = strPMNTS)
#'
dmnts <- function( x, st, subTS = NULL ){
  xadj <- (x-st$beta)/st$sigma
  return(dMultiStdNTS(xadj, st, subTS)/prod(st$sigma))
}

#' @export
#' @title copulaStdNTS
#' @description \code{copulaStdNTS} calculates the stdNTS copula values
#'
copulaStdNTS <- function(u, st, subTS = NULL){
  x <- matrix(nrow = length(u), ncol = 1)
  for (j in 1:length(u)){
    x[j] <- temStaR::qnts(u[j], c(st$alpha, st$theta, st$beta[j]))
  }
  pMultiStdNTS( as.numeric(x), st, subTS )
}

#' @export
#' @title dcopulaStdNTS
#' @description \code{dcopulaStdNTS} calculates
#' density of the stdNTS copula.
#'
dcopulaStdNTS <- function(u, st, subTS = NULL){
  x <- matrix(nrow = length(u), ncol = 1)
  y <- 1
  for (j in 1:length(u)){
    x[j] <- temStaR::qnts(u[j], c(st$alpha, st$theta, st$beta[j]))
    y <- y*temStaR::dnts(x[j], c(st$alpha, st$theta, st$beta[j]))
  }
  dMultiStdNTS( as.numeric(x), st, subTS )/y
}

#' @export
#' @title importantSamplining
#' @description \code{importantSamplining} do the important sampling for the TS Subordinator.
#'
importantSamplining <- function( alpha, theta ){
  u  <- c(
    seq(from=0, to = 0.009, length.out = 10),
    seq(from=0.01, to = 0.09, length.out = 9),
    seq(from=0.1, to = 0.9, length.out = 9),
    0.95, 0.99, 0.999, 0.9999, 0.99999, 0.999999
  )
  ti <- temStaR::qsubTS(u, c(alpha, theta))
  ti = c(ti, max(ti)*2)
  subtsi <- temStaR::dsubTS(ti, c(alpha,theta))
  subTS <- list( ti = ti, subtsi = subtsi)
}

dMultiNorm_Subord <- function( tVec, x, alpha, theta, beta, rhoMtx ){
  gamma <- as.numeric(sqrt(1-(2-alpha)/(2*theta)*beta^2))
  re <- matrix(nrow = length(tVec), ncol = 1)
  for (i in 1:length(tVec) ){
    t <- tVec[i]
    #mu <- beta*(t-1)
    #Sig <- matrix(nrow = 2, ncol = 2)
    #Sig[1,1] <- gamma[1]^2*t
    #Sig[2,2] <- gamma[2]^2*t
    #Sig[1,2] <- gamma[1]*gamma[2]*rhoMtx[1,2]*t
    #Sig[2,1] <- Sig[1,2]
    mu <- beta*(t-1)
    Sig <- t*(cbind(gamma)%*%rbind(gamma))*rhoMtx

    re[i] <- dmvnorm(x, mean=mu, sigma = Sig )
  }
  return(re)
}

func_indegrand <- function(t, x, st, ti, subtsi){
  fe <- dMultiNorm_Subord(t,
                          x = x,
                          alpha = st$alpha,
                          theta = st$theta,
                          beta = st$beta,
                          rhoMtx = st$Rho)
  ft <- pracma::pchip(ti, subtsi, t)
  #ft <- temStaR::dsubTS(t, c( st$alpha,  st$theta))
  return( fe*ft )
}

dMultiStdNTS <- function( x, st, subTS = NULL ){
  if ( is.null(subTS) ){
    subTS <- importantSamplining(st$alpha, st$theta)
  }
  ti <- subTS$ti
  subtsi <- subTS$subtsi

  d <- pracma::quad(
    functional::Curry(func_indegrand,
                      x = x,
                      st = st,
                      ti = ti,
                      subtsi = subtsi),
    xa = 0,
    xb = max(ti)
  )
  return(d)
}

pMultiNorm_Subord <- function( tVec, x, alpha, theta, beta, rhoMtx ){
  gamma <- as.numeric(sqrt(1-(2-alpha)/(2*theta)*beta^2))
  re <- matrix(nrow = length(tVec), ncol = 1)
  for (i in 1:length(tVec) ){
    t <- tVec[i]
    if (t == 0){
      re[i] <- 0
    }
    else{
      adjx <- (x-beta*(t-1))/(gamma*sqrt(t))
      re[i] <- pmvnorm(lower = c(-Inf, -Inf),
                       upper = adjx, mean = rep(0, length(x)), sigma = rhoMtx )
    }
  }
  return(re)
}

func_indegrand_cdf <- function(t, x, st, ti, subtsi){
  Gt <- pMultiNorm_Subord(t,
                          x = x,
                          alpha = st$alpha,
                          theta = st$theta,
                          beta = st$beta,
                          rhoMtx = st$Rho)
  ft <- pracma::pchip(ti, subtsi, t)
  #ft <- temStaR::dsubTS(t, c( st$alpha,  st$theta))
  return( Gt*ft )
}

pMultiStdNTS <- function( x, st, subTS = NULL ){
  if ( is.null(subTS) ){
    subTS <- importantSamplining(st$alpha, st$theta)
  }
  ti <- subTS$ti
  subtsi <- subTS$subtsi

  p <- pracma::quad(
    functional::Curry(func_indegrand_cdf,
                      x = x,
                      st = st,
                      ti = ti,
                      subtsi = subtsi),
    xa = 0,
    xb = max(ti)
  )
  p[p>1]=1
  return(p)
}
