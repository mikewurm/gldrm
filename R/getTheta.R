## Computes log(sum(exp(x))) with better precision
logSumExp <- function(x)
{
    i <- which.max(x)
    m <- x[i]
    lse <- log1p(sum(exp(x[-i]-m))) + m
    lse
}

## g function (logit transformation from appendix)
g <- function(mu, m, M) log(mu-m) - log(M-mu)

#' Control arguments for \eqn{\theta} update algorithm
#'
#' This function returns control arguments for the \eqn{\theta} update algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold for theta updates. Convergence is
#' evaluated separately for each observation. An observation has converged when
#' the difference between \eqn{b'(\theta)} and \eqn{\mu} is less than \code{epsTheta}.
#' @param maxiter Maximum number of iterations.
#' @param maxhalf Maximum number of half steps allowed per iteration if the
#' convergence criterion does not improve.
#' @param maxtheta Absolute value of theta is not allowed to exceed \code{maxtheta}.
#' @param logit Logical for whether logit transformation should be used. Use of
#' this stabilizing transformation appears to be faster in general. Default is TRUE.
#' @param logsumexp Logical argument for whether log-sum-exp trick should be used.
#' This may improve numerical stability at the expense of computational time.
#'
#' @return Object of S3 class "thetaControl", which is a list of control arguments.
#'
#' @export
theta.control <- function(eps=1e-10, maxiter=100, maxhalf=20, maxtheta=500,
                          logit=TRUE, logsumexp=FALSE)
{
    thetaControl <- as.list(environment())
    class(thetaControl) <- "thetaControl"
    thetaControl
}

#' getTheta
#' Updates theta. Vectorized but only updates observations that have not converged.
#'
#' @param spt Support of the observed response variable. (This is the set of
#' unique values observed, not the set of all possible values.)
#' @param f0 Values of the baseline distribution corresponding to the values of spt
#' @param mu The fitted mean for each observation. Note these values must lie
#' strictly within the range of the support.
#' @param sampprobs Matrix of sampling probabilities. The number of rows should
#' equal the number of observations, and the number of columns should equal
#' the number of unique observed support points.
#' @param ySptIndex Vector containing index of each obervation's response value
#' within the \code{spt} vector. This is only needed to calculate the log-likelihood
#' after each update.
#' @param thetaStart Vector of starting values. One value per observation. If
#' \code{NULL}, zero is used as the starting value for each observation.
#' @param thetaControl Object of class \code{thetaControl}, which is a list of
#' control arguments returned by the \code{thetaControl} function.
#'
#' @return List containing the following:
#' \itemize{
#' \item \code{theta} Updated values.
#' \item \code{fTilt} Matrix containing the exponentially tilted distribution for each
#' observation, i.e. f(y|X=x). Each column corresponds to an observation and sums to one.
#' \item \code{bPrime} Vector containing the mean of the exponentially tilted distribution
#' for each observation. Should match \code{mu} argument very closely.
#' \item \code{bPrime2} Vector containing the variance of the exponentially tilted
#' distribution for each observation.
#' \item \code{fTiltSW} Matrix containing the exponentially tilted distribution for each
#' observation, conditional on that observation being sampled, i.e. f(y|X=x, S=1).
#' If \code{sampprobs=NULL}, then \code{fTiltSW} matches \code{fTilt}.
#' \item \code{bPrimeSW} Vector containing the mean for each observation, conditional
#' on that observation being sampled. If \code{sampprobs=NULL}, then \code{bPrimeSW}
#' matches \code{bPrime}.
#' \item \code{bPrime2SW} Vector containing the variance for each observation, conditional
#' on that observation being sampled. If \code{sampprobs=NULL}, then \code{bPrime2SW}
#' matches \code{bPrime2}.
#' \item \code{llik} Semiparametric log-likelihood, evaluated at the current beta
#' and f0 values. If sampling weights are used, then the log-likelihood is conditional
#' on each observation being sampled.
#' \item \code{conv} Convergence indicator.
#' \item \code{iter} Number of iterations until convergence was reached.
#' }
#'
#' @keywords internal
getTheta <- function(spt, f0, mu, sampprobs, ySptIndex, thetaStart=NULL, thetaControl=theta.control())
{
    ## Extract control arguments
    if (class(thetaControl) != "thetaControl")
        stop("thetaControl must be an object of class \'thetaControl\' returned by
             thetaControl() function.")
    logit <- thetaControl$logit
    eps <- thetaControl$eps
    maxiter <- thetaControl$maxiter
    maxhalf <- thetaControl$maxhalf
    maxtheta <- thetaControl$maxtheta
    logsumexp <- thetaControl$logsumexp

    ## Define value from inputs
    sptN <- length(spt)
    m <- min(spt)
    M <- max(spt)
    n <- length(mu)

    ## Format arguments
    spt <- as.vector(spt)
    f0 <- as.vector(f0)
    mu <- as.vector(mu)
    if (!is.null(thetaStart)) {
        thetaStart <- as.vector(thetaStart)
    } else {
        thetaStart <- rep(0, n)
    }

    ## Argument checks
    if (length(f0) != sptN)
        stop("spt and f0 must be vectors of equal length.")
    if (any(f0 < 0))
        stop("f0 values cannot be negative.")
    if (min(mu)<m || max(mu)>M)
        stop("mu starting values must lie within the range of spt.")
    if (length(thetaStart) != n)
        stop("thetaStart must be a vector with length equal length(mu)")

    ## Value does not change
    gMu <- g(mu, m, M)

    ## Initialize values
    theta <- thetaStart  # initial values required
    thetaOld <- bPrimeErrOld <- rep(NA, n)
    conv <- rep(FALSE, n)
    maxedOut <- rep(FALSE, n)

    if (logsumexp) {
        logf0 <- log(f0)
        logfUnstd <- logf0 + tcrossprod(spt, theta)
        logb <- apply(logfUnstd, 2, logSumExp)
        fTilt <- exp(logfUnstd - rep(logb, each=sptN))
        normfact <- colSums(fTilt)
        normss <- which(normfact != 1)
        fTilt[, normss] <- fTilt[, normss, drop=FALSE] / rep(normfact[normss], each=sptN)
    } else {
        fUnstd <- f0 * exp(tcrossprod(spt, theta))  # |spt| x n matrix of tilted f0 values
        b <- colSums(fUnstd)
        fTilt <- fUnstd / rep(b, each=sptN)  # normalized
    }
    bPrime <- colSums(spt*fTilt)  # mean as a function of theta
    bPrime2 <- colSums(outer(spt, bPrime, "-")^2 * fTilt)  # variance as a function of theta
    bPrimeErr <- bPrime - mu  # used to assess convergence

    ## Update theta until convergence
    conv <- (abs(bPrimeErr) < eps) |
        (theta==maxtheta & bPrimeErr<0) |
        (theta==-maxtheta & bPrimeErr>0)
    s <- which(!conv)
    iter <- 0
    while(length(s)>0 && iter<maxiter) {
        iter <- iter + 1
        bPrimeErrOld[s] <- bPrimeErr[s]  # used to assess convergence

        ## 1) Update theta
        thetaOld[s] <- theta[s]
        if (logit) {
            tPrimeS <- (M-m) / ((bPrime[s]-m) * (M-bPrime[s])) * bPrime2[s]
                # t'(theta) from paper: temporary variable
                # only needed for the subset of observations that have not converged
            tPrimeS[is.na(tPrimeS) | tPrimeS==Inf] <- 0
                # If bPrime is on the boundary, then bPrime2 should be zero.
                # Exceptions are due to rounding error.
            thetaS <- theta[s] - (g(bPrime[s], m, M) - gMu[s]) / tPrimeS
        } else {
            thetaS <- theta[s] - bPrimeErr[s] / bPrime2[s]
        }
        thetaS[thetaS > maxtheta] <- maxtheta
        thetaS[thetaS < -maxtheta] <- -maxtheta
        theta[s] <- thetaS

        ## 2) Update fTilt, bPrime, and bPrime2 and take half steps if bPrimeErr not improved
        ss <- s
        nhalf <- 0
        while(length(ss)>0 && nhalf<maxhalf) {
            ## 2a) Update fTilt, bPrime, and bPrime2
            if (logsumexp) {
                logfUnstd[, ss] <- logf0 + tcrossprod(spt, theta[ss])
                logb[ss] <- apply(logfUnstd[, ss, drop=FALSE], 2, logSumExp)
                fTilt[, ss] <- exp(logfUnstd[, ss, drop=FALSE] - rep(logb[ss], each=sptN))
                normfact <- colSums(fTilt[, ss, drop=FALSE])
                normss <- which(normfact != 1)
                fTilt[, ss[normss]] <- fTilt[, ss[normss], drop=FALSE] / rep(normfact[normss], each=sptN)
            } else {
                fUnstd[, ss] <- f0*exp(tcrossprod(spt, theta[ss]))  # |spt| x n matrix of tilted f0 values
                b[ss] <- colSums(fUnstd[, ss, drop=FALSE])
                fTilt[, ss] <- fUnstd[, ss, drop=FALSE] / rep(b[ss], each=sptN)  # normalized
            }
            bPrime[ss] <- colSums(spt*fTilt[, ss, drop=FALSE])  # mean as a function of theta
            bPrime2[ss] <- colSums(outer(spt, bPrime[ss], "-")^2 * fTilt[, ss, drop=FALSE])  # variance as a function of theta
            bPrimeErr[ss] <- bPrime[ss] - mu[ss]  # used to assess convergence

            ## 2b) Take half steps if necessary
            ss <- ss[abs(bPrimeErr[ss]) > abs(bPrimeErrOld[ss])]
            if (length(ss) > 0) nhalf <- nhalf + 1
            theta[ss] <- (theta[ss] + thetaOld[ss]) / 2
        }
        ## If maximum half steps are exceeded, set theta to previous value
        maxedOut[ss] <- TRUE
        theta[ss] <- thetaOld[ss]

        ## 3) Check convergence
        conv[s] <- (abs(bPrimeErr[s]) < eps) |
            (theta[s]==maxtheta & bPrimeErr[s]<0) |
            (theta[s]==-maxtheta & bPrimeErr[s]>0)
        s <- s[!conv[s] & !maxedOut[s]]
    }

    ## Obtain conditional distribution of sampled observations, if sampling weights are included
    if (!(is.null(sampprobs))) {
        fTiltSW <- fTilt * t(sampprobs)
        fTiltSW <- fTiltSW / rep(colSums(fTiltSW), each=nrow(fTiltSW))
        bPrimeSW <- colSums(spt * fTiltSW)
        bPrime2SW <- colSums(outer(spt, bPrimeSW, "-")^2 * fTiltSW)
    } else {
        fTiltSW <- fTilt
        bPrimeSW <- bPrime
        bPrime2SW <- bPrime2
    }

    ## Calculate log-likelihood
    llik <- sum(log(fTiltSW[cbind(ySptIndex, seq_along(ySptIndex))]))

    list(theta=theta, fTilt=fTilt, bPrime=bPrime, bPrime2=bPrime2,
         fTiltSW=fTiltSW, bPrimeSW=bPrimeSW, bPrime2SW=bPrime2SW,
         llik=llik, conv=conv, iter=iter)
}
