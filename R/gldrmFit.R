#' Control arguments for \code{gldrm} algorithm
#'
#' This function returns control arguments for the \code{gldrm} algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.

#' @param eps Convergence threshold. The fitting algorithm has converged when the
#' relative change in log-likelihood between iterations is less than \code{eps}.
#' A single iteration consists of a \code{beta} update followed by an \code{f0}
#' update.
#' @param maxiter Maximum number of iterations allowed.
#' @param returnfTiltMatrix Logical. Return nonparametric fitted probabilities for
#' each observation. This is a matrix with nrow equal to the number of
#' observations and ncol equal to the number of unique response values observed.
#' @param returnf0ScoreInfo Logical. If \code{TRUE}, the score and information for
#' \code{log(f0)} are returned as components of the "gldrm" object.
#' @param betaStart Optional vector of starting values for \code{beta}. If the
#' call to gldrm contains a formula, the values of betaStart should correspond to
#' the columns of the model matrix.
#' @param f0Start Optional vector of starting values for \code{f0}. The length
#' of the vector should be the number of unique values in the response, and the
#' vector should correspond to these values sorted in increasing order. The starting
#' values will be scaled to sum to one and tilted to have mean \code{mu0}. All values
#' should be strictly positive.
#' @param print Logical. If \code{TRUE}, the relative change in the log-likelihood
#' will be printed after each iteration.
#'
#' @return Object of S3 class "gldrmControl", which is a list of control arguments.
#'
#' @export
gldrm.control <- function(eps=1e-10, maxiter=100, returnfTiltMatrix=TRUE,
                          returnf0ScoreInfo=FALSE, print=FALSE,
                          betaStart=NULL, f0Start=NULL)
{
    gldrmControl <- as.list(environment())
    class(gldrmControl) <- "gldrmControl"
    gldrmControl
}

#' Main optimization function
#'
#' This function is called by the main \code{gldrm} function.
#'
#' @keywords internal
gldrmFit <- function(x, y, linkfun, linkinv, mu.eta, mu0=NULL, offset=NULL, sampprobs=NULL,
                     gldrmControl=gldrm.control(), thetaControl=theta.control(),
                     betaControl=beta.control(), f0Control=f0.control())
{
    ## Extract control arguments
    if (class(gldrmControl) != "gldrmControl")
        stop("gldrmControl must be an object of class \'gldrmControl\' returned by
              gldrmControl() function.")
    eps <- gldrmControl$eps
    maxiter <- gldrmControl$maxiter
    returnHess <- gldrmControl$returnHess
    returnfTiltMatrix <- gldrmControl$returnfTiltMatrix
    returnf0ScoreInfo <- gldrmControl$returnf0ScoreInfo
    print <- gldrmControl$print
    betaStart <- gldrmControl$betaStart
    f0Start <- gldrmControl$f0Start

    ## Tabulation and summary of responses used in estimating f0
    n <- length(y)
    spt <- sort(unique(y))  # observed support
    ySptIndex <- match(y, spt)  # index of each y value within support
    sptFreq <- table(ySptIndex)
    attributes(sptFreq) <- NULL

    ## Check sampprobs
    if (!is.null(sampprobs)) {
        if (!(is.vector(sampprobs) || is.matrix(sampprobs)) || !is.numeric(sampprobs) || any(sampprobs < 0))
            stop("sampprobs must be a matrix or vector of nonnegative numeric values")

        if (is.vector(sampprobs)) {
            if (length(sampprobs) != length(spt))
                stop(paste0("sampprobs vector should have length equal to the ",
                            "number of unique observed values in the response."))
            sampprobs <- matrix(sampprobs, nrow=n, ncol=length(sampprobs), byrow=TRUE)
        } else {
            # sampprobs must be a matrix
            if (nrow(sampprobs) != n)
                stop(paste0("sampprobs matrix should have row dimension equal to ",
                            "the number of observations."))
            if (ncol(sampprobs) != length(spt))
                stop(paste0("sampprobs matrix should have column dimension equal ",
                            "to the number of unique observe values in the response."))
        }
    }

    ## Initialize offset
    if (is.null(offset))
        offset <- rep(0, n)

    ## Initialize mu0 if not provided by user
    if (is.null(mu0)) {
        mu0 <- mean(y)
        # mu0 <- linkinv(0)
        # mu0 <- mean(range(y))
        # mu0 <- mean(spt)
    } else if (mu0<=min(spt) || mu0>=max(spt)) {
        stop(paste0("mu0 must lie within the range of observed values. Choose a different ",
                    "value or set mu0=NULL to use the default value, mean(y)."))
    }

    ## Initialize f0
    if (is.null(f0Start)) {
        f0 <- sptFreq / n
        if (mu0 != mean(y))
            f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1, thetaStart=0,
                           thetaControl=thetaControl)$fTilt[, 1]
    } else {
        if (length(f0Start) != length(spt))
            stop("Length of f0Start should equal number of unique values in the response.")
        if (any(f0Start <= 0))
            stop("All values in f0Start should be strictly positive.")
        f0 <- f0Start / sum(f0Start)
        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1, thetaStart=0,
                       thetaControl=thetaControl)$fTilt[, 1]
    }

    ## Initialize beta
    ## The starting values returned by lm.fit guarantee that all mu values are
    ## within the support range, even if there is no intercept.
    ## Offset could still create problems.
    lmcoef <- stats::lm.fit(x, linkfun(mu0) - offset)$coef
    if (is.null(betaStart)) {
        beta <- lmcoef
    } else {
        if (length(betaStart) != ncol(x))
            stop("Length of betaStart should equal the number of columns in the model matrix.")
      beta <- betaStart
    }

    ## Drop coefficients if x is not full rank (add NA values back at the end)
    naID <- is.na(lmcoef)
    beta <- beta[!naID]
    x <- x[, !naID, drop=FALSE]
    eta <- c(x %*% beta + offset)
    mu <- linkinv(eta)
    if (ncol(x) >= n)
    stop("gldrm requires n > p.")
    if (any(mu<min(spt) | mu>max(spt)))
    stop("Unable to find beta starting values that do not violate convex hull condition.")

    ## Get initial theta and log likelihood
    th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
                   thetaStart=NULL, thetaControl=thetaControl)
    llik <- th$llik

    conv <- FALSE
    iter <- 0
    while (!conv && iter <= maxiter)
    {
        iter <- iter+1
        betaold <- beta
        f0old <- f0
        llikold <- llik

        ## update beta (mu) and theta, with fixed f0:
        bb <- getBeta(x=x, y=y, spt=spt, ySptIndex=ySptIndex, f0=f0,
                      linkinv=linkinv, mu.eta=mu.eta, offset=offset, sampprobs=sampprobs,
                      betaStart=beta, thStart=th,
                      thetaControl=thetaControl, betaControl=betaControl)
        th <- bb$th
        llik <- bb$llik
        mu <- bb$mu
        beta <- bb$beta

        ## update f0 and theta, with fixed beta (mu)
        ff <- getf0(y=y, spt=spt, ySptIndex=ySptIndex, sptFreq=sptFreq,
                    sampprobs=sampprobs, mu=mu, mu0=mu0, f0Start=f0, thStart=th,
                    thetaControl=thetaControl, f0Control=f0Control, trace=FALSE)
        th <- ff$th
        llik <- ff$llik
        f0 <- ff$f0

        ## Check convergence
        del <- abs((llik - llikold) / llik)
        if (llik == 0) del <- 0
        conv <- del < eps

        if (print) {
            cat("iteration ", iter,
                "\nrelative change in log-likelihood = ", del,
                "  (eps = ", eps, ")\n")
        }
    }

    ## Final values
    eta <- linkfun(mu)
    dmudeta <- mu.eta(eta)
    llik <- ff$llik
    theta <- th$theta
    bPrime <- th$bPrime
    bPrime2 <- th$bPrime2
    fTilt <- th$fTilt[cbind(ySptIndex, seq_along(ySptIndex))]

    ## Compute betaHat variance
    if (!is.null(sampprobs)) {
        q <- th$bPrime2SW / th$bPrime2
        w <- dmudeta^2 / th$bPrime2 * q
        wSqrt <- sqrt(w)
    } else {
        w <- dmudeta^2 / th$bPrime2
        wSqrt <- sqrt(w)
    }
    if (any(wSqrt == Inf)) {
        ## if any weights are infinite, return all standard errors as zero
        varbeta <- matrix(0, nrow=length(beta), ncol=length(beta))
    } else {
        wtdX <- wSqrt * x
        # varbeta <- chol2inv(qr.R(qr(wtdX)))  # not stable
        # varbeta <- solve(crossprod(wtdX))  # not stable
        varbeta <- tcrossprod(backsolve(qr.R(qr(wtdX)), diag(ncol(wtdX))))
    }

    ## Compute standard errors
    seBeta <- sqrt(diag(varbeta))
    seEta <- sqrt(pmax(0, apply(x, 1, function(xx) crossprod(xx, varbeta) %*% xx)))
    seMu <- dmudeta * seEta

    ## Add NA values back into beta vector and varbeta if covariate matrix is not full rank
    nBeta <- length(beta) + sum(naID)
    betaTemp <- seBetaTemp <- rep(NA, nBeta)
    betaTemp[!naID] <- beta
    seBetaTemp[!naID] <- seBeta
    beta <- betaTemp
    seBeta <- seBetaTemp
    varbetaTemp <- matrix(NA, nrow=nBeta, ncol=nBeta)
    varbetaTemp[!naID, !naID] <- varbeta
    varbeta <- varbetaTemp

    # Inference vs. null model
    containsIntercept <- any(apply(x, 2, function(xx) all(xx == xx[1])))
    if (!containsIntercept) {
        llikNull <- lr.stat <- lr.df <- lr.pval <- NA
    } else {
        xrank <- ncol(x)  # columns of x have already been dropped, if necessary to make x full rank
        lr.df <- c(max(xrank-1, 1), n-xrank)  # force df[1] >= 1 just in case the full model is a null model
        llikNull <- sum(sptFreq * log(sptFreq / n))
        lr.stat <- 2 * (llik - llikNull) / lr.df[1]
        lr.pval <- 1 - stats::pf(lr.stat, lr.df[1], lr.df[2])
    }

    ## Return gldrm object
    attributes(beta) <- NULL
    attributes(f0) <- NULL

    fit <- list(conv=conv, iter=iter, llik=llik,
                beta=beta, mu=mu, eta=eta, f0=f0, spt=spt, mu0=mu0,
                varbeta=varbeta, seBeta=seBeta, seMu=seMu, seEta=seEta,
                theta=theta, bPrime=bPrime, bPrime2=bPrime2, fTilt=fTilt, sampprobs=sampprobs,
                llikNull=llikNull, lr.stat=lr.stat, lr.df=lr.df, lr.pval=lr.pval)

    if (returnfTiltMatrix)
        fit$fTiltMatrix <- t(th$fTilt)

    if (returnf0ScoreInfo) {
        fit$score.logf0 <- ff$score.log
        fit$info.logf0 <- ff$info.log
    }

    class(fit) <- "gldrm"
    fit
}
