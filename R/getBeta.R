#' Control arguments for \eqn{\beta} update algorithm
#'
#' This function returns control arguments for the \eqn{\beta} update algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold. The update has converged when the relative
#' change in log-likelihood between iterations is less than \code{eps}.
#' Only applies if \code{maxiter>1}.
#' @param maxiter Maximum number of iterations allowed.
#' @param maxhalf Maximum number of half steps allowed per iteration if
#' log-likelihood does not improve.
#'
#' @return Object of S3 class "betaControl", which is a list of control arguments.
#'
#' @export
beta.control <- function(eps=1e-10, maxiter=1, maxhalf=10)
{
    betaControl <- as.list(environment())
    class(betaControl) <- "betaControl"
    betaControl
}

#' Beta optimization routing
#'
#' @param x Covariate matrix.
#' @param y Response vector.
#' @param spt Vector of unique observed support points in the response.
#' @param ySptIndex Index of each \code{y} value within the \code{spt} vector.
#' @param f0 Current values of f0.
#' @param linkinv Inverse link function.
#' @param mu.eta Deriviative of inverse link function.
#' @param offset Vector of known offset values to be added to the linear
#' combination (x' beta) for each observation. Mostly intended for likelihood ratio
#' and score confidence intervals.
#' @param sampprobs Optional matrix of sampling probabilities.
#' @param betaStart Starting values for beta (typically the estimates from the
#' previous iteration).
#' @param thStart Starting theta values. Needs to be a list of values matching
#' the output of the \code{getTheta} function.
#' @param thetaConrol A "thetaControl" object returned from the \code{theta.control}
#' function.
#' @param betaControl A "betaControl" object returned from the \code{beta.control}
#' function.
#'
#' @return A list containing the following:
#' \itemize{
#' \item \code{beta} Updated values.
#' \item \code{mu} Updated mean for each observation.
#' \item \code{th} Updated list returned from the \code{getTheta} function.
#' \item \code{llik} Updated log-likelihood.
#' \item \code{iter} Number of iterations until convergence. (Will always be
#' one unless \code{maxiter} is increased to something greater than one using the
#' \code{betaControl} argument.)
#' \item \code{conv} Convergence indicator. (Will always be FALSE unless
#' \code{maxiter} is increased to something greater than one using the
#' \code{betaControl} argument.)
#' }
#'
#' @keywords internal
getBeta <- function(x, y, spt, ySptIndex, f0, linkinv, mu.eta, offset, sampprobs,
                    betaStart, thStart,
                    thetaControl=theta.control(), betaControl=beta.control())
{
    ## Extract control arguments
    if (class(betaControl) != "betaControl")
      stop("betaControl must be an object of class betaControl returned by betaControl() function.")
    eps <- betaControl$eps
    maxiter <- betaControl$maxiter
    maxhalf <- betaControl$maxhalf

    sptMin <- min(spt)
    sptMax <- max(spt)
    beta <- betaStart
    th <- thStart
    llik <- th$llik

    conv <- FALSE
    maxhalfreached <- FALSE
    iter <- 0
    while (!conv && !maxhalfreached && iter < maxiter)
    {
        iter <- iter+1

        ## Update mean vector and related quantities
        eta <- c(x %*% beta + offset)
        mu <- linkinv(eta)
        dmudeta <- mu.eta(eta)
        betaold <- beta
        muold <- mu
        thold <- th
        llikold <- llik

        ## Compute weighted least squares update
        if (!is.null(sampprobs)) {
            q <- th$bPrime2SW / th$bPrime2
            w <- dmudeta^2 / th$bPrime2 * q
            ymm <- y - th$bPrimeSW
            r <- ymm / (q * dmudeta)
        } else {
            w <- dmudeta^2 / th$bPrime2
            ymm <- y - mu
            r <- ymm / dmudeta
        }
        yeqmu <- which(abs(ymm) < 1e-15)
        w[yeqmu] <- 0  # prevent 0/0
        r[yeqmu] <- 0  # prevent 0/0
        if (any(w==Inf)) break
        betastep <- unname(coef(lm.wfit(x, r, w)))
            # qr.coef(qr(wSqrt*x), wSqrt*r) no longer works for some singular matrices
        betastep[is.na(betastep)] <- 0
        ## Let q = b''*(theta) / b''(theta)
        ## W = diag{dmudeta^2 / b''(theta) * q}
        ## r = (y - b'*(theta)) / (q * dmudeta)
        ## We need to solve for beta such that I(betaHat) %*% beta = Score(betaHat),
        ## or equivalently, X'WX = X'Wr, or equivalently W^{1/2}X = W^{1/2}r.
        ## The linear system can be solved using qr.coef().

        ### Update beta and take half steps if log-likelihood does not improve
        beta <- beta + betastep
        eta <- c(x %*% beta + offset)
        mu <- linkinv(eta)
        if (min(mu)<sptMin || max(mu)>sptMax) {
            llik <- -Inf
        } else {
            th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
                           thetaStart=thold$theta, thetaControl=thetaControl)
            llik <- th$llik
        }

        nhalf <- 0
        while ((llik<llikold) && (nhalf<maxhalf)) {
            nhalf <- nhalf + 1
            beta <- (beta + betaold) / 2
            eta <- c(x %*% beta + offset)
            mu <- linkinv(eta)
            if (min(mu)<sptMin || max(mu)>sptMax) {
                llik <- -Inf
            } else {
                th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
                               thetaStart=thold$theta, thetaControl=thetaControl)
                llik <- th$llik
            }
        }

        if (llik < llikold) {
            beta <- betaold
            mu <- muold
            th <- thold
            llik <- llikold
            conv <- FALSE
            maxhalfreached <- TRUE
        } else {
            del <- (llik - llikold) / llik
            if (llik == 0) del <- 0  # consider converged if model fit is perfect
            conv <- del < eps
        }
    }

    return(list(beta=beta, mu=mu, th=th, llik=llik, iter=iter, conv=conv))
}
