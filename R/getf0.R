################################################################################
# These functions are only used for the 2 term approximate information
# (method = "approx2")
################################################################################

## Computes inverse of a 2x2 matrix
invert2by2 <- function(m) {
    matrix(c(m[4], -m[2], -m[3], m[1]), nrow=2) /
        (m[1] * m[4] - m[2] * m[3])
}

## Computes (A+BCB')^{-1}, where Ainv is available and B is rank 1 or 2
## Adiag is an indicator of whether A is a diagonal matrix
woodbury <- function(Ainv, Cinv, B) {
    AinvB <- Ainv %*% B
    mid <- Cinv + crossprod(B, AinvB)
    midinv <- invert2by2(mid)
    inv <- Ainv - AinvB %*% tcrossprod(midinv, AinvB)
    inv
}
################################################################################

#' Control arguments for f0 update algorithm
#'
#' This function returns control arguments for the \eqn{f_0} update algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold. The update has converged when the relative
#' change in log-likelihood between iterations is less than \code{eps}.
#' absolute change is less than \code{thesh}.
#' @param maxiter Maximum number of iterations allowed.
#' @param maxhalf Maximum number of half steps allowed per iteration if
#' log-likelihood does not improve between iterations.
#' @param maxlogstep Maximum optimization step size allowed on the
#' \code{log(f0)} scale.
#'
#' @return Object of S3 class "f0Control", which is a list of control arguments.
#'
#' @export
f0.control <- function(eps=1e-10, maxiter=1000, maxhalf=20, maxlogstep=2)
{
    f0Control <- as.list(environment())
    class(f0Control) <- "f0Control"
    f0Control
}

#' f0 optimization routine
#'
#' @param y Vector of response values.
#' @param spt Vector of unique observed support points in the response.
#' @param ySptIndex Index of each \code{y} value within \code{spt}.
#' @param sptFreq Vector containing frequency of each \code{spt} value.
#' @param sampprobs Optional matrix of sampling probabilities.
#' @param mu Fitted mean for each observation. Only used if \code{sampprobs=NULL}.
#' @param mu0 Mean constraing for f0.
#' @param f0Start Starting f0 values. (Typically the estimate from the previous
#' iteration.)
#' @param thStart Starting theta values. Needs to be a list of values matching
#' the output of the \code{getTheta} function.
#' @param thetaControl A "thetaControl" object returned from the \code{theta.control}
#' function.
#' @param f0Control An "f0Control" object returned from the \code{f0.control}
#' function.
#' trace Logical. If TRUE, then progress is printed to terminal at each iteration.
#'
#' @return A list containing the following:
#' \itemize{
#' \item \code{f0} Updated values.
#' \item \code{llik} Updated log-likelihood.
#' \item \code{th} Updated list returned from the \code{getTheta} function.
#' \item \code{conv} Convergence indicator.
#' \item \code{iter} Number of iterations until convergence.
#' \item \code{nhalf} The number of half steps taken on the last iteration if the
#' initial BFGS update did not improve the log-likelihood.
#' \item \code{score.log} Score function with respect to log(f0) at convergence.
#' \item \code{info.log} Information matrix with respect to log(f0) at convergence.
#' }
#'
#' @keywords internal
getf0 <- function(y, spt, ySptIndex, sptFreq, sampprobs, mu, mu0, f0Start, thStart,
                  thetaControl=theta.control(), f0Control=f0.control(), trace=FALSE)
{
    ## Extract theta control arguments
    if (class(f0Control) != "f0Control")
        stop("f0Control must be an object of class f0Control returned by f0Control() function.")
    eps <- f0Control$eps
    maxiter <- f0Control$maxiter
    maxhalf <- f0Control$maxhalf
    maxlogstep <- f0Control$maxlogstep

    f0 <- f0Start  # assumes sum(f0Start) = 1 and sum(f0Start * spt) = mu0
    th <- thStart
    llik <- th$llik
    score.log <- NULL
    if (is.null(sampprobs)) {
        smm <- outer(spt, mu, "-")
        ymm <- y - mu
        yeqmu <- which(abs(ymm) < 1e-15)
    }

    conv <- FALSE
    iter <- 0
    while (!conv && iter<maxiter) {
        iter <- iter + 1

        # Score calculation
        score.logOld <- score.log
        if (!is.null(sampprobs)) {
            smm <- outer(spt, th$bPrimeSW, "-")
            ymm <- y - th$bPrimeSW
            yeqmu <- which(abs(ymm) < 1e-15)
        }
        fTiltSWSums <- rowSums(th$fTiltSW)
        smmfTiltSW <- smm * th$fTiltSW
        ystd <- ymm / th$bPrime2SW
        ystd[yeqmu] <- 0  # prevent 0/0
        score.logT1 <- sptFreq
        score.logT2 <- fTiltSWSums
        score.logT3 <- c(smmfTiltSW %*% ystd)
        score.log <- score.logT1 - score.logT2 - score.logT3

        # Inverse info, score step, and f0 step are on the log scale (score is not)
        if (iter == 1) {
            d1 <- min(fTiltSWSums)  # max inverse diagonal of first information term, on log scale
            d2 <- max(abs(score.log)) / maxlogstep
            d <- max(d1, d2)
            infoinvBFGS.log <- diag(1/d, nrow=length(f0))
        } else {
            scorestep.log <- score.log - score.logOld
            f0step.log <- log(f0) - log(f0old)
            sy <- sum(f0step.log * scorestep.log)
            yiy <- c(crossprod(scorestep.log, infoinvBFGS.log %*% scorestep.log))
            iys <- tcrossprod(infoinvBFGS.log %*% scorestep.log, f0step.log)
            infoinvBFGS.log <- infoinvBFGS.log + ((yiy - sy) / sy^2) * tcrossprod(f0step.log) - (1 / sy) * (iys + t(iys))
        }
        logstep <- c(infoinvBFGS.log %*% score.log)

        # Cap log(f0) step size
        logstep.max <- max(abs(logstep))
        if (logstep.max > maxlogstep)
            logstep <- logstep * (maxlogstep / logstep.max)

        # Save values from previous iteration
        f0old <- f0
        thold <- th
        llikold <- llik

        # Take update step
        f0 <- exp(log(f0) + logstep)
        # Scale and tilt f0
        f0 <- f0 / sum(f0)
        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1,
                       thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
        # Update theta and likelihood
        thold <- th
        llikold <- llik
        th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
                       thetaStart=th$theta, thetaControl=thetaControl)
        llik <- th$llik
        conv <- abs((llik - llikold) / (llik + 1e-100)) < eps

        # If log-likelihood does not improve, change step direction to be along gradient
        # Take half steps until likelihood improves
        # Continue taking half steps until log likelihood no longer improves
        nhalf <- 0
        if (llik<llikold) {
            llikprev <- -Inf
            while ((llik<llikold || llik>llikprev) && nhalf<maxhalf) {
                nhalf <- nhalf + 1

                # Set previous values
                llikprev <- llik
                thprev <- th
                f0prev <- f0
                infoinvBFGS.logprev <- infoinvBFGS.log

                f0 <- exp((log(f0) + log(f0old)) / 2)
                f0 <- f0 / sum(f0)
                f0 <- getTheta(spt=spt, f0=f0, mu=mu0, sampprobs=NULL, ySptIndex=1,
                               thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
                th <- getTheta(spt=spt, f0=f0, mu=mu, sampprobs=sampprobs, ySptIndex=ySptIndex,
                               thetaStart=th$theta, thetaControl=thetaControl)
                llik <- th$llik
                infoinvBFGS.log <- infoinvBFGS.log / 2
            }

            if (llik < llikprev) {
                nhalf <- nhalf - 1
                llik <- llikprev
                th <- thprev
                f0 <- f0prev
                infoinvBFGS.log <- infoinvBFGS.logprev
            }

            conv <- abs((llik - llikold) / (llik + 1e-100)) < eps
        }

        if (llik < llikold) {
            f0 <- f0old
            th <- thold
            llik <- llikold
            conv <- TRUE
        }

        if (trace) {
            printout <- paste0("iter ", iter, ": llik=", llik)
            if (nhalf > 0)
                printout <- paste0(printout, "; ", nhalf, " half steps")
            cat(printout, "\n")
        }
    }

    # Final score calculation
    if (!is.null(sampprobs)) {
        smm <- outer(spt, th$bPrimeSW, "-")
        ymm <- y - th$bPrimeSW
        yeqmu <- which(abs(ymm) < 1e-15)
    }
    fTiltSWSums <- rowSums(th$fTiltSW)
    smmfTiltSW <- smm * th$fTiltSW
    ystd <- ymm / th$bPrime2SW
    ystd[yeqmu] <- 0  # prevent 0/0
    score.logT1 <- sptFreq
    score.logT2 <- fTiltSWSums
    score.logT3 <- c(smmfTiltSW %*% ystd)
    score.log <- score.logT1 - score.logT2 - score.logT3

    # Final info calculation
    info.logT1 <- diag(fTiltSWSums)
    info.logT2 <- tcrossprod(th$fTiltSW)
    info.logT3 <- tcrossprod(smmfTiltSW, smmfTiltSW * rep(ystd, each=nrow(smmfTiltSW)))
    info.log <- info.logT1 - info.logT2 - info.logT3

    list(f0=f0, llik=llik, th=th, conv=conv, iter=iter, nhalf=nhalf,
         score.log=score.log, info.log=info.log)
}
