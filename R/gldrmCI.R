#' Confidence intervals for gldrm coefficients
#'
#' Calculates a Wald, likelihood ratio, or score confidence interval for a single gldrm
#' coefficient. Also calculates upper or lower confidence bounds. Wald confidence
#' intervals and bounds are calculated from the standard errors which are available
#' from the gldrm model fit. For likelihood ratio and score intervals and bounds,
#' a bisection search method is used, which takes longer to run.
#'
#' @param gldrmFit A gldrm model fit. Must be an S3 object of class "gldrm",
#' returned from the \code{gldrm} function.
#' @param term Character string containing the name of the coefficient of interest.
#' The coefficient names are the names of the beta component of the fitted model
#' object. They can also be obtained from the printed model output. Usually the
#' names match the formula syntax, but can be more complicated for categorical
#' variables and interaction terms.
#' @param test Character string for the type confidence interval. Options are
#' "Wald", "LRT" (for likelihood ratio), and "Score".
#' @param level Confidence level of the interval. Should be between zero and one.
#' @param type Character string containing "2-sided" for a two-sided confidence interval,
#' "lb" for a lower bound, or "ub" for an upper bound.
#' @param eps Convergence threshold. Only applies for
#' \code{test = "LRT"} and \code{test = "Score"}.
#' Convergence is reached when likelihood ratio p-value is within \code{eps} of
#' the target p-value, based on the level of the test. For example, a two-sided
#' 95\% confidence interval has target p-value of 0.025 for both the upper and
#' lower bounds. A 95\% confidence bound has target p-value 0.05.
#' @param maxiter The maximum number of bisection method iterations for likelihood
#' ratio intervals or bounds. For two-sided intervals, \code{maxiter} iterations
#' are allowed for each bound.
#'
#' @return An S3 object of class 'gldrmCI', which is a list of the following items.
#'
#' \itemize{
#' \item \code{term} Coefficient name.
#' \item \code{test} Type of interval or bound - Wald or likelihood ratio.
#' \item \code{level} Confidence level.
#' \item \code{type} Type of interval or bound - two-sided, upper bound, or lower
#' bound.
#' \item \code{cilo}/\code{cihi} Upper and lower interval bounds. One one of the
#' two applies for confidence bounds.
#' \item \code{iterlo}/\code{iterhi} Number of bisection iterations used. Only
#' applies for likelihood ratio intervals and bounds.
#' \item \code{pvallo}/\code{pvalhi} For likelihood ratio intervals and bounds,
#' the p-value at convergence is reported.
#' \item \code{conv} Indicator for whether the confidence interval limit or bound
#' converged.
#' }
#'
#' @examples
#' data(iris, package="datasets")
#'
#' ### Fit gldrm with all variables
#' fit <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
#'              data=iris, link="log")
#'
#' ### Wald 95% confidence interval for Sepal.Width
#' ci <- gldrmCI(fit, "Sepal.Width", test="Wald", level=.95, type="2-sided")
#' ci
#'
#' @export
gldrmCI <- function(gldrmFit, term, test=c("Wald", "LRT", "Score"), level=.95,
                    type=c("2-sided", "lb", "ub"), eps=1e-10, maxiter=100)
{
    maxhalf <- 10  # hard-coded argument; no user input
    test <- match.arg(test)
    type <- match.arg(type)
    if (class(gldrmFit) != "gldrm")
        stop("gldrmFit should be an object of class gldrm, returned from the gldrm function.")
    if (!(term %in% names(gldrmFit$beta)))
        stop("term should be the name of a name from the fitted coefficient vector.")
    if (!(level>0 && level<1))
        stop("level should be between zero and one.")
    if (!(eps > 0))
        stop("eps should be a small value greater than zero.")
    if (!is.null(gldrmFit$sampprobs) && test=="Score")
        stop("Score confidence intervals with sampling weights are not currently supported.")

    mf <- model.frame(gldrmFit$formula, gldrmFit$data)
    x <- stats::model.matrix(attr(mf, "terms"), mf)
    attributes(x)[c("assign", "contrasts")] <- NULL
    y <- stats::model.response(mf, type="numeric")
    offset <- gldrmFit$object
    link <- gldrmFit$link
    mu0 <- gldrmFit$mu0
    f0 <- gldrmFit$f0
    offset <- gldrmFit$offset
    beta <- gldrmFit$beta
    seBeta <- gldrmFit$seBeta
    llik <- gldrmFit$llik
    df2 <- gldrmFit$lr.df[2]
    id <- match(term, names(beta))
    cilo <- cihi <- NA  # initialze to NA for one-sided intervals
    iterlo <- iterhi <- NA  # initialize to NA for Wald or one-sided intervals
    pvallo <- pvalhi <- NA  # initialize to NA for Wald or one-sided intervals
    convlo <- convhi <- NA  # initialize to NA for Wald or one-sided intervals

    if (type == "2-sided") {
        pvalTarget <- (1 - level) / 2
        waldstep <- stats::qt((level + 1) / 2, df2) * seBeta[id]  # length 2-sided of Wald CI
    } else {
        pvalTarget <- 1 - level
        waldstep <- stats::qt(level, df2) * unname(seBeta[id])  # length of 1-sided Wald CI
    }

    ### Wald test CI
    if (test == "Wald" && type %in% c("2-sided", "lb"))
        cilo <- unname(beta[id]) - waldstep
    if (test == "Wald" && type %in% c("2-sided", "ub"))
        cihi <- unname(beta[id]) + waldstep

    ### Likelihood ratio CI
    if (test %in% c("LRT", "Score") && type %in% c("2-sided", "lb")) {
        cilo.hi <- unname(beta[id])
        cilo.lo <- -Inf
        cilo.try <- unname(beta[id]) - waldstep / 2  # initial step size is multiplied by 2
        pvallo.hi <- 1
        pvallo.lo <- 0
        betaStartlo.hi <- beta[-id]
        f0Startlo.hi <- f0
        iterlo <- 0
        convlo <- FALSE
        while (!convlo && iterlo<maxiter) {
            iterlo <- iterlo + 1

            if (cilo.lo == -Inf) {
                cilo.try <- cilo.try - (beta[id] - cilo.try) * 2
                # Once hi and lo estimates are obtained, the average betas must
                # satisfy the convex hull condition
                betaStart <- betaStartlo.hi
                f0Start <- f0Startlo.hi
                fit <- NULL
                nhalf <- 0
                while (is.null(fit) && nhalf<maxhalf) {
                    # First try previous solution as starting values
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cilo.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStart))
                    }, error = function(e) NULL)
                    # Next, let gldrm function try to find starting values
                    if (is.null(fit)) {
                        fit <- tryCatch({
                            gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                                  mu0=mu0, offset=offset + cilo.try * x[, id],
                                  gldrmControl=gldrm.control(f0Start=f0Start, betaStart=NULL))
                        }, error = function(e) NULL)
                    }
                    # Finally, reduce the step size if no starting values can be found
                    if (is.null(fit)) {
                        nhalf <- nhalf + 1
                        cilo.try <- (cilo.try + cilo.hi) / 2
                    }
                }
            } else {
                cilo.try <- (cilo.lo + cilo.hi) / 2
                betaStart <- (betaStartlo.hi + betaStartlo.lo) / 2
                f0Start <- (f0Startlo.hi + f0Startlo.lo) / 2
                fit <- NULL

                # Analytically, betaStart should not violate the convex hull condition,
                # but it is possible due to numerical error
                fit <- tryCatch({
                    gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                          mu0=mu0, offset=offset + cilo.try * x[, id],
                          gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStart))
                }, error = function(e) NULL)

                # Next, try betaStart.lo
                if (is.null(fit)) {
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cilo.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStartlo.lo))
                    }, error = function(e) NULL)
                }

                # Next, try betaStart.hi
                if (is.null(fit)) {
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cilo.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStartlo.hi))
                    }, error = function(e) NULL)
                }

                # Next, try default starting values
                if (is.null(fit)) {
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cilo.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=NULL))
                    }, error = function(e) NULL)
                }
            }

            if (is.null(fit)) break

            if (test == "LRT") {
                pvallo.try <- 1 - stats::pf(2 * (llik - fit$llik), 1, df2)
            } else {
                # else test == "Score"
                # not designed to incorporate sampling weights
                bPrime2 <- fit$bPrime2
                mu <- fit$mu
                eta <- fit$eta
                dmudeta <- link$mu.eta(eta)
                wSqrt <- dmudeta / sqrt(bPrime2)
                r <- (y-mu) / dmudeta
                wtdx <- wSqrt * x
                wtdr <- wSqrt * r
                betastep <- qr.coef(qr(wtdx), wtdr)  # = (X'WX)^-1 (X'Wr)
                betastep[is.na(betastep)] <- 0
                score <- c(crossprod(wtdx, wtdr))
                scorestat <- c(crossprod(score, betastep))
                pvallo.try <- 1 - stats::pf(scorestat, 1, df2)
            }

            convlo <- abs(pvallo.try -  pvalTarget) < eps

            if (pvallo.try > pvalTarget) {
                cilo.hi <- cilo.try
                pvallo.hi <- pvallo.try
                betaStartlo.hi <- fit$beta
                f0Startlo.hi <- fit$f0
            } else {
                cilo.lo <- cilo.try
                pvallo.lo <- pvallo.try
                betaStartlo.lo <- fit$beta
                f0Startlo.lo <- fit$f0
            }

            if (cilo.hi == cilo.lo) break
        }

        if (convlo) {
            cilo <- cilo.try
            pvallo <- pvallo.try
        } else {
            # If not converged, use the narrow interval estimate
            cilo <- cilo.hi
            pvallo <- pvallo.hi
        }
    }

    if (test %in% c("LRT", "Score") && type %in% c("2-sided", "ub")) {
        cihi.hi <- Inf
        cihi.lo <- unname(beta[id])
        pvalhi.hi <- 0
        pvalhi.lo <- 1
        cihi.try <- unname(beta[id]) + waldstep / 2  # initial step size is multiplied by 2
        pvalhi.hi <- 0
        pvalhi.lo <- 1
        betaStarthi.lo <- beta[-id]
        f0Starthi.lo <- f0
        iterhi <- 0
        convhi <- FALSE
        while (!convhi && iterhi<maxiter) {
            iterhi <- iterhi + 1

            if (cihi.hi == Inf) {
                cihi.try <- cihi.try + (cihi.try - beta[id]) * 2
                betaStart <- betaStarthi.lo
                f0Start <- f0Starthi.lo
                fit <- NULL
                nhalf <- 0
                while (is.null(fit) && nhalf<maxhalf) {
                    # First try previous solution as starting values
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cihi.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStart))
                    }, error = function(e) NULL)
                    # Next, let gldrm function try to find starting values
                    if (is.null(fit)) {
                        fit <- tryCatch({
                            gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                                  mu0=mu0, offset=offset + cihi.try * x[, id],
                                  gldrmControl=gldrm.control(f0Start=f0Start, betaStart=NULL))
                        }, error = function(e) NULL)
                    }
                    # Finally, reduce the step size if no starting values can be found
                    if (is.null(fit)) {
                        nhalf <- nhalf + 1
                        cihi.try <- (cihi.try + cihi.lo) / 2
                    }
                }
            } else {
                cihi.try <- (cihi.lo + cihi.hi) / 2
                betaStart <- (betaStarthi.hi + betaStarthi.lo) / 2
                f0Start <- (f0Starthi.hi + f0Starthi.lo) / 2
                fit <- NULL

                # Analytically, betaStart should not violate the convex hull condition,
                # but it is possible due to numerical error
                fit <- tryCatch({
                    gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                          mu0=mu0, offset=offset + cihi.try * x[, id],
                          gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStart))
                }, error = function(e) NULL)

                if (is.null(fit)) {
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cihi.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStarthi.lo))
                    }, error = function(e) NULL)
                }

                if (is.null(fit)) {
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cihi.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStarthi.hi))
                    }, error = function(e) NULL)
                }

                if (is.null(fit)) {
                    fit <- tryCatch({
                        gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL, link=link,
                              mu0=mu0, offset=offset + cihi.try * x[, id],
                              gldrmControl=gldrm.control(f0Start=f0Start, betaStart=NULL))
                    }, error = function(e) NULL)
                }
            }

            if (is.null(fit)) break

            if (test == "LRT") {
                pvalhi.try <- 1 - stats::pf(2 * (llik - fit$llik), 1, df2)
            } else {  # test == "Score"
                bPrime2 <- fit$bPrime2
                mu <- fit$mu
                eta <- fit$eta
                dmudeta <- link$mu.eta(eta)
                wSqrt <- dmudeta / sqrt(bPrime2)
                r <- (y-mu) / dmudeta
                wtdx <- wSqrt * x
                wtdr <- wSqrt * r
                betastep <- qr.coef(qr(wtdx), wtdr)  # = (X'WX)^-1 (X'Wr)
                betastep[is.na(betastep)] <- 0
                score <- c(crossprod(wtdx, wtdr))
                scorestat <- c(crossprod(score, betastep))
                pvalhi.try <- 1 - stats::pf(scorestat, 1, df2)
            }

            convhi <- abs(pvalhi.try -  pvalTarget) < eps

            if (pvalhi.try > pvalTarget) {
                cihi.lo <- cihi.try
                pvalhi.lo <- pvalhi.try
                betaStarthi.lo <- fit$beta
                f0Starthi.lo <- fit$f0
            } else {
                cihi.hi <- cihi.try
                pvalhi.hi <- pvalhi.try
                betaStarthi.hi <- fit$beta
                f0Starthi.hi <- fit$f0
            }

            if (cihi.hi == cihi.lo) break
        }

        if (convhi) {
            cihi <- cihi.try
            pvalhi <- pvalhi.try
        } else {
            # If not converged, use the narrow interval estimate
            cihi <- cihi.lo
            pvalhi <- pvalhi.lo
        }
    }

    ci <- list(term=term, test=test, level=level, type=type, cilo=cilo, cihi=cihi,
               iterlo=iterlo, iterhi=iterhi, pvallo=pvallo, pvalhi=pvalhi,
               convlo=convlo, convhi=convhi)
    class(ci) <- "gldrmCI"
    ci
}

#' Print confidence interval
#'
#' Print method for gldrmCI objects.
#'
#' @param x An S3 object of class 'gldrmCI'.
#' @param digits Number of digits for rounding.
#' @param ... Not used. Additional arguments for print method.
#'
#' @export
print.gldrmCI <- function(x, digits=3, ...)
{
    level <- round(x$level * 100)
    fmt <- paste0("%.", digits, "f")
    cilo <- sprintf(fmt, x$cilo)
    cihi <- sprintf(fmt, x$cihi)

    if (x$test == "Wald") test <- "Wald"
    if (x$test == "LRT") test <- "likelihood ratio"
    if (x$test == "Score") test <- "score"

    if (x$type == "2-sided") {
        cat('\n', level, "% ", test, " confidence interval for ", x$term, ":\n", sep='')
        cat("    (", cilo, ", ", cihi, ")\n", sep='')
    }

    if (x$type == "lb") {
        cat('\n', level, "% ", test, " lower bound for ", x$term, ":\n", sep='')
        cat("    ", cilo, "\n", sep='')
    }

    if (x$type == "ub") {
        cat('\n', level, "% ", test, " upper bound for ", x$term, ":\n", sep='')
        cat("    ", cihi, "\n", sep='')
    }

    return(NULL)
}
