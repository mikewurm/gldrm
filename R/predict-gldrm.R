#' Predict method for a gldrm object
#'
#' Obtains predicted probabilities, predicted class, or linear predictors.
#'
#' @param object S3 object of class "gldrm", returned from the \code{gldrm} function.
#' @param newdata Optional data frame. If NULL, fitted values will be obtained
#' for the training data.
#' @param type The type of prediction required.  Type "link" returns the linear
#' predictor. Type "response" returns the fitted mean. Type "terms" returns
#' a matrix giving the fitted values of each term in the model formula on the
#' linear predictor scale. Type "fTilt" returns a matrix containing the
#' fitted nonparametric distribution for each observation. Each row of the matrix
#' corresponds to an observation in \code{newdata}, and each column corresponds
#' to a unique response value in the training data.
#' @param se.fit Logical. If TRUE, standard errors are also returned. Does not apply
#' for \code{type = "fTilt"}.
#' @param offset Optional offset vector. Only used if \code{newdata} is not NULL.
#' @param ... Not used. Additional predict arguments.
#'
#' @return The object returned depends on \code{type}.
#'
#' @export
predict.gldrm <- function(object, newdata=NULL,
                          type=c("link", "response", "terms", "fTilt"),
                          se.fit=FALSE, offset=NULL, ...)
{
    type <- match.arg(type)

    # Extract fitted model estimates and link function
    beta <- object$beta
    varbeta <- object$varbeta
    seBeta <- object$seBeta
    linkinv <- object$link$linkinv
    mu.eta <- object$link$mu.eta
    spt <- object$spt
    f0 <- object$f0
    mu0 <- object$mu0

    if (is.null(newdata)) {

        if (type == "link") {
            out <- object$eta
            if (se.fit) {
                out <- list(fit=out)
                out$se.fit <- object$seEta
            }
        }

        if (type == "response") {
            out <- object$mu
            if (se.fit) {
                out <- list(fit=out)
                out$se.fit <- object$seMu
            }
        }

        if (type == "terms") {
            mf <- stats::model.frame(object$formula, data=object$data)
            x <- stats::model.matrix(attr(mf, "terms"), mf)
            attributes(x)[c("assign", "contrasts")] <- NULL
            modelTerms <- x * rep(beta, each=nrow(x))
            out <- modelTerms
            if (se.fit) {
                seModelTerms <- abs(x) * rep(seBeta, each=nrow(x))
                out <- list(fit=modelTerms)
                out$se.fit <- seModelTerms
            }
        }

        if (type == "fTilt") {
            if (!is.null(object$fTiltMatrix)) {
                out <- object$fTiltMatrix
            } else {
                th <- getTheta(spt=spt, f0=f0, mu=object$mu, 
                               sampprobs=matrix(1, nrow=length(object$mu), ncol=length(f0)), 
                               ySptIndex=rep(1, length(object$mu)),
                               thetaStart=NULL)
                # ySptIndex values don't matter; only used to compute log-likelihood
                fTiltMatrix <- t(th$fTilt)
                out <- fTiltMatrix
            }
        }

        return(out)
    }

    # Otherwise, !is.null(newdata)

    # Create model matrix and offset
    mf <- model.frame(object$formula, newdata)
    x <- stats::model.matrix(attr(mf, "terms"), mf)
    attributes(x)[c("assign", "contrasts")] <- NULL
    if (is.null(offset)) offset <- rep(0, nrow(x))
    if (length(offset) != nrow(x))
        stop("offset should be NULL or a vector with length equal to the number of observations.")

    # Check convex hull condition for out of sample data
    eta <- c(x %*% beta + offset)
    mu <- linkinv(eta)
    if (any(mu<min(spt) | mu>max(spt)))
        warning(paste0("Out of sample data violate convex hull condition: some ",
                       "observations have predicted means that lie outside the ",
                       "range support values in the training data."))

    if (type == "link") {
        out <- eta
        if (se.fit) {
            seEta <- sqrt(apply(x, 1, function(xx) crossprod(xx, varbeta) %*% xx))
            out <- list(fit=eta)
            out$se.fit <- seEta
        }
    }

    if (type == "response") {
        out <- mu
        if (se.fit) {
            seEta <- sqrt(apply(x, 1, function(xx) crossprod(xx, varbeta) %*% xx))
            dmudeta <- mu.eta(eta)
            seMu <- dmudeta * seEta
            out <- list(fit=mu)
            out$se.fit <- seMu
        }
    }

    if (type == "terms") {
        modelTerms <- x * rep(beta, each=nrow(x))
        out <- modelTerms
        if (se.fit) {
            seModelTerms <- abs(x) * rep(seBeta, each=nrow(x))
            out <- list(fit=modelTerms)
            out$se.fit <- seModelTerms
        }
    }

    if (type == "fTilt") {
        th <- getTheta(spt=spt, f0=f0, mu=mu, 
                       sampprobs=matrix(1, nrow=length(mu), ncol=length(f0)),
                       ySptIndex=rep(1, length(mu)),
                       thetaStart=NULL)
        # ySptIndex values don't matter; only used to compute log-likelihood
        fTiltMatrix <- t(th$fTilt)
        out <- fTiltMatrix
    }

    return(out)
}
