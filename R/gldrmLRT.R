#' Likelihood ratio test for nested models
#'
#' Performs a likelihood ratio F-test between nested gldrm models.
#' The F-statistic is calculated as \eqn{2 \times (llik - llik_0) / r}, where
#' \eqn{r} is the difference is the number of parameters between the full and null
#' models. The F-statistic has degrees of freedom \eqn{r} and \eqn{n-p}, where
#' \eqn{n} is the number of observations and \eqn{p} is the number of parameters
#' in the full model.
#'
#' @param gldrmFit The full model. Must be an object of S3 class 'gldrm' returned from
#' the \code{gldrm} function.
#' @param gldrmNull The sub-model being tested under the null hypotheses.
#' Must be an object of S3 class 'gldrm' returned from the \code{gldrm} function.
#'
#' @return An S3 object of class 'gldrmLRT', containing numerator and denominator
#' degrees of freedom, an F-statistic, and a p-value.
#'
#' @examples
#' data(iris, package="datasets")
#'
#' ### Fit gldrm with all variables
#' fit <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
#'              data=iris, link="log")
#'
#' ### Fit gldrm without the categorical variable "Species"
#' fit0 <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width,
#'               data=iris, link="log")
#'
#' ### Likelihood ratio test for the nested models
#' lrt <- gldrmLRT(fit, fit0)
#' lrt
#'
#' @export
gldrmLRT <- function(gldrmFit, gldrmNull)
{
    beta <- gldrmFit$beta
    beta0 <- gldrmNull$beta
    p <- sum(!is.na(beta))
    p0 <- sum(!is.na(beta0))
    n <- length(gldrmFit$mu)
    n0 <- length(gldrmNull$mu)
    llik <- gldrmFit$llik
    llik0 <- gldrmNull$llik

    if (class(gldrmFit) != "gldrm")
        stop("gldrmFit must be an S3 object of class gldrm.")
    if (class(gldrmNull) != "gldrm")
        stop("gldrmNull must be an S3 object of class gldrm.")
    if (p0 >= p)
        stop("gldrmNull must be a sub-model of gldrmFit")
    if (n != n0)
        stop("gldrm and gldrmNull have a different number of observations.")
    if (!all(names(beta0) %in% names(beta)))
        warning(paste0("Coefficient names of the nested model are not a subset of names ",
                       "in the full model. Make sure the models are nested."))

    df <- c(p-p0, n-p)
    fstat <- 2 * (llik - llik0) / df[1]
    pval <- 1 - stats::pf(fstat, df[1], df[2])

    lrt <- list(df=df, fstat=fstat, pval=pval)
    class(lrt) <- "gldrmLRT"
    lrt
}

#' Print likelihood ratio test results
#'
#' Print method for gldrmLRT objects. Prints results of a likelihood ratio F-test
#' between nested models.
#'
#' @param x S3 object of class 'gldrmLRT', returned from the \code{gldrmLRT} function.
#' @param digits Number of digits for rounding.
#' @param ... Not used. Additional arguments for print method.
#'
#' @export
print.gldrmLRT <- function(x, digits=3, ...)
{
    if (x$pval < 2e-16) {
        pval <- "< 2e-16"
    } else {
        pval <- signif(x$pval, digits)
    }

    cat("\nLikelihood ratio test:\n\n")
    cat("                   F-statistic: ", signif(x$fstat, digits), "\n")
    cat("  Numerator degrees of freedom: ", x$df[1], "\n")
    cat("Denomicator degrees of freedom: ", x$df[2], "\n")
    cat("                       P-value: ", pval, "\n")

    return(NULL)
}
