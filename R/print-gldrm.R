#' Print summary of gldrm fit
#'
#' Prints fitted coefficients and standard errors, along with a likelihood ratio
#' test against the null model.
#'
#' @param x S3 object of class "gldrm", returned from the \code{gldrm} function.
#' @param digits Number of digits for rounding.
#' @param ... Unused. Additional arguments for print method.
#'
#' @export
print.gldrm <- function(x, digits=3, ...)
{
    # Extract items from gldrm object
    beta <- x$beta
    seBeta <- x$seBeta
    lr.stat <- signif(x$lr.stat, digits)
    lr.pval <- ifelse(x$lr.pval < 2e-16, "< 2e-16", signif(x$lr.pval, digits))
    lr.df <- x$lr.df

    # Create coefficient table
    wald.stat <- beta / seBeta
    wald.pval <- 2 * (1 - stats::pt(abs(wald.stat), lr.df[2]))
    # coefTable <- cbind(beta, seBeta, wald.stat, wald.pval)
    coefTable <- data.frame(beta, seBeta, wald.stat, wald.pval)
    rownames(coefTable) <- names(beta)
    colnames(coefTable) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

    # Print summary output
    cat("\nSummary of gldrm fit\n\n")
    cat("Coefficients:\n")
    printCoefmat(coefTable, digits=digits)
    cat("\nLikelihood ratio test against null model:\n")
    cat("                   F-statistic: ", lr.stat, "\n")
    cat("  Numerator degrees of freedom: ", lr.df[1], "\n")
    cat("Denominator degrees of freedom: ", lr.df[2], "\n")
    cat("                       P-value: ", lr.pval, "\n")

    NULL
}
