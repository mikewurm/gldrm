simData <- function(n, p, betaMax, link, ySim)
{
    beta <- c(0, seq(-betaMax, betaMax, length.out=p))
    x <- cbind(1, matrix(rnorm(n*p), nrow=n))
    eta <- as.vector(x %*% beta)
    mu <- link$linkinv(eta)
    y <- ySim(n, mu)

    list(x=x, y=y, eta=eta, mu=mu, beta=beta)
}

test_that("gldrm and gldrmFit match", {
    set.seed(100)
    lf <- stats::make.link("log")
    ySim <- function(n, mu) rgamma(n, 1, scale=mu)
    dat <- simData(n=100, p=5, betaMax=0, lf, ySim)
    m1 <- gldrm(dat$y~dat$x-1, data=NULL, link="log")
    m2 <- gldrm:::gldrmFit(x=dat$x, y=dat$y, linkfun=lf$linkfun, linkinv=lf$linkinv, mu.eta=lf$mu.eta)
    m1[c("formula", "data", "link", "offset")] <- NULL
    expect_equal(m1, m2, tolerance=1e-3, check.attributes=FALSE)
})

test_that("gldrm matches intercept-only (empirical distribution) model", {
    set.seed(100)
    ySim <- function(n, mu) rpois(n, 1)
    dat <- simData(n=100, p=0, betaMax=0, make.link("identity"), ySim)

    m1 <- as.vector(table(dat$y)) / length(dat$y)
    m2 <- gldrm(dat$y~dat$x-1, data=NULL, link="identity")  # link function doesn't matter with no covariates
    fcorr <- c(gldrm:::getTheta(spt=m2$spt, f0=m2$f0, mu=mean(dat$y), sampprobs=NULL, ySptIndex=1)$fTilt)  # tilt f0 so mean = mean(y)

    ## this is an intercept-only model, so all observations have fitted mean equal to mean(y)
    expect_equal(mean(dat$y), m2$beta, check.attributes=FALSE)
    ## f0 should match response frequency table
    expect_equal(m1, fcorr)
})

test_that("gldrm matches logistic regression", {
    set.seed(100)
    ySim <- function(n, mu) rbinom(n, 1, mu)
    dat <- simData(n=100, p=5, betaMax=.5, make.link("logit"), ySim)

    m1 <- with(dat, glm(y ~ -1 + x, family=binomial(link="logit")))
    m2 <- gldrm(dat$y~dat$x-1, data=NULL, link="logit")

    ## SPGLM should match logistic regression coefficient estimates
    ## (semiparametric model is identical to fully parametric in this case)
    expect_equal(as.vector(coef(m1)), m2$beta, tolerance=1e-7, check.attributes=FALSE)
})

test_that("Can handle muHat on boundary of spt", {
    n <- 10
    y <- rep(c(0, 1), each=n/2)
    x <- cbind(1, y)

    m1 <- gldrm(y~x-1, data=NULL, link="identity")
    expect_equal(m1$beta, c(0, 1), check.attributes=FALSE)

    m2 <- gldrm(y~x-1, data=NULL, link="logit")
    expect_equal(m2$mu, y)
})

test_that("Can handle singular covariate matrix", {
    n <- 10
    y <- rep(c(0, 1), each=n/2)
    x <- matrix(1, nrow=n, ncol=2)

    m1 <- gldrm(y~x-1, data=NULL, link="identity")
    expect_equal(m1$beta, c(.5, NA), check.attributes=FALSE)

    m2 <- gldrm(y~x-1, data=NULL, link="logit")
    expect_equal(m2$beta, c(0, NA), check.attributes=FALSE)
})
