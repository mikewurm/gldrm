test_that("getTheta returns correct theta.", {
    set.seed(100)
    n <- 100
    sptN <- 10
    sptMax <- 1
    spt <- seq(-sptMax, sptMax, length.out=sptN)
    f0 <- rep(1/sptN, sptN)
    mu <- seq(-sptMax, sptMax, length.out=n)  # includes mu on boundary on support
    th <- gldrm:::getTheta(spt=spt, f0=f0, mu=mu, sampprobs=NULL, ySptIndex=1, thetaStart=NULL)
    expect_equal(mu, th$bPrime)
})
