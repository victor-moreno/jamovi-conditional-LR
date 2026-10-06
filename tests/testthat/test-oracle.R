library(survival)

test_that("coefficients, SEs and CIs match survival::clogit", {
    d <- infertData()
    r <- runClr(d, covs = c("spontaneous", "induced"),
                blocks = list(list("spontaneous", "induced")), coefDetails = TRUE)
    cl <- clogit(case == "1" ~ spontaneous + induced + strata(stratum), data = d)
    tb <- coefDF(r)
    expect_equal(tb$est, unname(coef(cl)), tolerance = 1e-8)
    expect_equal(tb$se, unname(sqrt(diag(vcov(cl)))), tolerance = 1e-8)
    expect_equal(tb$or, unname(exp(coef(cl))), tolerance = 1e-8)
    ci <- confint(cl)
    expect_equal(tb$lower, unname(exp(ci[, 1])), tolerance = 1e-8)
    expect_equal(tb$upper, unname(exp(ci[, 2])), tolerance = 1e-8)
    r90 <- runClr(d, covs = c("spontaneous", "induced"),
                  blocks = list(list("spontaneous", "induced")), ciWidth = 90)
    expect_equal(coefDF(r90)$lower, unname(exp(confint(cl, level = 0.9)[, 1])), tolerance = 1e-8)
    expect_equal(tb$p, unname(summary(cl)$coefficients[, 5]), tolerance = 1e-8)
})

test_that("fit measures and overall test follow the conditional log-likelihood", {
    d <- infertData()
    r <- runClr(d, covs = c("spontaneous", "induced"),
                blocks = list(list("spontaneous", "induced")),
                modelTest = TRUE, bic = TRUE, pseudoR2 = c("r2mf", "r2cs", "r2n"))
    cl <- clogit(case == "1" ~ spontaneous + induced + strata(stratum), data = d)
    fit <- r$modelFit$asDF
    ll <- cl$loglik
    n <- nrow(d)
    expect_equal(fit$dev, -2 * ll[2], tolerance = 1e-8)
    expect_equal(fit$aic, AIC(cl), tolerance = 1e-8)
    expect_equal(fit$bic, BIC(cl), tolerance = 1e-8)
    expect_equal(fit$chi, unname(summary(cl)$logtest["test"]), tolerance = 1e-8)
    expect_equal(fit$df, 2)
    expect_equal(fit$r2mf, 1 - ll[2] / ll[1], tolerance = 1e-8)
    r2cs <- 1 - exp(2 * (ll[1] - ll[2]) / n)
    expect_equal(fit$r2cs, r2cs, tolerance = 1e-8)
    expect_equal(fit$r2n, r2cs / (1 - exp(2 * ll[1] / n)), tolerance = 1e-8)
})

test_that("model comparison is the LR test between nested clogit fits", {
    d <- infertData()
    r <- runClr(d, covs = c("spontaneous", "induced"),
                blocks = list(list("spontaneous"), list("induced")))
    c1 <- clogit(case == "1" ~ spontaneous + strata(stratum), data = d)
    c2 <- clogit(case == "1" ~ spontaneous + induced + strata(stratum), data = d)
    a <- anova(c1, c2)
    cmp <- r$modelComp$asDF
    expect_equal(cmp$chi, a$Chisq[2], tolerance = 1e-8)
    expect_equal(cmp$df, a$Df[2])
    expect_equal(cmp$p, a$`Pr(>|Chi|)`[2], tolerance = 1e-8)
})

test_that("omnibus LR tests equal refits without the term", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced", c("spont", "induced"))),
                refLevels = list(list(var = "case", ref = "0"), list(var = "spont", ref = "0")),
                omni = TRUE)
    lrt <- r$models[[1]]$lrt$asDF
    full <- clogit(case == "1" ~ spont * induced + strata(stratum), data = d)
    main <- clogit(case == "1" ~ spont + induced + strata(stratum), data = d)
    expect_equal(lrt$chi[3], 2 * (full$loglik[2] - main$loglik[2]), tolerance = 1e-8)
    expect_equal(lrt$df[3], 2)

    r1 <- runClr(d, covs = "induced", blocks = list(list("induced")),
                 omni = TRUE, modelTest = TRUE)
    expect_equal(r1$models[[1]]$lrt$asDF$chi, r1$modelFit$asDF$chi, tolerance = 1e-8)
})

test_that("factor coefficients follow the chosen reference levels", {
    d <- infertData()
    ref0 <- runClr(d, factors = "spont", blocks = list(list("spont")),
                   refLevels = list(list(var = "case", ref = "0"), list(var = "spont", ref = "0")))
    cl <- clogit(case == "1" ~ spont + strata(stratum), data = d)
    expect_equal(coefDF(ref0)$est[-1], unname(coef(cl)), tolerance = 1e-8)

    ref2 <- runClr(d, factors = "spont", blocks = list(list("spont")),
                   refLevels = list(list(var = "case", ref = "0"), list(var = "spont", ref = "2")))
    d$spont2 <- relevel(d$spont, ref = "2")
    cl2 <- clogit(case == "1" ~ spont2 + strata(stratum), data = d)
    expect_equal(coefDF(ref2)$est[-1], unname(coef(cl2)), tolerance = 1e-8)

    flip <- runClr(d, covs = "induced", blocks = list(list("induced")),
                   refLevels = list(list(var = "case", ref = "1")))
    base <- runClr(d, covs = "induced", blocks = list(list("induced")))
    expect_equal(coefDF(flip)$est, -coefDF(base)$est, tolerance = 1e-6)
})

test_that("an ordinal factor is coded by treatment contrasts, not polynomials", {
    d <- infertData()
    d$spontO <- factor(d$spontaneous, ordered = TRUE)
    r <- runClr(d, factors = "spontO", blocks = list(list("spontO")),
                refLevels = list(list(var = "case", ref = "0"), list(var = "spontO", ref = "0")))
    cl <- clogit(case == "1" ~ spont + strata(stratum), data = d)
    expect_equal(coefDF(r)$est[-1], unname(coef(cl)), tolerance = 1e-8)
})

test_that("approximate likelihoods equal the exact one with one case per set", {
    d <- infertData()
    args <- list(covs = "induced", blocks = list(list("induced")))
    ex <- do.call(runClr, c(list(d, ties = "exact"), args))
    ef <- do.call(runClr, c(list(d, ties = "efron"), args))
    br <- do.call(runClr, c(list(d, ties = "breslow"), args))
    expect_equal(coefDF(ef)$est, coefDF(ex)$est, tolerance = 1e-8)
    expect_equal(coefDF(br)$est, coefDF(ex)$est, tolerance = 1e-8)
})

test_that("several cases per set: exact matches clogit, efron differs", {
    set.seed(1)
    n <- 40
    d <- data.frame(stratum = factor(rep(seq_len(n), each = 5)),
                    case = factor(rep(c(1, 1, 0, 0, 0), n)),
                    x = rnorm(5 * n))
    d$x[d$case == "1"] <- d$x[d$case == "1"] + 0.7
    r <- runClr(d, covs = "x", blocks = list(list("x")))
    rEf <- runClr(d, covs = "x", blocks = list(list("x")), ties = "efron")
    cl <- clogit(case == "1" ~ x + strata(stratum), data = d)
    expect_equal(coefDF(r)$est, unname(coef(cl)), tolerance = 1e-8)
    expect_false(isTRUE(all.equal(coefDF(rEf)$est, coefDF(r)$est, tolerance = 1e-8)))
})

test_that("collinearity statistics match car::vif on the same fit", {
    skip_if_not_installed("car")
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced")), collin = TRUE)
    cl <- clogit(case == "1" ~ spont + induced + strata(stratum), data = d)
    v <- suppressWarnings(car::vif(cl))
    col <- r$models[[1]]$assump$collin$asDF
    expect_equal(col$vif, unname(v[1:2, 1]), tolerance = 1e-8)
    expect_equal(col$tol, 1 / col$vif)
})
