library(survival)

test_that("covariate scaling: OR per 10 units, per SD of the model rows, labels", {
    d <- infertData()
    d$induced[c(3, 30)] <- NA
    base <- runClr(d, covs = c("age", "spontaneous", "induced"),
                   blocks = list(list("spontaneous", "induced")))
    sc <- runClr(d, covs = c("age", "spontaneous", "induced"),
                 blocks = list(list("spontaneous", "induced")),
                 covScales = list(list(var = "spontaneous", scale = "ten"),
                                  list(var = "induced", scale = "sd")))
    b <- coefDF(base); s <- coefDF(sc)
    expect_equal(s$or[1], b$or[1]^10, tolerance = 1e-8)
    expect_equal(s$upper[1], b$upper[1]^10, tolerance = 1e-8)
    sdInd <- sd(d$induced[! is.na(d$induced)])
    expect_equal(s$or[2], b$or[2]^sdInd, tolerance = 1e-8)
    expect_equal(s$p, b$p, tolerance = 1e-8)
    expect_equal(s$level, c("per 10 units", sprintf("per SD (%s)", format(signif(sdInd, 3)))))
    expect_equal(sc$modelFit$asDF$dev, base$modelFit$asDF$dev, tolerance = 1e-8)
})

test_that("scaling reaches the univariable ORs, interactions and the forest plot", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced", c("spont", "induced"))),
                covScales = list(list(var = "induced", scale = "hundred")),
                uniMulti = TRUE, forest = TRUE,
                refLevels = list(list(var = "case", ref = "0"), list(var = "spont", ref = "0")))
    tb <- coefDF(r)
    uni <- clogit(case == "1" ~ induced + strata(stratum), data = d)
    expect_equal(tb$or_u[4], unname(exp(100 * coef(uni))), tolerance = 1e-8)
    expect_equal(tb$level[4:6], c("per 100 units", "1 × per 100 units", "2 × per 100 units"))
    st <- r$models[[1]]$forest$state
    expect_true("induced (per 100 units)" %in% st$label)
})

test_that("no covScales entry means per unit", {
    d <- infertData()
    r <- runClr(d, covs = "induced", blocks = list(list("induced")))
    expect_equal(coefDF(r)$level, "per unit")
})
