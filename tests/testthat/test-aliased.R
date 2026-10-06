library(survival)

test_that("a matching variable alone gives no estimate and a named note", {
    d <- infertData()
    r <- runClr(d, covs = "age", blocks = list(list("age")), modelTest = TRUE)
    tb <- coefDF(r)
    expect_true(is.na(suppressWarnings(as.numeric(tb$est[1]))))
    notes <- vapply(r$models[[1]]$coef$notes, function(n) n$note, "")
    expect_true(any(grepl("Not estimable: age", notes, fixed = TRUE)))
    expect_false(any(grepl("warning", notes)))
    fit <- r$modelFit$asDF
    expect_equal(fit$aic, fit$dev)
    expect_equal(fit$df, 0)
    expect_equal(fit$r2mf, 0)
})

test_that("aliasing does not depend on rounding noise or the units", {
    set.seed(3)
    sets <- factor(rep(1:30, each = 3))
    level <- rnorm(30, 40, 8)[sets]
    X <- cbind(age = level, ageNoise = level + rnorm(90, 0, 1e-12), ageDays = level * 365.25,
               x = rnorm(90))
    expect_equal(clrAliased(X, sets), c(TRUE, TRUE, TRUE, FALSE))
})

test_that("a predictor collinear with another only within sets is aliased", {
    set.seed(4)
    sets <- factor(rep(1:40, each = 3))
    x1 <- rnorm(120)
    x2 <- 2 * x1 + rnorm(40)[sets]
    expect_equal(clrAliased(cbind(x1, x2, x3 = rnorm(120)), sets), c(FALSE, TRUE, FALSE))
})

test_that("interaction with a matching variable is estimated as clogit does", {
    d <- infertData()
    r <- runClr(d, covs = c("induced", "age"),
                blocks = list(list("induced", "age", c("induced", "age"))))
    cl <- clogit(case == "1" ~ induced + induced:age + strata(stratum), data = d)
    tb <- coefDF(r)
    est <- suppressWarnings(as.numeric(tb$est))
    expect_equal(est[c(1, 3)], unname(coef(cl)), tolerance = 1e-8)
    expect_true(is.na(est[2]))
})
