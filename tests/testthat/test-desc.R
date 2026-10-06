library(survival)


test_that("descriptive table: dependent variable, covariates and factor levels", {
    d <- infertData()
    d$education <- factor(d$education)
    r <- runClr(d, covs = "age", factors = "education",
                blocks = list(list("education")), descTable = TRUE)
    tb <- r$desc$asDF
    expect_equal(tb$var[1], "case")
    expect_equal(tb$ctrl[1], sprintf("%d (%.1f%%)", sum(d$case == "0"), 100 * mean(d$case == "0")))
    expect_equal(tb$total[1], sprintf("%d (100.0%%)", nrow(d)))
    expect_equal(tb$cases[2], sprintf("%.2f ± %.2f", mean(d$age[d$case == "1"]), sd(d$age[d$case == "1"])))
    lev <- levels(d$education)
    expect_equal(tb$level[3:5], lev)
    n <- table(d$education[d$case == "1"])[lev]
    expect_equal(tb$cases[3:5], sprintf("%d (%.1f%%)", as.integer(n), 100 * as.vector(n) / sum(d$case == "1")))
})

test_that("descriptive table counts only the analysed (complete) rows", {
    d <- infertData()
    d$induced[1:5] <- NA
    r <- runClr(d, covs = c("induced", "spontaneous"),
                blocks = list(list("spontaneous")), descTable = TRUE)
    expect_equal(r$desc$asDF$total[1], sprintf("%d (100.0%%)", nrow(d) - 5))
})

test_that("odds-ratio-only table gives OR, CI and p, p last", {
    d <- infertData()
    r <- runClr(d, covs = c("spontaneous", "induced"),
                blocks = list(list("spontaneous", "induced")), orTable = TRUE)
    tbl <- r$models[[1]]$coef
    visible <- vapply(tbl$columns, function(c) isTRUE(c$visible), logical(1))
    expect_equal(names(tbl$columns)[visible], c("term", "odds", "oddsLower", "oddsUpper", "pOR"))
    cl <- survival::clogit(case == "1" ~ spontaneous + induced + survival::strata(stratum),
                           data = transform(d, one = 1))
    tb <- tbl$asDF
    expect_equal(tb$odds, unname(exp(coef(cl))), tolerance = 1e-8)
    expect_equal(tb$pOR, tb$p)

    full <- runClr(d, covs = "induced", blocks = list(list("induced")))$models[[1]]$coef
    visible <- vapply(full$columns, function(c) isTRUE(c$visible), logical(1))
    expect_equal(names(full$columns)[visible], c("term", "est", "se", "z", "p"))
})

test_that("the bundled data sets load and fit", {
    files <- list.files(test_path("..", "..", "data"), pattern = "\\.csv$", full.names = TRUE)
    if (length(files) == 0)
        files <- list.files(system.file("data", package = "condLogReg"), pattern = "\\.csv$", full.names = TRUE)
    expect_setequal(basename(files), c("infert.csv", "bdendo.csv", "nested.csv"))
    bd <- read.csv(files[basename(files) == "bdendo.csv"], stringsAsFactors = TRUE, na.strings = "")
    bd$case <- factor(bd$case)
    r <- condLogReg(data = bd, dep = "case", strata = "set", factors = c("est", "gall"),
                    blocks = list(list("est", "gall")), OR = TRUE)
    cl <- survival::clogit(case == "1" ~ est + gall + survival::strata(set), data = bd)
    expect_equal(r$models[[1]]$coef$asDF$odds[c(2, 4)], unname(exp(coef(cl))), tolerance = 1e-8)
})

test_that("a nominal integer covariate (factor with values, as jamovi sends it) works", {
    d <- infertData()
    d$induced[3] <- NA
    lev <- sort(unique(d$age))
    ageF <- factor(d$age, levels = lev)
    attr(ageF, "values") <- as.integer(lev)
    d$ageF <- ageF
    r <- runClr(d, covs = c("ageF", "induced"), blocks = list(list("induced")), descTable = TRUE)
    keep <- ! is.na(d$induced)
    expect_equal(r$desc$asDF$total[2], sprintf("%.2f ± %.2f", mean(d$age[keep]), sd(d$age[keep])))
})

test_that("descriptive table with factors only, or covariates only", {
    d <- infertData()
    r <- runClr(d, factors = "spont", blocks = list(list("spont")), descTable = TRUE)
    expect_equal(nrow(r$desc$asDF), 1 + nlevels(d$spont))
    r <- runClr(d, covs = "induced", blocks = list(list("induced")), descTable = TRUE)
    expect_equal(nrow(r$desc$asDF), 2)
})
