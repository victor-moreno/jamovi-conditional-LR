library(survival)

refs <- function(...) c(list(list(var = "case", ref = "0")), list(...))

test_that("prognosis layout: reference row OR 1, level-only labels, default columns", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced")),
                refLevels = refs(list(var = "spont", ref = "0")))
    tbl <- r$models[[1]]$coef
    tb <- tbl$asDF
    expect_equal(tb$var, c("spont", "spont", "spont", "induced"))
    expect_equal(tb$level, c("0", "1", "2", "per unit"))
    expect_equal(tb$or[1], 1)
    visible <- vapply(tbl$columns, function(c) isTRUE(c$visible), logical(1))
    expect_equal(names(tbl$columns)[visible], c("var", "level", "or", "lower", "upper", "p"))

    r2 <- runClr(d, covs = "induced", blocks = list(list("induced")),
                 uniMulti = TRUE, trend = TRUE, coefDetails = TRUE)
    tbl <- r2$models[[1]]$coef
    visible <- vapply(tbl$columns, function(c) isTRUE(c$visible), logical(1))
    expect_equal(names(tbl$columns)[visible],
                 c("var", "level", "or_u", "lower_u", "upper_u", "p_u", "ptrend_u",
                   "or", "lower", "upper", "p", "ptrend", "est", "se", "z"))
})

test_that("reference row p: LR test of the whole factor (3+ levels), uni and multi", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced")), uniMulti = TRUE,
                refLevels = refs(list(var = "spont", ref = "0")))
    tb <- coefDF(r)
    full <- clogit(case == "1" ~ spont + induced + strata(stratum), data = d)
    red <- clogit(case == "1" ~ induced + strata(stratum), data = d)
    expect_equal(tb$p[1], pchisq(2 * (full$loglik[2] - red$loglik[2]), 2, lower.tail = FALSE),
                 tolerance = 1e-8)
    uni <- clogit(case == "1" ~ spont + strata(stratum), data = d)
    expect_equal(tb$p_u[1], unname(summary(uni)$logtest["pvalue"]), tolerance = 1e-8)
    expect_equal(tb$or_u[2:3], unname(exp(coef(uni))), tolerance = 1e-8)
})

test_that("univariable ORs use the rows of the model", {
    d <- infertData()
    d$induced[c(2, 10, 40)] <- NA
    r <- runClr(d, covs = c("spontaneous", "induced"),
                blocks = list(list("spontaneous", "induced")), uniMulti = TRUE)
    cc <- d[! is.na(d$induced), ]
    uni <- clogit(case == "1" ~ spontaneous + strata(stratum), data = cc)
    expect_equal(coefDF(r)$or_u[1], unname(exp(coef(uni))), tolerance = 1e-8)
    expect_equal(coefDF(r)$lower_u[1], unname(exp(confint(uni)[1, 1])), tolerance = 1e-8)
})

test_that("trend test: level order as scores, whatever the reference level", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced")), trend = TRUE, uniMulti = TRUE,
                refLevels = refs(list(var = "spont", ref = "2")))
    tb <- coefDF(r)
    expect_equal(tb$level[1], "2")
    d$score <- as.integer(d$spont)
    multi <- clogit(case == "1" ~ score + induced + strata(stratum), data = d)
    uni <- clogit(case == "1" ~ score + strata(stratum), data = d)
    expect_equal(tb$ptrend[1], summary(multi)$coefficients["score", "Pr(>|z|)"], tolerance = 1e-8)
    expect_equal(tb$ptrend_u[1], summary(uni)$coefficients["score", "Pr(>|z|)"], tolerance = 1e-8)
})

test_that("no trend or whole-factor p for a two-level factor", {
    d <- infertData()
    d$two <- factor(ifelse(d$spontaneous > 0, "yes", "no"))
    r <- runClr(d, factors = "two", blocks = list(list("two")), trend = TRUE)
    tb <- coefDF(r)
    expect_true(is.na(tb$p[1]))
    expect_true(is.na(tb$ptrend[1]))
})

test_that("interaction rows are labelled by their levels", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced", c("spont", "induced"))),
                refLevels = refs(list(var = "spont", ref = "0")))
    tb <- coefDF(r)
    expect_equal(tb$var[5:6], rep("spont × induced", 2))
    expect_equal(tb$level[5:6], c("1", "2"))
    cl <- clogit(case == "1" ~ spont * induced + strata(stratum), data = d)
    expect_equal(tb$or[5:6], unname(exp(coef(cl)[4:5])), tolerance = 1e-8)
})

test_that("notes: explanations only with showNotes, warnings always", {
    d <- infertData()
    notes <- function(tbl) vapply(tbl$notes, function(n) n$note, "")
    r <- runClr(d, covs = c("induced", "age"), blocks = list(list("induced", "age")))
    expect_equal(names(r$models[[1]]$coef$notes), "singular")
    expect_length(r$modelFit$notes, 0)
    expect_length(r$sets$notes, 0)
    r <- runClr(d, covs = c("induced", "age"), blocks = list(list("induced", "age")),
                showNotes = TRUE)
    expect_setequal(names(r$models[[1]]$coef$notes), c("est", "singular"))
    expect_true(grepl("within the matched sets", notes(r$models[[1]]$coef)[["est"]]))
    expect_length(r$modelFit$notes, 1)

    d$induced[d$case == "1" & d$stratum == "1"] <- NA
    r <- runClr(d, covs = "induced", blocks = list(list("induced")))
    expect_equal(names(r$sets$notes), "uninf")
})

test_that("forest plot state: one line per table row plus headers, same numbers", {
    d <- infertData()
    r <- runClr(d, covs = "induced", factors = "spont",
                blocks = list(list("spont", "induced")), forest = TRUE, uniMulti = TRUE,
                refLevels = refs(list(var = "spont", ref = "0")))
    st <- r$models[[1]]$forest$state
    tb <- coefDF(r)
    expect_equal(nrow(st), nrow(tb) + 1)
    expect_equal(trimws(st$label), c("spont", "0", "1", "2", "induced"))
    expect_equal(st$ref, c(FALSE, TRUE, FALSE, FALSE, FALSE))
    expect_equal(st$or[3:5], tb$or[2:4])
    expect_equal(st$or_u[3:5], tb$or_u[2:4])
    p <- clrForestPlot(st, uniMulti = TRUE)
    expect_s3_class(p, "ggplot")
})
