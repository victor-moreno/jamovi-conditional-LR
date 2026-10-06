library(survival)


test_that("matched-set table counts the composition of the sets", {
    d <- infertData()
    r <- runClr(d, covs = "induced", blocks = list(list("induced")))
    sets <- r$sets$asDF
    expect_equal(sets$cases, c(1L, 1L))
    expect_equal(sets$controls, c(1L, 2L))
    expect_equal(sets$nSets, c(1L, 82L))
    expect_equal(sum(sets$nObs), nrow(d))
})

test_that("uninformative sets and missing values change nothing but the counts", {
    d <- infertData()
    extra <- d[d$stratum == "1" & d$case == "0", ]
    extra$stratum <- "999"
    d2 <- rbind(d, extra)
    d2$stratum <- factor(d2$stratum)
    d2$induced[5] <- NA
    r <- runClr(d2, covs = "induced", blocks = list(list("induced")))
    cl <- survival::clogit(case == "1" ~ induced + survival::strata(stratum), data = d2)
    expect_equal(coefDF(r)$est, unname(coef(cl)), tolerance = 1e-8)
    sets <- r$sets$asDF
    expect_true(any(sets$cases == 0))
    expect_equal(sum(sets$nObs), nrow(d2) - 1)
})

test_that("a matching factor gets no coefficient and no test", {
    d <- infertData()
    d$education <- factor(d$education)
    r <- runClr(d, covs = "induced", factors = "education",
                blocks = list(list("induced", "education")), omni = TRUE)
    tb <- coefDF(r)
    expect_true(all(is.na(suppressWarnings(as.numeric(tb$est[3:4])))))
    expect_true(is.na(suppressWarnings(as.numeric(r$models[[1]]$lrt$asDF$chi[2]))))
})

test_that("no informative set is an error, more than two levels too", {
    d <- infertData()
    d$stratum <- factor(ifelse(d$case == "1", "a", "b"))
    expect_error(runClr(d, covs = "induced", blocks = list(list("induced"))), "No matched set")

    d <- infertData()
    d$case <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
    expect_error(runClr(d, covs = "induced", blocks = list(list("induced"))), "more than two levels")
})

test_that("results serialise as jamovi sends them", {
    skip_if_not_installed("RProtoBuf")
    d <- infertData()
    opts <- condLogRegOptions$new(dep = "case", strata = "stratum", covs = "induced",
                                  factors = "spont", blocks = list(list("induced", "spont")),
                                  omni = TRUE, collin = TRUE)
    an <- condLogRegClass$new(options = opts, data = d, datasetId = "d",
                              analysisId = "1", revision = 1)
    an$init()
    an$run()
    jmvcore:::initProtoBuf()
    raw <- RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)
    expect_gt(length(raw), 0)
})
