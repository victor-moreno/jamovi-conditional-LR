infertData <- function() {
    d <- datasets::infert
    d$case <- factor(d$case)
    d$stratum <- factor(d$stratum)
    d$spont <- factor(d$spontaneous)
    d
}

runClr <- function(data, ...) {
    condLogReg(data = data, dep = "case", strata = "stratum", ...)
}

coefDF <- function(results, model = 1) {
    results$models[[model]]$coef$asDF
}
