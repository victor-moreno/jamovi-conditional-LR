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
    tb <- results$models[[model]]$coef$asDF
    # asDF holds the visible columns only: est (log OR) is hidden by default
    if (is.null(tb$est) && ! is.null(tb$or))
        tb$est <- suppressWarnings(log(as.numeric(tb$or)))
    tb
}
