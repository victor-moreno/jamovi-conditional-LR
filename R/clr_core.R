# Pure compute for conditional logistic regression: no jmvcore, no R6, no file
# access, so every function here is testable on its own (tests/testthat).
#
# The conditional likelihood of a matched case-control study is the partial
# likelihood of a Cox model stratified by matched set, with every subject given
# the same "time" and cases as the events (Breslow et al. 1978); this is what
# survival::clogit() does, with method = "exact" by default.

# Fit one model. `rhs` is the right-hand side (terms only, may be empty);
# `data` holds the response columns .TIME/.EVENT and the set id .STRATA.
# Warnings (e.g. infinite estimates under separation) are returned, not raised.
clrFit <- function(rhs, data, ties = "exact") {
    rhs <- if (length(rhs) == 0 || rhs == "") "" else paste(rhs, "+")
    formula <- stats::as.formula(
        paste("survival::Surv(.TIME, .EVENT) ~", rhs, "survival::strata(.STRATA)"))
    warnings <- character()
    fit <- withCallingHandlers(
        survival::coxph(formula, data = data, method = ties, model = TRUE),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    fit$clrWarnings <- unique(warnings)
    fit
}

# Log-likelihoods of the null (beta = 0) and fitted model, and the number of
# estimated (non-aliased) coefficients. A model with no terms has one loglik.
clrLogLik <- function(fit) {
    ll <- fit$loglik
    k <- if (is.null(fit$coefficients)) 0L else sum(!is.na(fit$coefficients))
    list(ll0 = ll[1], ll1 = ll[length(ll)], k = k)
}

# Fit measures, mirroring jmv's logRegBin: deviance = -2 logLik, AIC, BIC with
# the number of cases as the sample size (as stats::BIC does for coxph), the
# overall likelihood ratio test against beta = 0, and McFadden, Cox & Snell and
# Nagelkerke pseudo-R² with n = observations analysed.
clrFitMeasures <- function(fit) {
    l <- clrLogLik(fit)
    n <- fit$n
    dev <- -2 * l$ll1
    chi <- 2 * (l$ll1 - l$ll0)
    r2cs <- 1 - exp(-chi / n)
    list(
        dev = dev,
        aic = dev + 2 * l$k,
        bic = dev + l$k * log(fit$nevent),
        r2mf = 1 - l$ll1 / l$ll0,
        r2cs = r2cs,
        r2n = r2cs / (1 - exp(2 * l$ll0 / n)),
        chi = chi,
        df = l$k,
        p = stats::pchisq(chi, l$k, lower.tail = FALSE))
}

# Likelihood ratio test between two nested models fitted to the same rows.
clrCompare <- function(fit1, fit2) {
    l1 <- clrLogLik(fit1)
    l2 <- clrLogLik(fit2)
    chi <- 2 * (l2$ll1 - l1$ll1)
    df <- l2$k - l1$k
    list(chi = chi, df = df, p = stats::pchisq(chi, df, lower.tail = FALSE))
}

# The model's terms without the strata, with their position among all term
# labels: the model matrix's "assign" numbers terms counting strata(), which R
# orders with the main effects, i.e. before any interaction.
clrTerms <- function(fit) {
    all <- attr(stats::terms(fit), "term.labels")
    index <- which(! grepl("strata(", all, fixed = TRUE))
    list(labels = all[index], index = index)
}

# Omnibus (type III) likelihood ratio test of each term: refit with that term's
# columns removed from the model matrix, the coding of the others unchanged.
# car::Anova cannot be used here: it substitutes Wald tests for any model with
# strata, and stops on aliased coefficients (e.g. a matching factor).
clrTermTests <- function(fit, ties = "exact") {
    X <- stats::model.matrix(fit)
    assign <- attr(X, "assign")
    terms <- clrTerms(fit)
    labels <- terms$labels
    mf <- stats::model.frame(fit)
    y <- mf[[1]]
    sets <- mf[[grep("strata(", names(mf), fixed = TRUE)]]
    est <- ! is.na(fit$coefficients)
    llFull <- clrLogLik(fit)$ll1

    out <- data.frame(term = labels, chi = NA_real_, df = NA_integer_, p = NA_real_)
    for (j in seq_along(labels)) {
        drop <- assign == terms$index[j]
        df <- sum(drop & est)
        if (df == 0)
            next
        keep <- ! drop & est
        if (any(keep)) {
            Xr <- X[, keep, drop = FALSE]
            red <- survival::coxph(y ~ Xr + survival::strata(sets), method = ties)
            llRed <- clrLogLik(red)$ll1
        } else {
            llRed <- clrLogLik(fit)$ll0
        }
        out$chi[j] <- 2 * (llFull - llRed)
        out$df[j] <- df
        out$p[j] <- stats::pchisq(out$chi[j], df, lower.tail = FALSE)
    }
    out
}

# Generalised variance inflation factor of each term (Fox & Monette 1992),
# from the correlation matrix of the estimated coefficients; the same formula
# as car::vif, without counting the strata as a term. For a one-column term it
# is the ordinary VIF.
clrVIF <- function(fit) {
    est <- ! is.na(fit$coefficients)
    assign <- attr(stats::model.matrix(fit), "assign")[est]
    terms <- clrTerms(fit)
    labels <- terms$labels
    R <- stats::cov2cor(as.matrix(stats::vcov(fit))[est, est, drop = FALSE])
    detR <- det(R)
    vif <- vapply(seq_along(labels), function(j) {
        idx <- which(assign == terms$index[j])
        if (length(idx) == 0)
            return(NA_real_)
        if (length(idx) == ncol(R))
            return(1)
        det(R[idx, idx, drop = FALSE]) * det(R[-idx, -idx, drop = FALSE]) / detR
    }, numeric(1))
    stats::setNames(vif, labels)
}

# Composition of the matched sets: one row per (cases, controls) pattern with
# the number of sets and observations. Sets with no case or no control add
# nothing to the conditional likelihood (uninformative).
clrSetComposition <- function(event, strata) {
    cases <- tapply(event, strata, sum)
    size <- tapply(event, strata, length)
    controls <- size - cases
    key <- paste(cases, controls)
    tab <- data.frame(cases = as.integer(cases), controls = as.integer(controls))
    tab <- tab[! duplicated(key), ]
    tab <- tab[order(tab$cases, tab$controls), ]
    k <- paste(tab$cases, tab$controls)
    tab$nSets <- as.integer(table(key)[k])
    tab$nObs <- as.integer(tab$nSets * (tab$cases + tab$controls))
    tab$informative <- tab$cases > 0 & tab$controls > 0
    rownames(tab) <- NULL
    tab
}
