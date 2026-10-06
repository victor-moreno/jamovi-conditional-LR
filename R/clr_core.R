# Pure compute for conditional logistic regression: no jmvcore, no R6, no file
# access, so every function here is testable on its own (tests/testthat).
#
# The conditional likelihood of a matched case-control study is the partial
# likelihood of a Cox model stratified by matched set, with every subject given
# the same "time" and cases as the events (Breslow et al. 1978); this is what
# survival::clogit() does, with method = "exact" by default.

# Fit one model. `rhs` is the right-hand side (terms only, may be empty);
# `data` holds the response columns .TIME/.EVENT and the set id .STRATA.
#
# The model matrix is built here and its aliased columns are found by
# clrAliased() before fitting, instead of leaving that to coxph(): coxph's
# check is relative to the information matrix itself, so a predictor that is
# constant within every set (zero information up to rounding) was dropped on
# macOS but estimated on Linux (jamovi in docker) as a meaningless coefficient
# with SE ~1e6 when it was the only predictor. Columns are coded as coxph would
# (treatment contrasts, intercept removed). Warnings (e.g. infinite estimates
# under separation) are returned, not raised.
clrFit <- function(rhs, data, ties = "exact") {
    y <- survival::Surv(data$.TIME, data$.EVENT)
    sets <- data$.STRATA
    if (length(rhs) == 0 || rhs == "") {
        X <- matrix(numeric(0), nrow = nrow(data), ncol = 0)
        assign <- integer(0)
        labels <- character(0)
    } else {
        tt <- stats::terms(stats::as.formula(paste("~", rhs)))
        mm <- stats::model.matrix(tt, data)
        X <- mm[, -1, drop = FALSE]
        assign <- attr(mm, "assign")[-1]
        labels <- attr(tt, "term.labels")
    }
    aliased <- clrAliased(X, sets)

    warnings <- character()
    cox <- withCallingHandlers(
        clrCox(X[, ! aliased, drop = FALSE], y, sets, ties),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })

    coef <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
    var <- matrix(NA_real_, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X)))
    if (any(! aliased)) {
        coef[! aliased] <- cox$coefficients
        var[! aliased, ! aliased] <- cox$var
    }
    ll <- cox$loglik
    list(coefficients = coef, var = var, loglik = c(ll[1], ll[length(ll)]),
         n = cox$n, nevent = cox$nevent, X = X, assign = assign, labels = labels,
         aliased = aliased, y = y, sets = sets, ties = ties,
         clrWarnings = unique(warnings))
}

# The stratified Cox fit on a numeric matrix (possibly with no columns, the
# null model, whose loglik is a single value).
clrCox <- function(X, y, sets, ties) {
    if (ncol(X) == 0)
        return(survival::coxph(y ~ survival::strata(sets), method = ties))
    survival::coxph(y ~ X + survival::strata(sets), method = ties)
}

# Columns of the model matrix that the conditional likelihood cannot estimate:
# no variation within any matched set (e.g. a matching variable), or a linear
# combination of other columns after removing the set means. Judged on the
# within-set deviations scaled by the column's overall spread, so the test does
# not depend on rounding or on the units of the variable.
clrAliased <- function(X, sets, tol = 1e-7) {
    aliased <- rep(FALSE, ncol(X))
    if (ncol(X) == 0)
        return(aliased)
    within <- apply(X, 2, function(x) x - stats::ave(x, sets))
    within <- matrix(within, nrow = nrow(X))
    spread <- sqrt(colSums(sweep(X, 2, colMeans(X))^2))
    ratio <- sqrt(colSums(within^2)) / spread
    aliased <- unname(! is.finite(ratio) | ratio < tol)
    ok <- which(! aliased)
    if (length(ok) > 1) {
        q <- qr(sweep(within[, ok, drop = FALSE], 2, spread[ok], "/"), tol = tol)
        if (q$rank < length(ok))
            aliased[ok[q$pivot[-seq_len(q$rank)]]] <- TRUE
    }
    aliased
}

# Estimates, standard errors, Wald z and p for every column (NA if aliased).
clrCoefTable <- function(fit) {
    est <- fit$coefficients
    se <- sqrt(diag(fit$var))
    z <- est / se
    data.frame(est = est, se = se, z = z, p = 2 * stats::pnorm(-abs(z)),
               row.names = names(est))
}

# Wald confidence intervals, as confint() gives for a coxph/clogit fit.
clrConfint <- function(fit, level = 0.95) {
    tab <- clrCoefTable(fit)
    q <- stats::qnorm(1 - (1 - level) / 2)
    cbind(lower = tab$est - q * tab$se, upper = tab$est + q * tab$se)
}

# The terms whose columns are (at least partly) aliased.
clrAliasedTerms <- function(fit) {
    fit$labels[sort(unique(fit$assign[fit$aliased]))]
}

# Log-likelihoods of the null (beta = 0) and fitted model, and the number of
# estimated (non-aliased) coefficients. A model with no terms has one loglik.
clrLogLik <- function(fit) {
    list(ll0 = fit$loglik[1], ll1 = fit$loglik[2], k = sum(! is.na(fit$coefficients)))
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

# Omnibus (type III) likelihood ratio test of each term: refit with that term's
# columns removed from the model matrix, the coding of the others unchanged.
# car::Anova cannot be used here: it substitutes Wald tests for any model with
# strata, and stops on aliased coefficients (e.g. a matching factor).
clrTermTests <- function(fit, ties = fit$ties) {
    est <- ! is.na(fit$coefficients)
    llFull <- fit$loglik[2]

    labels <- fit$labels
    out <- data.frame(term = labels, chi = NA_real_, df = NA_integer_, p = NA_real_)
    for (j in seq_along(labels)) {
        drop <- fit$assign == j
        df <- sum(drop & est)
        if (df == 0)
            next
        keep <- ! drop & est
        red <- clrCox(fit$X[, keep, drop = FALSE], fit$y, fit$sets, ties)
        llRed <- red$loglik[length(red$loglik)]
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
    assign <- fit$assign[est]
    labels <- fit$labels
    R <- stats::cov2cor(fit$var[est, est, drop = FALSE])
    detR <- det(R)
    vif <- vapply(seq_along(labels), function(j) {
        idx <- which(assign == j)
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
