# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container: coefficients against
# the clogit fit (exact stratified Cox) on the infert matched case-control data.
d <- datasets::infert
d$case <- factor(d$case)
r <- condLogReg::condLogReg(
    data = d, dep = "case", strata = "stratum",
    covs = c("spontaneous", "induced"),
    blocks = list(list("spontaneous", "induced"))
)
tbl <- r$models[[1]]$coef$asDF
d$one <- 1
oracle <- coef(survival::coxph(
    survival::Surv(one, case == "1") ~ spontaneous + induced +
        survival::strata(stratum), data = d, method = "exact"))
stopifnot(isTRUE(all.equal(tbl$est, unname(oracle), tolerance = 1e-6)))
cat(sprintf("   smoke test passed: spontaneous = %.4f, induced = %.4f\n",
            tbl$est[1], tbl$est[2]))
