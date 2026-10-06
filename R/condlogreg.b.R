# Cloned from jmv's logRegBin (R/logregbin.b.R). The changes: a matching
# variable defines the matched sets, the models are conditional logistic
# regressions (Cox models stratified by matched set, R/clr_core.R) instead of
# glm(), there is no intercept, and the parts that need absolute probabilities
# (marginal means, prediction, saved values) are gone. The coefficient table,
# the univariable/multivariable option, the trend test, the forest plot and
# the notes switch follow the prognosis module's Cox regression instead of jmv.
#' @importFrom jmvcore .
condLogRegClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "condLogRegClass",
    inherit = condLogRegBase,
    #### Active bindings ----
    active = list(
        dataProcessed = function() {
            if (is.null(private$.dataProcessed))
                private$.dataProcessed <- private$.cleanData()

            return(private$.dataProcessed)
        },
        formulas = function() {
            if (is.null(private$.formulas))
                private$.formulas <- private$.getFormulas()

            return(private$.formulas)
        },
        models = function() {
            if (is.null(private$.models))
                private$.models <- private$.computeModels()

            return(private$.models)
        },
        nModels = function() {
            if (is.null(private$.nModels))
                private$.nModels <- length(self$options$blocks)

            return(private$.nModels)
        },
        fitMeasures = function() {
            if (is.null(private$.fitMeasures))
                private$.fitMeasures <- lapply(self$models, clrFitMeasures)

            return(private$.fitMeasures)
        },
        lrtModelTerms = function() {
            if (is.null(private$.lrtModelTerms))
                private$.lrtModelTerms <- lapply(
                    self$models, clrTermTests, ties = self$options$ties)

            return(private$.lrtModelTerms)
        },
        VIF = function() {
            if (is.null(private$.VIF))
                private$.VIF <- lapply(self$models, clrVIF)

            return(private$.VIF)
        },
        refLevels = function() {
            if (is.null(private$.refLevels)) {
                factors <- c(self$options$dep, self$options$factors)
                refLevels <- getReferenceLevels(
                    self$data, factors, self$options$refLevels
                )
                private$.refLevels <- refLevels$refLevels

                if (length(refLevels$changedVars) > 0)
                    setRefLevelWarning(self, refLevels$changedVars)
            }

            return(private$.refLevels)
        }
    ),
    private = list(
        #### Member variables ----
        .dataProcessed = NULL,
        .formulas = NULL,
        .models = NULL,
        .nModels = NULL,
        .modelTerms = NULL,
        .fitMeasures = NULL,
        .lrtModelTerms = NULL,
        .VIF = NULL,
        .levelsDep = NULL,
        .refLevels = NULL,
        .uniFits = list(),
        .coefValues = list(),
        .blankCoef = list(or_u='', lower_u='', upper_u='', p_u='', ptrend_u='',
                          or='', lower='', upper='', p='', ptrend='', est='', se='', z=''),

        #### Init + run functions ----
        .init = function() {
            if (! is.null(self$options$dep))
                private$.errorCheck()

            private$.initDescTable()
            private$.initModelFitTable()
            private$.initModelCompTable()
            private$.initModelSpec()
            private$.initLrtTables()
            private$.initCoefTables()
            private$.initCollinearityTable()
        },
        .run = function() {
            if (is.null(self$options$dep) || is.null(self$options$strata))
                return()

            private$.errorCheck()
            private$.checkSets()
            private$.populateSetsTable()

            if (self$nModels < 1 || length(self$options$blocks[[1]]) == 0)
                return()

            private$.populateModelFitTable()
            private$.populateModelCompTable()
            private$.populateLrtTables()
            private$.populateCoefTables()
            private$.populateForest()
            private$.populateCollinearityTable()
        },

        #### Compute results ----
        .computeModels = function() {
            data <- self$dataProcessed

            models <- list()
            for (i in seq_along(self$formulas))
                models[[i]] <- clrFit(self$formulas[[i]], data, self$options$ties)

            return(models)
        },
        # univariable fit of one variable, on the rows of the models; cached
        # across the models of the blocks
        .uniFit = function(var) {
            if (is.null(private$.uniFits[[var]]))
                private$.uniFits[[var]] <- clrFit(
                    jmvcore::toB64(var), self$dataProcessed, self$options$ties)
            private$.uniFits[[var]]
        },
        # level order (1, 2, 3, ...) of a factor on the rows of the models, as
        # the user ordered the levels (dataProcessed has the reference first)
        .trendScores = function(factor) {
            rows <- match(rownames(self$dataProcessed), rownames(self$data))
            x <- self$data[[factor]]
            match(as.character(x), levels(x))[rows]
        },

        #### Init tables/plots functions ----
        .initDescTable = function() {
            if ( ! self$options$descTable || is.null(self$options$dep))
                return()

            # added to jmv's logRegBin. The rows depend on the values (a
            # Missing row only for a variable with missing values), and
            # .init() gets only the data header from jamovi, so the dataset
            # is read here (as the prognosis module does); the table is
            # filled here as well, all in one place
            if (nrow(self$data) == 0)
                private$.data <- self$readDataset()
            data <- self$data
            table <- self$results$desc
            dep <- self$options$dep
            strata <- self$options$strata

            depLevels <- private$.getLevelsDep()
            table$getColumn('ctrl')$setTitle(depLevels$ref)
            table$getColumn('cases')$setTitle(depLevels$other)
            table$getColumn('ctrl')$setSuperTitle(dep)
            table$getColumn('cases')$setSuperTitle(dep)

            # rows with the case/control status and the matched set known;
            # convert whole columns before subsetting, which would drop the
            # "values" attribute toNumeric needs for a nominal integer column
            keep <- ! is.na(data[[dep]])
            if ( ! is.null(strata))
                keep <- keep & ! is.na(data[[strata]])
            event <- as.integer(as.character(data[[dep]]) != depLevels$ref)[keep]

            covs <- list()
            for (cov in self$options$covs)
                covs[[cov]] <- jmvcore::toNumeric(data[[cov]])[keep]
            factors <- list()
            for (factor in self$options$factors)
                factors[[factor]] <- factor(data[[factor]], levels=levels(data[[factor]]))[keep]

            table$addRow(rowKey='.dep', values=c(list(var=dep, level=.('n (%)')),
                                                  clrDepFrequencies(event)))

            # never an empty level cell: jamovi renders these columns as row
            # headers and merges an empty one into the cell above, which
            # shifts the rest of the row
            desc <- clrDescriptives(event, covs, factors)
            for (k in seq_len(nrow(desc))) {
                level <- switch(desc$type[k],
                                cov=.('Mean \u00b1 SD'), missing=.('Missing'), desc$level[k])
                table$addRow(rowKey=desc$key[k], values=list(
                    var=desc$var[k], level=level,
                    ctrl=desc$ctrl[k], cases=desc$cases[k], total=desc$total[k]))
            }

            private$.note(table, "desc", paste(
                .("Mean \u00b1 SD for covariates; n (%) for factor levels, with percentages of the known values within each column; for the dependent variable and missing values, of the column total."),
                jmvcore::format(
                    .("All {n} rows with the dependent and matching variables known; the models use only the rows complete for every variable."),
                    n=sum(keep))))
        },
        .initModelFitTable = function() {
            table <- self$results$modelFit

            for (i in seq_along(private$.getModelTerms()))
                table$addRow(rowKey=i, values=list(model = i))
        },
        .initModelCompTable = function() {
            table <- self$results$modelComp
            terms <- private$.getModelTerms()

            if (length(terms) <= 1) {
                table$setVisible(visible = FALSE)
                return()
            }

            for (i in 1:(length(terms)-1))
                table$addRow(rowKey=i, values=list(model1 = i, model2 = as.integer(i+1)))

        },
        .initModelSpec = function() {
            groups <- self$results$models

            for (i in seq_along(private$.getModelTerms())) {
                groups$addItem(key=i)
                group <- groups$get(key=i)
                group$setTitle(paste(.("Model"),i))
            }
        },
        .initLrtTables = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$lrt
                terms <- termsAll[[i]]

                for (j in seq_along(terms))
                    table$addRow(rowKey=paste0(terms[[j]]), values=list(term = jmvcore::stringifyTerm(terms[j])))
            }
        },
        .initCoefTables = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            ciTitle <- jmvcore::format(.('{ciWidth}% CI'), ciWidth=self$options$ciWidth)

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$coef
                if (self$options$uniMulti) {
                    for (col in c('or', 'lower', 'upper', 'p', 'ptrend'))
                        table$getColumn(col)$setSuperTitle(.('Multivariable'))
                } else {
                    table$getColumn('lower')$setSuperTitle(ciTitle)
                    table$getColumn('upper')$setSuperTitle(ciTitle)
                }

                layout <- private$.coefLayout(i)
                for (row in layout)
                    table$addRow(rowKey=row$key, values=list(var=row$var, level=row$level))

                # forest plot: one line per table row, plus a header line per
                # factor and per interaction
                headers <- sum(vapply(termsAll[[i]], function(t)
                    length(t) > 1 || t %in% self$options$factors, logical(1)))
                groups$get(key=i)$forest$setSize(700, 90 + 26 * (length(layout) + headers))
            }
        },
        .initCollinearityTable = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {

                table <- groups$get(key=i)$assump$collin
                terms <- termsAll[[i]]

                if (length(terms) < 1)
                    terms <- ''

                for (i in seq_along(terms))
                    table$addRow(rowKey=i, values=list(term = jmvcore::stringifyTerm(terms[i])))
            }
        },

        #### Populate tables functions ----
        .populateSetsTable = function() {
            if ( ! self$options$setSummary)
                return()

            # rows depend on the data, which .init() does not see (header only)
            table <- self$results$sets
            data <- self$dataProcessed
            comp <- clrSetComposition(data$.EVENT, data$.STRATA)

            for (i in seq_len(nrow(comp))) {
                table$addRow(rowKey=i, values=list(
                    cases=comp$cases[i], controls=comp$controls[i],
                    nSets=comp$nSets[i], nObs=comp$nObs[i]))
            }

            nSets <- sum(comp$nSets)
            nInf <- sum(comp$nSets[comp$informative])
            private$.note(table, "sets", jmvcore::format(
                .("{sets} matched sets, {inf} informative. Sets without a case or without a control do not contribute to the conditional likelihood."),
                sets=nSets, inf=nInf))
            # a warning, shown even with the notes off
            table$setNote("uninf", if (nInf < nSets && ! self$options$showNotes)
                jmvcore::format(.("{n} uninformative sets (no case or no control left)"),
                                n=nSets - nInf))
        },
        .populateModelFitTable = function() {
            table <- self$results$modelFit
            for (i in 1:self$nModels) {
                fm <- self$fitMeasures[[i]]
                row <- list()
                row[["r2mf"]] <- fm$r2mf
                row[["r2cs"]] <- fm$r2cs
                row[["r2n"]] <- fm$r2n
                row[["dev"]] <- fm$dev
                row[["aic"]] <- fm$aic
                row[["bic"]] <- fm$bic
                row[["chi"]] <- fm$chi
                row[["df"]] <- fm$df
                row[["p"]] <- fm$p

                table$setRow(rowNo=i, values = row)
            }

            private$.note(table, "n", jmvcore::format(
                .("Models estimated using N={n} observations in {sets} matched sets"),
                n=self$models[[1]]$n,
                sets=length(unique(self$dataProcessed$.STRATA))))
        },
        .populateModelCompTable = function() {
            if (self$nModels <= 1)
                return()

            table <- self$results$modelComp

            for (i in 1:(self$nModels - 1)) {
                r <- clrCompare(self$models[[i]], self$models[[i + 1]])
                table$setRow(rowNo=i, values = list(chi=r$chi, df=r$df, p=r$p))
            }
        },
        .populateLrtTables = function() {
            if ( ! self$options$omni)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            lrTests <- self$lrtModelTerms

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$lrt

                terms <- termsAll[[i]]
                termsB64 <- lapply(terms, jmvcore::toB64)
                lrt <- lrTests[[i]]
                rowTerms <- jmvcore::decomposeTerms(lrt$term)

                for (j in seq_along(terms)) {
                    term <- termsB64[[j]]

                    # check which rows have the same length + same terms
                    index <- which(length(term) == sapply(rowTerms, length) &
                                       sapply(rowTerms, function(x) all(term %in% x)))

                    row <- list()
                    row[["chi"]] <- lrt$chi[index]
                    row[["df"]] <- lrt$df[index]
                    row[["p"]] <- lrt$p[index]

                    table$setRow(rowKey=paste0(terms[[j]]), values = row)

                    if (is.na(lrt$chi[index]))
                        table$setNote(
                            "singular",
                            .("A term with no estimable coefficient (e.g. a matching factor) is not tested")
                        )
                }
            }
        },
        .populateCoefTables = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            factors <- self$options$factors
            level <- self$options$ciWidth / 100

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$coef
                fit <- self$models[[i]]
                coef <- clrCoefTable(fit)
                ci <- exp(clrConfint(fit, level=level))
                termsB64 <- jmvcore::composeTerms(lapply(termsAll[[i]], jmvcore::toB64))
                values <- list()

                for (row in private$.coefLayout(i)) {
                    term <- termsAll[[i]][[row$term]]
                    main <- length(term) == 1
                    isFactor <- main && term %in% factors
                    nLevels <- if (isFactor) length(levels(self$data[[term]])) else 0
                    vals <- list()

                    if (row$ref) {
                        # reference level: OR 1; for a factor with 3+ levels
                        # the p of the whole factor and of the trend
                        vals$or <- 1
                        if (nLevels >= 3) {
                            lrt <- self$lrtModelTerms[[i]]
                            vals$p <- lrt$p[lrt$term == jmvcore::toB64(term)]
                            if (self$options$trend)
                                vals$ptrend <- clrTrendP(termsB64, jmvcore::toB64(term),
                                                         private$.trendScores(term),
                                                         self$dataProcessed, self$options$ties)
                        }
                        if (self$options$uniMulti) {
                            vals$or_u <- 1
                            if (nLevels >= 3) {
                                vals$p_u <- clrFitMeasures(private$.uniFit(term))$p
                                if (self$options$trend)
                                    vals$ptrend_u <- clrTrendP(jmvcore::toB64(term), jmvcore::toB64(term),
                                                               private$.trendScores(term),
                                                               self$dataProcessed, self$options$ties)
                            }
                        }
                    } else {
                        # aliased columns (clrAliased) stay NA and their cells empty
                        k <- private$.findCoef(row$parts, rownames(coef))
                        if (length(k) == 1 && ! is.na(coef$est[k])) {
                            vals$or <- exp(coef$est[k])
                            vals$lower <- ci[k, 1]
                            vals$upper <- ci[k, 2]
                            vals$p <- coef$p[k]
                            vals$est <- coef$est[k]
                            vals$se <- coef$se[k]
                            vals$z <- coef$z[k]
                        }
                        if (self$options$uniMulti && main) {
                            uni <- private$.uniFit(term)
                            uCoef <- clrCoefTable(uni)
                            uCi <- exp(clrConfint(uni, level=level))
                            k <- private$.findCoef(row$parts, rownames(uCoef))
                            if (length(k) == 1 && ! is.na(uCoef$est[k])) {
                                vals$or_u <- exp(uCoef$est[k])
                                vals$lower_u <- uCi[k, 1]
                                vals$upper_u <- uCi[k, 2]
                                vals$p_u <- uCoef$p[k]
                            }
                        }
                    }

                    vals <- Filter(function(v) length(v) == 1 && ! is.na(v), vals)
                    # '' rather than nothing: jamovi shows an unset number as '.'
                    table$setRow(rowKey=row$key, values=utils::modifyList(private$.blankCoef, vals))
                    values[[row$key]] <- vals
                }
                private$.coefValues[[i]] <- values

                private$.note(table, "est", private$.coefNote(i))

                # name the terms that cannot be estimated, rather than only
                # leaving their rows empty (a warning: always shown)
                aliased <- clrAliasedTerms(fit)
                table$setNote("singular", if (length(aliased) > 0)
                    jmvcore::format(
                        .("Not estimable: {terms}. A predictor that does not vary within the matched sets (e.g. a matching variable) or that is collinear with others cannot be estimated by conditional logistic regression."),
                        terms=private$.termNames(aliased)))

                # e.g. "Loglik converged before variable 1; coefficient may be
                # infinite": separation within the matched sets
                warnings <- fit$clrWarnings
                table$setNote("fit", if (length(warnings) > 0)
                    paste(.("Model fit warning:"), paste(warnings, collapse="; ")))
            }
        },
        .populateForest = function() {
            if ( ! self$options$forest)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                image <- groups$get(key=i)$forest
                if ( ! is.null(image$state))
                    next

                # the rows of the coefficient table, with a header line for
                # each factor and interaction; levels indented under it
                values <- private$.coefValues[[i]]
                get <- function(vals, name) if (is.null(vals[[name]])) NA_real_ else vals[[name]]
                lines <- list()
                shown <- character()
                for (row in private$.coefLayout(i)) {
                    term <- termsAll[[i]][[row$term]]
                    header <- length(term) > 1 || term %in% self$options$factors
                    if (header && ! row$var %in% shown) {
                        lines[[length(lines) + 1]] <- data.frame(
                            label=row$var, ref=FALSE, or=NA_real_, lo=NA_real_, hi=NA_real_,
                            p=NA_real_, or_u=NA_real_, lo_u=NA_real_, hi_u=NA_real_)
                        shown <- c(shown, row$var)
                    }
                    v <- values[[row$key]]
                    lines[[length(lines) + 1]] <- data.frame(
                        label=if (header) paste0("    ", row$level) else row$var,
                        ref=row$ref, or=if (row$ref) 1 else get(v, "or"),
                        lo=get(v, "lower"), hi=get(v, "upper"), p=get(v, "p"),
                        or_u=get(v, "or_u"), lo_u=get(v, "lower_u"), hi_u=get(v, "upper_u"))
                }
                image$setState(do.call(rbind, lines))
            }
        },
        .forestPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            print(clrForestPlot(image$state, uniMulti=self$options$uniMulti, .=self$translate))
            TRUE
        },
        .populateCollinearityTable = function() {
            if ( ! self$options$collin)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$assump$collin
                terms <- lapply(termsAll[[i]], jmvcore::toB64)
                VIF <- self$VIF[[i]]
                rowTerms <- jmvcore::decomposeTerms(names(VIF))

                for (j in seq_along(terms)) {
                    # check which rows have the same length + same terms
                    index <- which(length(terms[[j]]) == sapply(rowTerms, length) &
                                       sapply(rowTerms, function(x) all(terms[[j]] %in% x)))

                    row <- list()
                    row[["vif"]] <- as.numeric(VIF[index])
                    row[["tol"]] <- 1 / as.numeric(VIF[index])
                    table$setRow(rowNo=j, values=row)
                }

                if (any(vapply(termsAll[[i]], function(t) any(t %in% self$options$factors), logical(1))))
                    private$.note(table, "gvif", .("For terms with more than one coefficient the generalised VIF is shown"))
            }
        },

        #### Helper functions ----
        # an explanatory note: shown only with "Show notes" (warnings call
        # setNote directly and are always shown)
        .note = function(table, key, text) {
            table$setNote(key, if (self$options$showNotes) text)
        },
        .coefNote = function(i) {
            dep <- self$options$dep
            depLevels <- private$.getLevelsDep()
            terms <- private$.getModelTerms()[[i]]
            hasFactor <- any(vapply(terms, function(t) length(t) == 1 && t %in% self$options$factors, logical(1)))
            hasCov <- any(vapply(terms, function(t) length(t) == 1 && t %in% self$options$covs, logical(1)))
            # one literal per sentence, so each is a whole msgid for translators
            paste(c(
                jmvcore::format(
                    .('Odds ratios (OR) of "{dep} = {level1}" vs. "{dep} = {level2}" within the matched sets, with {ciWidth}% confidence intervals.'),
                    dep=dep, level1=depLevels$other, level2=depLevels$ref,
                    ciWidth=self$options$ciWidth),
                if (hasFactor)
                    .("Factors: OR against the reference level (OR = 1); on the reference row of a factor with three or more levels, p of the likelihood ratio test of the whole factor."),
                if (hasCov)
                    .("Covariates: OR per one-unit increase."),
                if (length(terms) > 1)
                    .("Multivariable ORs are adjusted for the other variables in the model."),
                if (self$options$uniMulti)
                    .("Univariable: each variable alone, on the same rows as the model."),
                if (self$options$trend)
                    .("p trend: the levels scored 1, 2, 3, ... in their order (Wald test); the multivariable trend keeps the other variables of the model."),
                if (any(lengths(terms) > 1))
                    .("With interactions, main-effect ORs apply at the reference level (or 0) of the interacting variable.")),
                collapse=" ")
        },
        # Rows of the coefficient table of model i, from the options and the
        # level names alone (so they exist in .init()): for a factor, its
        # reference level (OR 1) and one row per other level; for a covariate,
        # one row; for an interaction, one row per coefficient. `parts` are
        # the model-matrix names of the variables in a coefficient.
        .coefLayout = function(i) {
            terms <- private$.getModelTerms()[[i]]
            factors <- self$options$factors
            rows <- list()
            for (j in seq_along(terms)) {
                t <- terms[[j]]
                var <- paste(t, collapse=" \u00d7 ")
                isf <- t %in% factors
                if (length(t) == 1 && isf)
                    rows[[length(rows) + 1]] <- list(key=paste0('.ref', t), term=j, var=var,
                                                     level=private$.getLevels(t)$ref,
                                                     ref=TRUE, parts=NULL)
                # non-reference levels crossed in model-matrix order (the
                # first variable varies fastest, as expand.grid)
                grid <- lapply(t, function(v) if (v %in% factors) private$.getLevels(v)$other else NA)
                g <- expand.grid(grid, stringsAsFactors=FALSE)
                for (r in seq_len(nrow(g))) {
                    lv <- unlist(g[r, ], use.names=FALSE)
                    parts <- ifelse(isf, paste0(jmvcore::toB64(t), jmvcore::toB64(lv)), jmvcore::toB64(t))
                    level <- if (any(isf)) paste(lv[isf], collapse=" \u00d7 ") else .("per unit")
                    rows[[length(rows) + 1]] <- list(key=paste(parts, collapse=':'), term=j,
                                                     var=var, level=level, ref=FALSE, parts=parts)
                }
            }
            rows
        },
        # the coefficient whose model-matrix name holds the same parts (R may
        # order the variables of an interaction differently)
        .findCoef = function(parts, names) {
            which(vapply(strsplit(names, ':', fixed=TRUE), function(x)
                length(x) == length(parts) && all(parts %in% x), logical(1)))
        },
        .termNames = function(labels) {
            # model-matrix term labels (base64 names, ':' for interactions)
            # back to the names the user sees
            names <- vapply(labels, function(l) {
                jmvcore::stringifyTerm(jmvcore::fromB64(jmvcore::decomposeTerm(l)))
            }, '')
            paste(names, collapse=', ')
        },
        .getModelTerms = function() {
            if (is.null(private$.modelTerms)) {
                blocks <- self$options$blocks

                terms <- list()
                if (is.null(blocks)) {
                    terms[[1]] <- c(self$options$covs, self$options$factors)
                } else {
                    for (i in seq_along(blocks)) {
                        terms[[i]] <- unlist(blocks[1:i], recursive = FALSE)
                    }
                }
                private$.modelTerms <- terms
            }

            return(private$.modelTerms)
        },
        .getLevelsDep = function() {
            if (is.null(private$.levelsDep) && ! is.null(self$options$dep))
                private$.levelsDep <- private$.getLevels(self$options$dep)

            return(private$.levelsDep)
        },
        .getLevels = function(var) {
            levels <- levels(self$data[[var]])
            refLevels <- self$refLevels
            refVars <- sapply(refLevels, function(x) x$var)
            ref <- refLevels[[which(var == refVars)]][['ref']]
            # workaround ... client currently initially sends null
            # before issuing a second request with the correct value
            if (is.null(ref))
                ref <- levels[1]
            other <- levels[-which(ref == levels)]

            return(list("ref"=ref, "other"=other))
        },
        .getFormulas = function() {
            terms <- private$.getModelTerms()

            # right-hand sides only: clrFit() adds the response and the strata
            formulas <- list()
            for (i in seq_along(terms)) {
                termsB64 <- lapply(terms[[i]], jmvcore::toB64)
                composedTerms <- jmvcore::composeTerms(termsB64)
                formulas[[i]] <- paste0(composedTerms, collapse ="+")
            }

            return(formulas)
        },
        .errorCheck = function() {
            dep <- self$options$dep
            column <- self$data[[dep]]

            if (length(levels(column)) > 2) {
                jmvcore::reject(
                    jmvcore::format(
                        .("The dependent variable '{}' has more than two levels; conditional logistic regression can only be performed on dependent variables with two levels (case and control)."),
                        dep
                    ),
                    code=''
                )
            }
        },
        .checkSets = function() {
            data <- self$dataProcessed
            comp <- clrSetComposition(data$.EVENT, data$.STRATA)
            if ( ! any(comp$informative)) {
                jmvcore::reject(
                    .("No matched set contains both a case and a control, so the conditional likelihood is empty. Check the matching variable."),
                    code=''
                )
            }
        },
        .cleanData = function(naOmit=TRUE) {
            dep <- self$options$dep
            strata <- self$options$strata
            covs <- self$options$covs
            factors <- self$options$factors
            refLevels <- self$refLevels

            dataRaw <- self$data

            data <- list()

            refVars <- sapply(refLevels, function(x) x$var)

            for (factor in c(dep, factors)) {
                ref <- refLevels[[which(factor == refVars)]][['ref']]
                column <- factor(
                    dataRaw[[factor]],
                    ordered = FALSE,
                    levels = levels(dataRaw[[factor]])
                )
                levels(column) <- jmvcore::toB64(levels(column))
                column <- stats::relevel(column, ref = jmvcore::toB64(ref))

                data[[jmvcore::toB64(factor)]] <- column
            }

            for (cov in covs)
                data[[jmvcore::toB64(cov)]] <- jmvcore::toNumeric(dataRaw[[cov]])

            # the response of the stratified Cox model: one common "time" and
            # the cases (non-reference level of dep) as events
            depB64 <- data[[jmvcore::toB64(dep)]]
            data[[".EVENT"]] <- ifelse(is.na(depB64), NA, as.integer(as.integer(depB64) != 1))
            data[[".TIME"]] <- rep(1, length(depB64))
            data[[".STRATA"]] <- as.character(dataRaw[[strata]])

            attr(data, 'row.names') <- rownames(self$data)
            attr(data, 'class') <- 'data.frame'

            if (naOmit)
                data <- data[stats::complete.cases(data), , drop = FALSE]

            data[[".STRATA"]] <- factor(data[[".STRATA"]])

            return(data)
        })
)
