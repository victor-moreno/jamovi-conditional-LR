# Cloned from jmv's logRegBin (R/logregbin.b.R). The changes: a matching
# variable defines the matched sets, the models are conditional logistic
# regressions (Cox models stratified by matched set, R/clr_core.R) instead of
# glm(), there is no intercept, and the parts that need absolute probabilities
# (marginal means, prediction, saved values) are gone.
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
        CICoefEst = function() {
            if (is.null(private$.CICoefEst))
                private$.CICoefEst <- private$.computeCICoefEst()

            return(private$.CICoefEst)
        },
        CICoefEstOR = function() {
            if (is.null(private$.CICoefEstOR))
                private$.CICoefEstOR <- private$.computeCICoefEst(type="OR")

            return(private$.CICoefEstOR)
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
        .CICoefEst = NULL,
        .CICoefEstOR = NULL,
        .rowNamesModel = NULL,
        .levelsDep = NULL,
        .refLevels = NULL,

        #### Init + run functions ----
        .init = function() {
            if (! is.null(self$options$dep))
                private$.errorCheck()

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
        .computeCICoefEst = function(type="LOR") {
            if (type == "OR")
                level <- self$options$ciWidthOR/100
            else
                level <- self$options$ciWidth/100

            ci <- list()
            for (i in seq_along(self$models)) {
                ci[[i]] <- stats::confint(self$models[[i]], level=level)
                if (type == "OR")
                    ci[[i]] <- exp(ci[[i]])
            }
            return(ci)
        },

        #### Init tables/plots functions ----
        .initModelFitTable = function() {
            table <- self$results$modelFit

            for (i in seq_along(private$.getModelTerms()))
                table$addRow(rowKey=i, values=list(model = i))

            table$setNote(
                "n",
                jmvcore::format(
                    .("Models estimated using N={n} observations in {sets} matched sets"),
                    n="...", sets="..."
                )
            )
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
            factors <- self$options$factors
            termsAll <- private$.getModelTerms()

            dep <- self$options$dep
            if (is.null(dep) ) {
                note <- paste0(.("Estimates represent the log odds of \u2026 within matched sets"))
            } else {
                depLevels <- private$.getLevelsDep()
                note <- jmvcore::format(
                    .('Estimates represent the log odds of "{dep} = {level1}" vs. "{dep} = {level2}", conditional on the matched sets'),
                    dep = dep,
                    level1 = depLevels$other,
                    level2 = depLevels$ref
                )
            }

            rowNamesModel <- private$.getRowNamesModel()

            ciWidthTitleString <- .('{ciWidth}% Confidence Interval')

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$coef

                ciWidth <- self$options$ciWidth
                ciWidthTitle <- jmvcore::format(ciWidthTitleString, ciWidth=ciWidth)
                table$getColumn('lower')$setSuperTitle(ciWidthTitle)
                table$getColumn('upper')$setSuperTitle(ciWidthTitle)

                ciWidthOR <- self$options$ciWidthOR
                ciWidthORTitle <- jmvcore::format(ciWidthTitleString, ciWidth=ciWidthOR)
                table$getColumn('oddsLower')$setSuperTitle(ciWidthORTitle)
                table$getColumn('oddsUpper')$setSuperTitle(ciWidthORTitle)

                coefTerms <- rowNamesModel[[i]]

                terms <- termsAll[[i]]

                # no intercept in a conditional model: coefTerms starts at the
                # first predictor (jmv's iter started at 1 for the intercept)
                iter <- 0
                for (j in seq_along(terms)) {
                    if (any(terms[[j]] %in% factors)) {

                        table$addRow(rowKey=terms[[j]],
                                     values=list(term = paste0(jmvcore::stringifyTerm(terms[[j]]), ':'),
                                                 est='', se='', odds='', z='', p='',
                                                 lower='', upper='', oddsLower='', oddsUpper=''))

                        contrastNames <- private$.contrastsCoefTable(terms[[j]])

                        for (k in seq_along(contrastNames)) {
                            iter <- iter + 1
                            rowKey <- jmvcore::composeTerm(coefTerms[[iter]])
                            table$addRow(rowKey=rowKey, values=list(term = contrastNames[[k]]))
                            table$addFormat(rowKey=rowKey, col=1, jmvcore::Cell.INDENTED)
                        }
                    } else {
                        iter <- iter + 1
                        rowKey <- jmvcore::composeTerm(jmvcore::toB64(terms[[j]]))
                        table$addRow(rowKey=rowKey, values=list(term = jmvcore::stringifyTerm(terms[[j]])))
                    }
                }

                table$setNote("est", note)
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

            nInf <- sum(comp$nSets[comp$informative])
            table$setNote(
                "sets",
                jmvcore::format(
                    .("{sets} matched sets, {inf} informative. Sets without a case or without a control do not contribute to the conditional likelihood."),
                    sets=sum(comp$nSets), inf=nInf
                )
            )
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

            table$setNote(
                "n",
                jmvcore::format(
                    .("Models estimated using N={n} observations in {sets} matched sets"),
                    n=self$models[[1]]$n,
                    sets=length(unique(self$dataProcessed$.STRATA))
                )
            )
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
            termsAll <- private$.getRowNamesModel()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$coef

                model <- summary(self$models[[i]])
                CI <- self$CICoefEst[[i]]
                CIOR <- self$CICoefEstOR[[i]]
                coef<- model$coefficients
                terms <- termsAll[[i]]
                rowTerms <- jmvcore::decomposeTerms(rownames(coef))

                for (j in seq_along(terms)) {
                    term <- terms[[j]]

                    # check which rows have the same length + same terms
                    index <- which(length(term) == sapply(rowTerms, length) &
                                       sapply(rowTerms, function(x) all(term %in% x)))

                    # coxph keeps aliased coefficients as NA rows, e.g. for a
                    # predictor that is constant within every matched set
                    if (length(index) != 1 || is.na(coef[index, 'coef'])) {
                        table$setNote(
                            "singular",
                            .("Not all coefficients could be estimated: a predictor is constant within the matched sets (e.g. a matching factor) or collinear with others")
                        )
                        next
                    }

                    row <- list()
                    row[["est"]] <- coef[index, 'coef']
                    row[["se"]] <- coef[index, 'se(coef)']
                    row[["odds"]] <- exp(coef[index, 'coef'])
                    row[["z"]] <- coef[index, 'z']
                    row[["p"]] <- coef[index, 'Pr(>|z|)']
                    row[["lower"]] <- CI[index, 1]
                    row[["upper"]] <- CI[index, 2]
                    row[["oddsLower"]] <- CIOR[index, 1]
                    row[["oddsUpper"]] <- CIOR[index, 2]

                    table$setRow(rowKey=jmvcore::composeTerm(term), values = row)
                }

                # e.g. "Loglik converged before variable 1; coefficient may be
                # infinite": separation within the matched sets
                warnings <- self$models[[i]]$clrWarnings
                if (length(warnings) > 0)
                    table$setNote("fit", paste(.("Model fit warning:"), paste(warnings, collapse="; ")))
            }
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
                    table$setNote("gvif", .("For terms with more than one coefficient the generalised VIF is shown"))
            }
        },

        #### Helper functions ----
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
        .getRowNamesModel = function() {
            if (is.null(private$.rowNamesModel)) {
                factors <- self$options$factors
                termsAll <- private$.getModelTerms()

                rowNamesAll <- list()
                for (i in seq_along(termsAll)) {
                    rowNames <- list()

                    terms <- termsAll[[i]]
                    for (term in terms) {
                        if (any(term %in% factors))
                            rowNames <- c(rowNames, private$.contrastModel(term))
                        else
                            rowNames[[length(rowNames) + 1]] <- jmvcore::toB64(term)
                    }
                    rowNamesAll[[i]] <- rowNames
                }
                private$.rowNamesModel <- rowNamesAll
            }

            return(private$.rowNamesModel)
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
        },
        .contrastModel = function(terms) {
            # Returns the names of the factor contrasts as they are
            # defined by the coxph function

            factors <- self$options$factors

            contrast <- list()
            for (term in terms) {
                if (term %in% factors) {
                    contrast[[term]] <- paste0(
                        jmvcore::toB64(term),
                        jmvcore::toB64(private$.getLevels(term)$other)
                    )
                } else {
                    contrast[[term]] <- jmvcore::toB64(term)
                }
            }

            contrastGrid <- expand.grid(contrast)
            contrastTerms <- list()
            for (i in 1:nrow(contrastGrid))
                contrastTerms[[i]] <- as.character(unlist(contrastGrid[i,]))

            return(contrastTerms)
        },
        .contrastsCoefTable = function(terms) {
            # Returns the names of the factor contrasts as they are
            # displayed in the jmv coef table

            factors <- self$options$factors
            refLevels <- self$refLevels
            refVars <- sapply(refLevels, function(x) x$var)

            levels <- list()
            for (factor in factors)
                levels[[factor]] <- levels(self$data[[factor]])

            contrast <- list()
            for (term in terms) {
                if (term %in% factors) {
                    ref <- refLevels[[which(term == refVars)]][['ref']]
                    refNo <- which(ref == levels[[term]])

                    contrLevels <- levels[[term]][-refNo]
                    refLevel <- levels[[term]][refNo]

                    c <- paste(contrLevels, refLevel, sep = ' \u2013 ')
                    if (length(terms) > 1)
                        c <- paste0('(', c, ')')

                    contrast[[term]] <- c
                } else {
                    contrast[[term]] <- term
                }
            }
            contrastGrid <- expand.grid(contrast)
            contrastNames <- apply(contrastGrid, 1, jmvcore::stringifyTerm)

            return(contrastNames)
        })
)
