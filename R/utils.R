# Helpers copied unchanged from jmv (R/utilsreg.R, R/utils.R), which
# logRegBin's reference-level handling relies on.

getReferenceLevels = function(data, vars, refLevels) {
    updatedRefLevels <- list()
    changedVars <- c()

    # Create a named list from the refLevels input for easier access
    refLevelsList <- stats::setNames(
        lapply(refLevels, function(ref) ref$ref),
        sapply(refLevels, function(ref) ref$var)
    )

    for (var in vars) {
        factorLevels <- levels(data[[var]])
        refLevel <- refLevelsList[[var]]

        # If no refLevel is provided or the provided level is invalid, use the first level
        if (is.null(refLevel) || ! (refLevel %in% factorLevels)) {
            refLevel <- factorLevels[1]
            changedVars <- c(changedVars, var)
        }

        updatedRefLevels[[ length(updatedRefLevels) + 1 ]] <- list(
            var = var, ref = refLevel
        )
    }

    return(list(refLevels=updatedRefLevels, changedVars=changedVars))
}

setRefLevelWarning = function(self, changedVars) {
    message <- jmvcore::format(
        .("The specified reference level was not found for the following variable(s): {vars}. Defaulting to the first available level. To use a custom reference level, ensure the defined reference level is available in the data."),
        vars=listItems(self, changedVars)
    )

    setAnalysisNotice(
        self, message, name="refLevelWarning", type=jmvcore::NoticeType$STRONG_WARNING
    )
}

listItems <- function(self, items, quote = '\'') {

    items <- vapply(items, function(x) paste0(quote, x, quote), '', USE.NAMES=FALSE)

    if (length(items) == 1) {

        l <- items

    } else if (length(items) == 2) {

        l <- jmvcore::format(.('{item1} and {item2}'), item1=items[1], item2=items[2])

    } else {

        l <- items[1]
        for (i in seq(2, length(items)-1))
            l <- jmvcore::format(.('{list}, {nextItem}'), list=l, nextItem=items[i])
        l <- jmvcore::format(.('{list}, and {lastItem}'), list=l, lastItem=items[length(items)])
    }

    return(l)
}

setAnalysisNotice = function(self, message, name, type) {
    notice <- jmvcore::Notice$new(
        options = self$options,
        name = 'warningMessage',
        type = jmvcore::NoticeType$STRONG_WARNING
    )
    notice$setContent(message)
    self$results$insert(1, notice)
}
