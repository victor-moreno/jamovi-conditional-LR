# Forest plot of odds ratios, adapted from the prognosis module's Cox forest
# plot (R/coxutils.R there): one line per table row, reference levels as an
# open diamond at 1, the estimates and their intervals as text on the right,
# and univariable beside multivariable estimates when both are shown.
# `r` has one row per plot line: label, ref (reference level), or, lo, hi, p
# and, for univariable estimates, or_u, lo_u, hi_u.

clrFmtOR <- function(or, lo, hi, p, ref, . = identity) {
    ifelse(ref, .("Reference"),
    ifelse(is.na(or), "",
           sprintf("%.2f (%.2f-%.2f)   %s", or, lo, hi,
                   ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))))
}

clrForestPlot <- function(r, uniMulti = FALSE, col = "#1B4F8A", . = identity) {
    header <- if (uniMulti) .("OR uni | multi") else sprintf("%-19s %s", .("OR (95% CI)"), "p")
    text <- if (uniMulti)
        ifelse(r$ref, .("Reference"), ifelse(is.na(r$or), "",
               paste0(ifelse(is.na(r$or_u), "  -  ", sprintf("%5.2f", r$or_u)),
                      " | ", sprintf("%.2f", r$or))))
    else clrFmtOR(r$or, r$lo, r$hi, r$p, r$ref, . = .)
    n <- nrow(r)
    y <- rev(seq_len(n))
    lab <- data.frame(y = c(n + 1, y), label = c("", r$label), text = c(header, text))

    uLab <- .("Univariable"); mLab <- .("Multivariable")
    if (uniMulti) {
        m <- data.frame(y = y - 0.15, or = r$or, lo = r$lo, hi = r$hi, model = mLab, ref = r$ref)
        u <- data.frame(y = y + 0.15, or = r$or_u, lo = r$lo_u, hi = r$hi_u, model = uLab, ref = r$ref)
        pts <- rbind(u, m)
        pts$model <- factor(pts$model, c(uLab, mLab))
        cols <- stats::setNames(c("grey55", col), c(uLab, mLab))
    } else {
        pts <- data.frame(y = y, or = r$or, lo = r$lo, hi = r$hi, model = mLab, ref = r$ref)
        cols <- stats::setNames(col, mLab)
    }
    refs <- data.frame(y = y[r$ref], or = rep(1, sum(r$ref)))
    est <- pts[! pts$ref & ! is.na(pts$or), ]

    p <- ggplot2::ggplot(est, ggplot2::aes(.data$or, .data$y, colour = .data$model)) +
        ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
        ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                                linewidth = 0.8, na.rm = TRUE) +
        ggplot2::geom_point(shape = 15, size = 3, na.rm = TRUE) +
        ggplot2::geom_point(data = refs, ggplot2::aes(.data$or, .data$y), inherit.aes = FALSE,
                            shape = 23, size = 2.5, fill = "white", colour = col) +
        ggplot2::scale_colour_manual(values = cols, name = NULL) +
        ggplot2::scale_x_log10(labels = function(x) formatC(x, format = "fg")) +
        ggplot2::scale_y_continuous(breaks = lab$y, labels = lab$label, limits = c(0.4, n + 1.4),
                                    expand = ggplot2::expansion(add = 0),
                                    sec.axis = ggplot2::dup_axis(labels = lab$text, name = NULL)) +
        ggplot2::labs(x = .("Odds ratio (log scale)"), y = NULL) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.text.y.left = ggplot2::element_text(hjust = 0, colour = "black"),
                       axis.text.y.right = ggplot2::element_text(hjust = 0, colour = "black",
                                                                 family = "mono", size = 11))
    p + if (uniMulti) ggplot2::theme(legend.position = "top")
        else ggplot2::theme(legend.position = "none")
}
