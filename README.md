# condLogReg — conditional logistic regression for jamovi

A [jamovi](https://www.jamovi.org) module for **matched case-control studies**: conditional
logistic regression, which estimates odds ratios by comparing each case with its own controls.
It is jmv's *Binomial Logistic Regression* with a **Matching Variable** added, fitted by the
conditional likelihood (a Cox model stratified by matched set, as `survival::clogit`).

Menu: **Regression ▸ Logistic Regression ▸ Conditional (Matched)**. English, Spanish and Catalan.

## Features

- Model builder with blocks and interactions, reference levels, and **covariate scaling**
  (OR per unit, per 1 SD, per 10, 100 or 1000 units).
- **Matched sets** table, and a **descriptive table** by case/control status with missing values.
- **Odds ratios** with confidence intervals, one row per level (reference level OR = 1),
  univariable beside multivariable, **trend test**, and **forest plot**.
- Model fit (deviance, AIC, BIC, pseudo R²), model comparisons, likelihood ratio tests,
  collinearity.
- Exact conditional likelihood, or Efron's or Breslow's approximation.

Example data in jamovi's Data Library: infertility after abortion (`datasets::infert`), endometrial
cancer and oestrogens (Breslow & Day 1980), and a simulated nested case-control study.

## Tutorial

[victor-moreno.github.io/jamovi-conditional-LR](https://victor-moreno.github.io/jamovi-conditional-LR/):
the method, how each result is computed, and three use cases step by step.

## Install

Download the `.jmo` for your jamovi from the [Releases](../../releases) page and, in jamovi,
choose **Modules (⊕) ▸ Sideload**:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (R 4.6.0) | `condLogReg_<version>_current_R4.6.0_arm64.jmo` | `condLogReg_<version>_current_R4.6.0_x64.jmo` |
| **solid** (R 4.5.0) | `condLogReg_<version>_solid_R4.5.0_arm64.jmo` | `condLogReg_<version>_solid_R4.5.0_x64.jmo` |

The same file works on macOS, Windows and Linux (see Help ▸ About for your R version).

## Acknowledgements

The analysis is a clone of `logRegBin` from [jmv](https://github.com/jamovi/jmv) (R. Selker,
J. Love, D. Dropmann; GPL ≥ 2); the results layout follows the
[prognosis](https://github.com/victor-moreno/jamovi-prognosis) module. Model fitting uses
[survival](https://github.com/therneau/survival) (T. M. Therneau). The endometrial cancer data
come from the R package [Epi](https://cran.r-project.org/package=Epi) (GPL-2). Designed and built
together with Claude Code.

This module has been developed with support of the Instituto de Salud Carlos III (ISCIII),
“Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number
FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health
(CIBERESP), action Genrisk.

## Licence

GPL (≥ 2).
