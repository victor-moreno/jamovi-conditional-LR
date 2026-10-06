# condLogReg — conditional logistic regression for jamovi

A [jamovi](https://www.jamovi.org) module for **matched case-control studies** (and other
binary outcomes in matched sets or strata): conditional logistic regression, which estimates
odds ratios from the comparison of each case with its own controls.

It is a clone of jmv's *Binomial Logistic Regression* with one addition to the panel, a
**Matching Variable** that identifies the matched sets. The models are fitted by the
conditional likelihood, as a Cox model stratified by matched set (`survival::coxph(...,
method = "exact")`), the same fit as `survival::clogit`. It needs no packages beyond those
bundled with jamovi: `jmvcore`, `R6` and `survival`.

Appears in jamovi's menu as **Regression ▸ Logistic Regression ▸ Conditional (Matched)**.

## What it does

Kept from jmv's logistic regression, and working the same way:

- the **Model Builder** with blocks of predictors and interactions, and **Reference Levels**
  for the dependent variable (the reference level marks the controls) and the factors;
- **Model Fit**: deviance (−2 conditional log-likelihood), AIC, BIC, McFadden's, Cox & Snell's
  and Nagelkerke's pseudo R², and the overall likelihood ratio test;
- **Model Comparisons** between blocks (likelihood ratio tests);
- **Model Coefficients**: omnibus likelihood ratio tests per term, estimates (log odds ratios)
  and odds ratios with confidence intervals;
- **Collinearity statistics** (VIF, or the generalised VIF for factors).

Added:

- **Matching Variable**, and a **Matched Sets** table: the composition of the sets (cases and
  controls per set), with the uninformative sets (no case or no control left after removing
  missing values) counted in a note;
- **Likelihood**: the exact conditional likelihood (default), or Efron's or Breslow's
  approximation for sets with many cases. All three agree when each set has one case;
- a **Descriptive table** (*Descriptives* box): the frequencies of the dependent variable,
  mean ± SD of each covariate and n (%) of each factor level, in controls, cases and in
  total, for the rows that enter the models;
- **Odds ratios only** (*Model Coefficients*): the coefficient table reduced to the odds
  ratio, its confidence interval and p, p last, as in a published table.

Removed, because a conditional model has no intercept and so gives no absolute
probabilities: Tjur's R², estimated marginal means, the prediction section (classification
table, ROC curve, AUC, cut-off plot) and the saved predicted values, residuals and Cook's
distances. Data weights are ignored, since a weight would break the matched sets.

Two parts are computed by the module instead of `car`, which cannot handle the strata: the
omnibus tests are likelihood ratio tests from refitting without each term (`car::Anova`
substitutes Wald tests for any model with strata), and the GVIF does not count the strata as
a term.

A variable that is constant within every matched set (a matching factor) cancels from the
conditional likelihood. Its rows stay empty, and a note names it. Its interactions with
other predictors are estimable. The module finds such columns itself (no variation within
sets, or collinear with others after removing the set means) before fitting, because
`coxph`'s own check depends on rounding: on Linux it estimated a matching variable entered
alone as a meaningless coefficient.

Three data sets ship with the module, in jamovi's Data Library:

- **Infertility after abortion** (R's `datasets::infert`, Trichopoulos et al. 1976): 83 sets
  of a case and up to two controls;
- **Endometrial cancer and oestrogens** (Breslow & Day 1980, as `Epi::bdendo`): 63 sets of a
  case and four controls, with oestrogen use, its duration and dose, gallbladder disease,
  hypertension, obesity and other drugs;
- **Nested case-control study (simulated)**: 200 cases with up to four controls, smoking
  confounding the effect of coffee, missing values; its true odds ratios are given in the
  tutorial.

## Tutorial

A [website on conditional logistic regression](https://victor-moreno.github.io/jamovi-conditional-LR/)
(in `docs/`) covers why matched data need a conditional analysis (with a simulation of the
alternatives), how every number in the output is computed, and two use cases worked step by
step in jamovi, with their data: the `infert` study, the endometrial cancer study, and a
simulated 1:4 nested case-control study with confounding, a variable matching ratio and
missing values. Every table in it is
computed by this module. The site is rendered with Quarto from a source kept outside this
repository.

## Install

Prebuilt `.jmo` files are attached to the [Releases](../../releases) page — four per version, one
for each combination of jamovi series and CPU:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (bundles R 4.6.0) | `condLogReg_<version>_current_R4.6.0_arm64.jmo` | `condLogReg_<version>_current_R4.6.0_x64.jmo` |
| **solid** (bundles R 4.5.0) | `condLogReg_<version>_solid_R4.5.0_arm64.jmo` | `condLogReg_<version>_solid_R4.5.0_x64.jmo` |

The same file works on macOS, Windows and Linux: jamovi checks the R version and the CPU, not the
operating system (see Help → About). In jamovi, go to Modules (⊕) → Sideload and choose the
`.jmo`.

## Build and install (developers)

```sh
bash tools/install.sh                 # desktop and Docker (container "jamovi")
bash tools/install.sh desktop         # jamovi.app, via jmvtools::install()
bash tools/install.sh docker [name]   # jamovi >= 28.4 container, via ../jamovi-src's compiler
bash tools/build-jmo.sh current       # release: dist/*_current_R4.6.0_{x64,arm64}.jmo
bash tools/build-jmo.sh solid         # release: dist/*_solid_R4.5.0_{x64,arm64}.jmo
bash tools/release.sh                 # build all four and publish one GitHub release
```

After each install, `tools/smoke.R` checks the installed module against the `clogit` fit of
`infert`. The tests (`Rscript -e 'testthat::test_local(".")'`) check every table against
`survival::clogit`, nested-model `anova()` and `car::vif`. They also cover the reference
levels, ordinal factors, the three likelihoods, uninformative sets, a matching factor and
serialisation of the results as jamovi sends them.

## Acknowledgements

The analysis is a clone of `logRegBin` from [jmv](https://github.com/jamovi/jmv), by Ravi
Selker, Jonathon Love and Damian Dropmann (GPL ≥ 2). The panel, its JavaScript and most of the
R backend are theirs. Model fitting relies on [survival](https://github.com/therneau/survival)
(T. M. Therneau). The endometrial cancer data are Breslow & Day's, as distributed in the R
package [Epi](https://cran.r-project.org/package=Epi) (B. Carstensen, M. Plummer, E. Laara,
M. Hills; GPL-2).

This module was designed and built together with Claude Code.

Methods: Breslow et al. (1978) and Breslow & Day (1980) for the conditional analysis of
matched case-control studies; Gail, Lubin & Rubinstein (1981) for the exact likelihood with
several cases per set; Fox & Monette (1992) for the generalised VIF.

This module has been developed with support of the Instituto de Salud Carlos III (ISCIII),
“Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number
FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health
(CIBERESP), action Genrisk.

## Versions

- **0.2.0** — descriptive table, odds-ratios-only coefficient table, and two more data sets
  (endometrial cancer, simulated nested study).
- **0.1.1** — aliased predictors (e.g. a matching variable) are detected by the module, so
  they are reported the same on every platform; the note names them.
- **0.1.0** — first version.

## Licence

GPL (≥ 2), as jmv.
