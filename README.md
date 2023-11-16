
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesQuantify

<!-- badges: start -->

[![R-CMD-check](https://github.com/liusihan/BayesQuantify/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/liusihan/BayesQuantify/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/liusihan/BayesQuantify/branch/master/graph/badge.svg)](https://app.codecov.io/gh/liusihan/BayesQuantify?branch=master)
<!-- badges: end -->

The ACMG/AMP guidelines have undergone continuous review and refinement
for different rules, genes, and diseases, driving optimization and
enhancing variant interpretation standards in genetic testing. In 2018,
the ClinGen Sequence Variant Interpretation Working Group has proposed a
Bayesian Classification Framework to model the ACMG/AMP guidelines. This
framework has successfully quantified the thresholds for applying PM5
and PP3/BP4. However, there are challenges for clinicians in utilizing
the Bayesian Classification Framework, as tools and software for
convincingly calculating the positive likelihood ratio are lacking.

The [BayesQuantify](#bayesquantify) R Package provide a comprehensive
functions to define the thresholds for each evidence strength level

## Installation

You can install the development version of BayesQuantify like so:

``` r
devtools::install_github("liusihan/BayesQuantify")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BayesQuantify)
data("VCI_data")
data("ClinVar2020_AJHG_Pejaver_data")
```
