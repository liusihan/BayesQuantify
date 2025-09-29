
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesQuantify <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/liusihan/BayesQuantify/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/liusihan/BayesQuantify/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/liusihan/BayesQuantify/branch/master/graph/badge.svg)](https://app.codecov.io/gh/liusihan/BayesQuantify?branch=master)
<!-- badges: end -->

## Overview
The ACMG/AMP guidelines have undergone continuous review and refinement
for different rules, genes, and diseases, driving optimization and
enhancing variant interpretation standards in genetic testing. In 2018,
the ClinGen Sequence Variant Interpretation Working Group has proposed a
Bayesian Classification Framework to model the ACMG/AMP guidelines. This
framework has successfully quantified the thresholds for applying PM5
and PP3/BP4. However, existing software and tools designed for quantifying 
the evidence strength and establishing corresponding thresholds to refine 
the ACMG/AMP criteria are lacking.

The `BayesQuantify` R Package provides users with a unified resource 
for quantifying the strength of evidence for ACMG/AMP criteria using a 
naive Bayes classifier. The functions included in BayesQuantify consist 
of five main steps (as shown in Figure 1):

![Figure 1](https://github.com/liusihan/BayesQuantify/blob/master/man/figures/Figure1.png?raw=true)
<p align="center"> Figure 1. Schematic overview of BayesQuantify. </p>

## Installation

You can install `BayesQuantify` from GitHub:

``` r
library(devtools)
devtools::install_github("liusihan/BayesQuantify")
```

## Required Input Data
The package requires one essential input: a data.frame containing variants classification results, including the variants, ACMG/AMP evidence, variant classification, and tested features.

Two distinct datasets: the ClinGen Curated Variants dataset and the ClinVar 2019 dataset have been included in `BayesQuantify`. Here is an example input file：

|#Variation|Assertion|Applied Evidence Codes (Met)|HGNC Gene Symbol|Expert Panel|...|
|---|---|---|---|---|---|
|NM_000277.2:c.1A>G|Pathogenic|PS3, PM3, PP4_Moderate, PM2|PAH|Phenylketonuria VCEP|...|
|NM_000314.6:c.209+3A>T|Uncertain Significance|PP3|PTEN|PTEN VCEP|...|
|NM_004999.3:c.2836C>T|Likely Benign|BS1|MYO6|Hearing Loss VCEP|...|
|...|...|...|...|...|...|


## Usage
The full usage of BayesQuantify can be found [here](https://github.com/liusihan/BayesQuantify/blob/79848f80b975984615efea4879d8aecaecb7115a/man/figures/BayesQuantify_1.0.0.pdf).

You can find a step-by-step example of how to use it below:
``` r
> library(BayesQuantify)

## Step 1: Calculation of OP and Post_P for each evidence strength with user-defined Prior_P.
# Prior_P=0.1 was recommended by ClinGen SVI Working Group
> auto_select_postp(0.1)
  Evidence Strength Odds of pathogenicity Posterior probability of pathogenicity and benignity
1               PVS                   351                                                0.975
2                PS      18.7349939951952                                     0.67550020016016
3                PM      4.32839392791312                                    0.324749849931156
4                PP      2.08047925438182                                    0.187760764369383
5                BP     0.480658481882885                                    0.949301150041276
6                BM     0.231032576205959                                    0.974972184931784
7                BS    0.0533760512683624                                    0.994104293142569
8               BVS   0.00284900284900285                                    0.999683544303797

## Step 2: Variants of US (VUS) subclassification
# This step will categorize VUS into six confidence tiers (hot, warm, tepid, cool, cold, and ice cold) according to
# the Association for Clinical Genomic Science (ACGS) Best Practice Guidelines.
> data("ClinGen_dataset")
> data <- add_info(ClinGen_dataset, "Assertion")
> data <- VUS_classify(data, "Assertion", "Applied Evidence Codes (Met)")

## Step3: Calculating LR or lr
# To ensure a robust and reliable truth set for evidence calibration, we recommend excluding hot, warm, tepid VUS.
# For binary variables: 
> data<-discrete_cutoff(data, "Applied Evidence Codes (Met)", criteria = "PM2_Supporting")
> data<-discrete_cutoff(data, "Applied Evidence Codes (Met)", criteria = "PM2")
> library(dplyr)
> truth_set <- data %>% filter(VUS_class %in% c("IceCold","Cold","Cool","")) %>% filter(!grepl("BA1",`Applied Evidence Codes (Met)`)) 
> LR_result<-LR(truth_set, 28, 29)

# For quantitative variables: 
> data("ClinVar_2019_dataset")
> data <- add_info(ClinVar_2019_dataset, "clnsig")
> local_bootstrapped_lr(data, "PrimateAI_score", "Pathogenic",0.0441, 10000, 100, 0.01, "test_dir")
> postp_list <- c(0.100, 0.211, 0.608, 0.981)
> lr_CI_result <- lr_CI(10000, "test_dir")
> get_lr_threshold(lr_CI_result, postp_list, "Pathogenic")

## Step 4: Variant classification
# BayesQuantify offers two functions for variant classification
> data("ClinGen_dataset")
> ACMG_Classification(ClinGen_dataset, "Applied Evidence Codes (Met)")
> BCF(ClinGen_dataset, "Applied Evidence Codes (Met)", 0.1, 350)

## Step 5: Visualization
# Display variant characteristics
> data("ClinGen_dataset")
> ClinGen_dataset <- add_info(ClinGen_dataset, "Assertion")
> ClinGen_dataset <- VUS_classify(ClinGen_dataset, "Assertion", "Applied Evidence Codes (Met)")
> multi_plot(ClinGen_dataset, "Assertion", "HGNC Gene Symbol")

# For LR
> data("LR_result")
> op_list <- c(2.08, 4.33, 18.70, 351)
> heatmap_LR(LR_result, "Pathogenic", op_list)

# For lr
> data("lr_CI_result")
> postp_list <- c(0.100, 0.211, 0.608, 0.981)
> plot_lr(lr_CI_result, "Pathogenic", postp_list)

```

Furthermore, BayesQuantify could play a crucial role in addressing the issue of double-counting evidence. BayesQuantify has incorporated correlation checks to reduce this bias:
``` r
> library(BayesQuantify)
> data("ClinGen_dataset")
> evidence_corplot(ClinGen_dataset, "Applied Evidence Codes (Met)")
> evaluation_variant(ClinGen_dataset, "Applied Evidence Codes (Met)")
```
The `evaluation_variant` function can automatically identify variants where high inter-evidence correlations suggest potential double-counting risk. 
These variants are flagged for user re-evaluation. The `evidence_corplot` will generate a correlation plot like below:

![Evidence correlation](https://github.com/liusihan/BayesQuantify/blob/master/man/figures/evidence_corrplot.jpg?raw=true)


## Citation
If you use BayesQuantify, please cite our paper (thanks!):
> Liu S, Feng X, Wu Y, et alCalibration and refinement of ACMG/AMP criteria for variant classification with BayesQuantifyJournal of Medical Genetics Published Online First: 19 September 2025. doi: 10.1136/jmg-2025-110863.


## Getting help
If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub. For questions and other discussions, please contact Sihan Liu (liusihan@wchscu.cn).
