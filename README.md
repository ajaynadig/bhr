# BHR
*Burden Heritability Regression* (BHR) is a method to estimate the heritability explained by mutational burden in each gene. BHR is described in detail in [here](https://www.medrxiv.org/content/10.1101/2022.07.06.22277335v1). 

`BHR` is a utility implemented in R for estimating burden heritability and derived quantities, including gene set enrichments and genetic correlations.

**Setting up `BHR` on your machine**

This repository contains the `BHR` code. To download into an R session, run:

`devtools::install_github("ajaynadig/bhr")`

We have also published detailed scripts in the `example` folder, including scripts for generating gene-level summary statistics from variant-level summary statistics, scripts for running `BHR`, scripts for generating the figures in the manuscript. Scripts for running the simulations in the manuscript are in the `MATLAB` folder.

## Elements of a `BHR` analysis

1) *Gene-level summary statistics*

Overview: A text file, with a row per gene and column per required variable (note mandatory column names). The gene-level summary statistics are frequently generated from aggregating across variant-level summary statistics.

*Required columns*

a) **Gene name** (required column name: `gene`): Any gene naming convention is valid (i.e. ENSEMBL ID), as long as the convention is consistent with that used in Baseline-BHR (see below)

b) **Chromosome** (required column name: `chromosome`): The chromosome of the gene. Note that these values are only used to order genes to divide them into jackknife blocks, so minor variations due to genome build, TSS vs midpoint, etc., should not change results substantially.

c) **Gene position in base pairs** (required column name: `gene_position`): The position of the gene in base pairs. Note that these values are only used to order genes to divide them into jackknife blocks, so minor variations due to genome build, TSS vs midpoint, etc., should not change results substantially.

d) **Phenotype sample size** (required column name: `N`): Phenotype sample size in the association study

e) **Variance variances** (required column name: `variant_variances`): This is a list of variant variances for variants in the gene. The variance variance = 2*p*(1-p), where p it the frequency of the minor allele in the association study. Since this is a list of the variant variances for each variant in the gene, the length of the list should be equal to the number of variants in the gene; note that summary statistics are often stratified by variant frequency and function (i.e., ultra-rare loss-of-function variants), and it is the number of variants in that frequency-function class which are included in the summary statistics.

f) **Variant per-allele effect sizes** (required column name: `betas`): This is a list of per-allele effect sizes for variants in the gene. The length of the list for each gene is the same as the length of the list of `variant_variances` (see above)

2) *Baseline-BHR*

Overview: A text file, with a row per gene and a column per gene set annotation, with elements equal to 1 to denote gene set membership and 0 otherwise. A *Baseline-BHR* file is required for `BHR`, as failure to control for LD-dependent architecture can lead to bias in heritability estimates (analogous to motivation for baseline model in LD Score Regression). We provide *Baseline-BHR* files with annotations corresponding to 1/5th of the observed/expected loss-of-function distribution (see manuscript and reference_files in this repository). `BHR` will also estimate genetic architecture parameters for annotations in the baseline model.

*Required variables*

a) **Gene name** (required column name: `gene`): Same gene name convention as in the *Gene-level summary statistics* file

b) **Gene membership annotations** (required column names: no restrictions): 1 or 0 to denote presence/absence of gene in gene set

3) *Gene set annotations*

Overview: `BHR` can accept an arbitrary number of gene sets, in addition to the *Baseline-BHR* annotations. `BHR` will estimate genetic architecture parameters for these gene sets. 

*Required variables*

a) **Gene name** (required column name: `gene`): Same gene name convention as in the *Gene-level summary statistics* file

b) **Gene membership annotations** (required column names: no restrictions): 1 or 0 to denote presence/absence of gene in gene set

## Overview of `BHR` usage

`BHR` can be used in different modes, depending on desired outputs. The three primary modes are:

1) Univariate: estimate heritability and genetic architecture for a single phenotype
2) Bivariate: estimate cross trait genetic correlation and genetic architecture

**Univariate `BHR` analysis**

Sample command:

```
BHR(mode = "univariate", 
    trait1_sumstats = sumstats,
      annotations = list(baseline, annotation_1))
```

*Required flags:*

1) `mode`: For univariate analysis, select "univariate"
2) `trait1_sumstats`: The gene-level summary statisics file described above, filtered to the phenotype of interest
3) `annotations`: A list of gene annotations, including the baseline file (required), and any additional gene set annotations

The output of `univariate` is an R object. Of interest to most users will be:

1) `output_name$mixed_model$heritabilities`: Reports the burden h2 and burden h2 standard error for each annotation and total heritability
2) `output_name$mixed_model$enrichments`: Reports the burden h2 enrichment and burden h2 enrichment standard error for each annotation

**Bivariate `BHR` analysis (Genetic Correlation)**

Sample command:

```
BHR(mode = "bivariate", 
    trait1_sumstats = sumstats_1,
    trait2_sumstats = sumstats_2
      annotations = list(baseline, annotation_1))
```

*Required flags:*

1) `mode`: For bivariate analysis, select "bivariate"
2) `trait1_sumstats`: The gene-level summary statisics file described above, filtered to phenotype 1
3) `trait2_sumstats`: The gene-level summary statisics file described above, filtered to phenotype 2
4)  `annotations`: A list of gene annotations, including the baseline file (required), and any additional gene set annotations

The output of `bivariate` is an R object. Of interest to most users will be:

1) `output_name$rg$rg_mixed` and `output_name$rg$rg_mixed_se`: Reports the burden h2 rg and burden h2 rg standard error for the trait pair
