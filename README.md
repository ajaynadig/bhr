[![DOI](https://zenodo.org/badge/511173381.svg)](https://zenodo.org/badge/latestdoi/511173381)
# BHR
*Burden Heritability Regression* (BHR) is a method to estimate the heritability explained by mutational burden in each gene. BHR is described in detail in [here](https://www.medrxiv.org/content/10.1101/2022.07.06.22277335v1). 

`BHR` is a utility implemented in R for estimating burden heritability and derived quantities, including gene set enrichments and genetic correlations.

**Setting up `BHR` on your machine**

This repository contains the `BHR` code. To download into an R session, run:

`devtools::install_github("ajaynadig/bhr")`

We have also published detailed scripts in the `example` folder, including scripts for generating gene-level summary statistics from variant-level summary statistics, scripts for running `BHR`, scripts for generating the figures in the manuscript. Scripts for running the simulations in the manuscript are in the `MATLAB` folder.

## Elements of a `BHR` analysis

1) *Variant-level summary statistics*

Overview: A text file, with a row per variant and column per required variable (note mandatory column names). For example: if you are analyzing 15,000 genes, each with 10 variants, this file will have 150,000 lines, excluding the header.

*Required columns*

a) **Gene name** (required column name: `gene`): Any gene naming convention is valid (i.e. ENSEMBL ID), as long as the convention is consistent with that used in Baseline-BHR (see below)

b) **Chromosome** (required column name: `chromosome`): The chromosome of the gene. 

c) **Gene position in base pairs** (required column name: `gene_position`): The position of the gene in base pairs. Note that these values are only used to order genes to divide them into jackknife blocks, so minor variations due to genome build, TSS vs midpoint, etc., should not meaningfully change results.

d) **Phenotype sample size** (required column name: `N`): Phenotype sample size in the association study

e) **Variant per-allele effect sizes** (required column name: `beta`): The per-allele effect size of the variant

f) **Allele frequency** (required column name: `AF`): The frequency of the allele. Note: users may also provide the variance of the allele instead of the allele frequency, with a column named `variant_variance` and setting custom_variant_variances = TRUE. 

g) **Phenotype name** (required column name: `phenotype_key`): Phenotype name, any string

*Example:*

```
             gene chromosome gene_position      N      beta         AF phenotype_key
  ENSG00000000419         20      49563248 375630  0.014129 5.1037e-06          50NA
  ENSG00000000419         20      49563248 375630  0.445650 1.2695e-06          50NA
  ENSG00000000419         20      49563248 375630 -0.640410 1.2691e-06          50NA
  ENSG00000000419         20      49563248 375630  0.250800 2.5382e-06          50NA
  ENSG00000000419         20      49563248 375630 -0.382830 1.2746e-06          50NA
  ENSG00000000419         20      49563248 375630  0.288030 2.5496e-06          50NA

```

2) *Baseline-BHR*

Overview: A text file, with a row per gene and a column per gene set annotation, with elements equal to 1 to denote gene set membership and 0 otherwise. A *Baseline-BHR* file is required for `BHR`, as failure to control for LD-dependent architecture can lead to bias in heritability estimates (analogous to motivation for baseline model in LD Score Regression). We provide *Baseline-BHR* files with annotations corresponding to 1/5th of the observed/expected loss-of-function distribution (see manuscript and reference_files in this repository). `BHR` will also estimate genetic architecture parameters for annotations in the baseline model.

*Required variables*

a) **Gene name** (required column name: `gene`): Same gene name convention as in the *Gene-level summary statistics* file

b) **Gene membership annotations** (required column names: no restrictions): 1 or 0 to denote presence/absence of gene in gene set. Note: if intercept = TRUE, the union of baseline annotations must not span all genes to avoid colinearity. As seen in example below, we avoid colinearity by omitting a single baseline annotation from the regression.

*Example:*
```
             gene baseline_oe1 baseline_oe2 baseline_oe3 baseline_oe4
  ENSG00000000419            0            0            1            0
  ENSG00000000457            0            1            0            0
  ENSG00000000460            0            0            1            0
  ENSG00000000938            0            1            0            0
  ENSG00000000971            0            1            0            0
  ENSG00000001036            0            0            1            0
```

3) *Gene set annotations*

Overview: `BHR` can accept an arbitrary number of gene sets, in addition to the *Baseline-BHR* annotations. `BHR` will estimate genetic architecture parameters for these gene sets. 

*Required variables*

a) **Gene name** (required column name: `gene`): Same gene name convention as in the *Gene-level summary statistics* file

b) **Gene membership annotations** (required column names: no restrictions): 1 or 0 to denote presence/absence of gene in gene set

*Example:*

```
             gene gene_set_1
  ENSG00000187634          0
  ENSG00000188976          0
  ENSG00000187961          0
  ENSG00000187583          1
  ENSG00000187642          0
  ENSG00000188290          0
```

## Overview of `BHR` usage

`BHR` can be used in different modes, depending on desired outputs. The three primary modes are:

1) Univariate: estimate heritability and genetic architecture for a single phenotype
2) Bivariate: estimate cross trait genetic correlation and genetic architecture

**Important note about frequency-function filtering**

Variants profiled through exome sequencing span functional consequence (e.g., predicted loss-of-function, missense, synonymous) and orders of magnitude of allele frequency. When running `BHR`, we recommend partitioning variants into frequency-function bins, where an individual frequency-function bin is defined by a frequency (e.g., AF < 1e-5) and a function (predicted loss-of-function). We recommend this approach for multiple reasons:

1) Improved interpretability
2) Attenuation bias when analyzing a wide allele frequency range together, due to effect-size - frequency dependent architecture
3) Attenuation bias from jointly analyzing variants with different effect sizes (e.g., predicted loss-of-function with synonymous).

In the manuscript, we define frequency bins by order of magnutide (e.g., 1e-5 < AF < 1e-4).

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
