# BHR
*Burden Heritability Regression* (BHR) is a method to estimate the heritability explained by mutational burden in each gene. BHR is described in detail in here (manuscript link). 

`BHR` is a utility implemented in R for estimating burden heritability and derived quantities, including gene set enrichments and genetic correlations.

**Setting up `BHR` on your machine**

This repository contains the `BHR` code. To download:

`git clone https://github.com/danjweiner/rv_h2.git`

## Elements of a `BHR` analysis

1) *Gene-level summary statistics*

Overview: A text file, with a row per gene and column per required variable (note mandatory column names). The gene-level summary statistics are frequently generated from aggregating across variant-level summary statistics -- see below in (section) for overview and scripts. The provided variant-to-gene script automatically generates the following columns:

*Required variables*

a) **Gene name** (required column name: `gene`): Any gene naming convention is valid (i.e. ENSEMBL ID), as long as the convention is consistent with that used in Baseline-BHR (see below)

b) **Gene burden score** (required column name: `burden_score`): Sum of the variances of variants in a gene (often restricted to variants of a particular functional class and/or frequency). This is the indepedent variable in the heritability regression.

c) **Sum of variance-weighted per-allele effects** (required column name: `w_t_beta`): This is the sum of the effects of variants in a gene, where each per allele effect has been multiplied by the frequency variance of the variant (2p(1-p)). This becomes the dependent variable in the heritability regression after a modest under-the-hood `BHR` transformation (squaring and normalizing by cumulative variant frequencies).

d) **Overdispersion** (required column name: `overdispersion`): Weighted mean heterozygosity (weighted by the heterozygosity)

2) *Baseline-BHR*

Overview: A text file, with a row per gene and a column per gene set annotation, with elements equal to 1 to denote gene set membership and 0 otherwise. A *Baseline-BHR* file is required for `BHR`, as failure to control for LD-dependent architecture can lead to bias in heritability estimates (analogous to motivation for baseline model in LD Score Regression). We provide *Baseline-BHR* files with annotations corresponding to 1/5th of the observed/expected loss-of-function distribution (see manuscript). `BHR` will also estimate genetic architecture parameters for annotations in the baseline model.

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
3) Aggregate: estimate average heritability across traits, summed across variant groups (for example, across loss-of-function and missense variants)

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

*Optional flags*

1) `num_blocks`: Number of jackknife blocks (default = 100)
