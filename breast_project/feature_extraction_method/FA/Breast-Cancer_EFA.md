BREAST Cancer: Exploratory Factor Analysis
================

### Load libraries

``` r
library(EFAtools)
```

### Load data

``` r
sub_feat_breast_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data//BREAST_T.Feat.Select.csv")
#dim(sub_feat_breast_data)
```

# Exploratory Factor Analysis (EFA) on the BREAST dataset

As the dataset is small, a feature selection step was applied before
performing the exploratory factor analysis, in order to meet the
condition that the sample size must be at least 4 times greater than the
number of variables. To do this, The Kolmogorov-Smirnov (KS) test was
used to keep the features with most significant difference between both
populations.

Now, we want to explore factors that may exist in the dataset.

## Step 1: Appropriateness of factor analysis

### Check there is no missing value

``` r
sum(is.na(sub_feat_breast_data))
```

    ## [1] 0

### Calculation of Kaiser–Meyer–Olkin (KMO) index

We compute the KMO index to measure the sampling adequacy, it indicates
the proportion of variance in your variables that may be caused by
underlying factors. The expected value of KMO index must be greater than
0.7, a value smaller than 0.5 is unacceptable as in this case, the FA
will not be useful.

``` r
EFAtools::KMO(sub_feat_breast_data[,-1])
```

    ## ℹ 'x' was not a correlation matrix. Correlations are found from entered raw data.

    ## 
    ## ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
    ## 
    ## ✔ The overall KMO value for your data is meritorious.
    ##   These data are probably suitable for factor analysis.
    ## 
    ##   Overall: 0.812
    ## 
    ##   For each variable:
    ##              X2.aminoadipic.acid          X2.hydroxyglutaric.acid 
    ##                            0.723                            0.871 
    ##           X2.hydroxyvaleric.acid              X3.phosphoglycerate 
    ##                            0.840                            0.636 
    ## X5..deoxy.5..methylthioadenosine                          adenine 
    ##                            0.863                            0.852 
    ##                          alanine          alpha.ketoglutaric.acid 
    ##                            0.874                            0.665 
    ##                 arachidonic.acid                       asparagine 
    ##                            0.725                            0.763 
    ##                    benzylalcohol                     beta.alanine 
    ##                            0.852                            0.796 
    ##                      citric.acid                       citrulline 
    ##                            0.620                            0.747 
    ##                       creatinine    cytidine.5.monophosphate.NIST 
    ##                            0.801                            0.845 
    ##                        dodecanol                     fumaric.acid 
    ##                            0.883                            0.807 
    ##                          glucose                  glucuronic.acid 
    ##                            0.802                            0.765 
    ##                    glutamic.acid                        glutamine 
    ##                            0.876                            0.827 
    ##                    glyceric.acid                         glycerol 
    ##                            0.761                            0.821 
    ##         glycerol.alpha.phosphate          glycerol.beta.phosphate 
    ##                            0.800                            0.784 
    ##                    glycolic.acid                    hydroxylamine 
    ##                            0.918                            0.825 
    ##                     hypoxanthine                 idonic.acid.NIST 
    ##                            0.845                            0.870 
    ##         inosine.5..monophosphate                           malate 
    ##                            0.865                            0.830 
    ##                          maltose                methanolphosphate 
    ##                            0.597                            0.659 
    ##           N.acetyl.D.mannosamine                  N.methylalanine 
    ##                            0.620                            0.635 
    ##                     nicotinamide                      oxalic.acid 
    ##                            0.828                            0.824 
    ##                   p.hydroquinone              phosphoethanolamine 
    ##                            0.849                            0.848 
    ##                  phosphoric.acid                          proline 
    ##                            0.843                            0.863 
    ##                   pseudo.uridine     pyrazine.2.5.dihydroxy..NIST 
    ##                            0.733                            0.894 
    ##                           serine                    threonic.acid 
    ##                            0.800                            0.715 
    ##             tocopherol.beta.NIST                        trehalose 
    ##                            0.662                            0.589 
    ##                           uracil                             urea 
    ##                            0.848                            0.602 
    ##         uridine.5..monophosphate                         xanthine 
    ##                            0.842                            0.773

### Calculation of Bartlett’s test of Spericity

We examine the hypothesis that variables are uncorrelated in the
population. If variables are uncorrelated, this means FA is not
appropriate (will not make sense). Expected value of r should be smaller
than 0.05: means data does not produce an identity matrix.

``` r
EFAtools::BARTLETT(sub_feat_breast_data[,-1])
```

    ## ℹ 'x' was not a correlation matrix. Correlations are found from entered raw data.

    ## 
    ## ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
    ##   These data are probably suitable for factor analysis.
    ## 
    ##   𝜒²(1326) = 8368.7, p < .001

**Conclusion** : The data are suitable for FA.

## Step 2: Decision according to the number of factors

Expected number of factors: Cumulative factors explaining 60% - 70% of
variance should be retained in the model.

``` r
N_FACTORS(sub_feat_breast_data[,-1], method = "ULS")
```

    ##                                                                                                                                                                   🏃 ◯ ◯ ◯ ◯ ◯ ◯ Running CD                                                                                                                                                                 ◉ 🏃 ◯ ◯ ◯ ◯ ◯ Running EKC                                                                                                                                                                 ◉ ◉ 🏃 ◯ ◯ ◯ ◯ Running HULL                                                                                                                                                                 ◉ ◉ ◉ 🏃 ◯ ◯ ◯ Running KGC                                                                                                                                                                 ◉ ◉ ◉ ◉ 🏃 ◯ ◯ Running PARALLEL                                                                                                                                                                 ◉ ◉ ◉ ◉ ◉ 🏃 ◯ Running SCREE                                                                                                                                                                 ◉ ◉ ◉ ◉ ◉ ◉ 🏃  Running SMT                                                                                                                                                                 ◉ ◉ ◉ ◉ ◉ ◉ ◉ Done!

    ## 
    ## ── Tests for the suitability of the data for factor analysis ───────────────────
    ## 
    ## Bartlett's test of sphericity
    ## 
    ## ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
    ##   These data are probably suitable for factor analysis.
    ## 
    ##   𝜒²(1326) = 8368.7, p < .001
    ## 
    ## Kaiser-Meyer-Olkin criterion (KMO)
    ## 
    ## ✔ The overall KMO value for your data is meritorious with 0.812.
    ##   These data are probably suitable for factor analysis.
    ## 
    ## ── Number of factors suggested by the different factor retention criteria ──────
    ## 
    ## ◌ Comparison data: 10
    ## ◌ Empirical Kaiser criterion: 7
    ## ◌ Hull method with CAF: 1
    ## ◌ Hull method with CFI: 7
    ## ◌ Hull method with RMSEA: 7
    ## ◌ Kaiser-Guttman criterion with PCA: 14
    ## ◌ Kaiser-Guttman criterion with SMC: 10
    ## ◌ Kaiser-Guttman criterion with EFA: 7
    ## ◌ Parallel analysis with PCA: 9
    ## ◌ Parallel analysis with SMC: 20
    ## ◌ Parallel analysis with EFA: 10
    ## ◌ Sequential 𝜒² model tests: 26
    ## ◌ Lower bound of RMSEA 90% confidence interval: 17
    ## ◌ Akaike Information Criterion: 23

![](Breast-Cancer_EFA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](Breast-Cancer_EFA_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](Breast-Cancer_EFA_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

The most frequently recurring numbers of factors are: 7 and 10.
