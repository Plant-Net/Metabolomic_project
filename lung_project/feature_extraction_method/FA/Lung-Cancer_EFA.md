LUNG Cancer: Exploratory Factor Analysis
================

### Load libraries

``` r
library(EFAtools)
```

### Load data

``` r
sub_feat_lung_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG_T.Feat.Select.csv")
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
sum(is.na(sub_feat_lung_data))
```

    ## [1] 0

### Calculation of Kaiser–Meyer–Olkin (KMO) index

We compute the KMO index to measure the sampling adequacy, it indicates
the proportion of variance in your variables that may be caused by
underlying factors. The expected value of KMO index must be greater than
0.7, a value smaller than 0.5 is unacceptable as in this case, the FA
will not be useful.

``` r
EFAtools::KMO(sub_feat_lung_data[,-1])
```

    ## ℹ 'x' was not a correlation matrix. Correlations are found from entered raw data.

    ## 
    ## ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
    ## 
    ## ✔ The overall KMO value for your data is marvellous.
    ##   These data are probably suitable for factor analysis.
    ## 
    ##   Overall: 0.903
    ## 
    ##   For each variable:
    ##  mz.59.01287386...RT.71.42227731   mz.73.02849529...RT.0.00599967 
    ##                            0.842                            0.937 
    ##  mz.93.03373494...RT.237.7226903   mz.99.00795993...RT.293.162389 
    ##                            0.431                            0.909 
    ##  mz.100.0031536...RT.247.7813984   mz.117.0552304...RT.90.5310804 
    ##                            0.887                            0.940 
    ##  mz.118.0586918...RT.72.32408512   mz.120.0453103...RT.209.387502 
    ##                            0.854                            0.856 
    ##  mz.122.0248884...RT.22.88221777  mz.125.0967611...RT.266.6393549 
    ##                            0.840                            0.838 
    ##  mz.130.0659344...RT.255.2208808  mz.131.0697087...RT.55.08024311 
    ##                            0.815                            0.936 
    ##  mz.131.0710157...RT.151.8947143  mz.141.0917095...RT.286.8086252 
    ##                            0.891                            0.905 
    ##  mz.142.0657154...RT.247.8825273  mz.152.0428684...RT.225.1072904 
    ##                            0.896                            0.606 
    ##  mz.153.0188573...RT.159.2787529  mz.174.0556436...RT.255.4197604 
    ##                            0.920                            0.826 
    ##  mz.178.0506869...RT.46.75697413  mz.187.0972355...RT.66.65890112 
    ##                            0.802                            0.844 
    ##  mz.190.0510309...RT.73.88959713  mz.191.0190272...RT.27.95505778 
    ##                            0.895                            0.919 
    ##  mz.192.0227173...RT.23.53824173  mz.206.0479066...RT.97.50586957 
    ##                            0.932                            0.854 
    ##  mz.206.0458584...RT.207.6729159  mz.208.0646181...RT.94.75749426 
    ##                            0.918                            0.872 
    ##  mz.209.0792391...RT.266.6850542  mz.210.0804896...RT.266.4169753 
    ##                            0.917                            0.923 
    ##  mz.215.1283928...RT.19.04282835   mz.216.0875212...RT.64.7802733 
    ##                            0.922                            0.863 
    ##   mz.219.027369...RT.226.1705257  mz.222.1132032...RT.326.6518114 
    ##                            0.746                            0.906 
    ##  mz.223.0616652...RT.50.61336354  mz.225.1126774...RT.319.7106863 
    ##                            0.859                            0.907 
    ##  mz.227.1284347...RT.360.5070576  mz.229.0537497...RT.319.1535666 
    ##                            0.555                            0.708 
    ##  mz.229.1439476...RT.335.2732307  mz.230.0566323...RT.319.4638675 
    ##                            0.902                            0.706 
    ##   mz.231.033041...RT.303.4775614   mz.232.0370222...RT.302.956958 
    ##                            0.655                            0.658 
    ##  mz.234.1132916...RT.300.2999706  mz.236.0922707...RT.230.1987384 
    ##                            0.791                            0.801 
    ##  mz.237.1240578...RT.257.9099072  mz.239.0916448...RT.302.1258351 
    ##                            0.743                            0.918 
    ##   mz.239.091892...RT.336.8878266    mz.239.99519...RT.103.8639763 
    ##                            0.908                            0.845 
    ##  mz.241.1076748...RT.212.6572024   mz.243.0773345...RT.247.723442 
    ##                            0.870                            0.897 
    ##  mz.244.0804605...RT.247.7182844  mz.245.0487473...RT.240.1849505 
    ##                            0.922                            0.717 
    ##  mz.245.0927751...RT.245.0846723  mz.246.0521938...RT.240.0325276 
    ##                            0.922                            0.725 
    ##  mz.246.0955269...RT.245.0999688  mz.247.0480515...RT.82.09026545 
    ##                            0.870                            0.936 
    ##  mz.248.0958022...RT.240.6692675   mz.249.098146...RT.240.6458873 
    ##                            0.730                            0.711 
    ##  mz.253.1191851...RT.114.9111053  mz.276.0007389...RT.140.8707702 
    ##                            0.910                            0.939 
    ##  mz.279.1275229...RT.328.0227537  mz.293.1161117...RT.149.6636498 
    ##                            0.445                            0.875 
    ##  mz.303.0178574...RT.150.6830591  mz.304.0207077...RT.150.6612564 
    ##                            0.712                            0.730 
    ##  mz.308.0984878...RT.18.92752141   mz.310.102618...RT.65.41692807 
    ##                            0.860                            0.907 
    ##  mz.310.9692754...RT.6.663987882   mz.311.070589...RT.145.0482519 
    ##                            0.967                            0.933 
    ##  mz.311.9722409...RT.7.177236383  mz.320.1078038...RT.140.2623736 
    ##                            0.972                            0.873 
    ##  mz.326.9424254...RT.25.48383577  mz.352.0856921...RT.313.2502606 
    ##                            0.946                            0.858 
    ##  mz.353.0897307...RT.313.0080811  mz.368.9722179...RT.6.413331425 
    ##                            0.866                            0.774 
    ##  mz.378.1016025...RT.336.1376535  mz.380.1137514...RT.333.5800212 
    ##                            0.875                            0.822 
    ##  mz.391.0273783...RT.121.3899302  mz.405.0277483...RT.27.75357078 
    ##                            0.929                            0.915 
    ##   mz.406.0326611...RT.28.0647121  mz.413.0429781...RT.144.2454436 
    ##                            0.918                            0.901 
    ##  mz.421.0002752...RT.27.82324105    mz.430.91982...RT.24.97892162 
    ##                            0.948                            0.958 
    ##  mz.435.9577544...RT.28.18022409  mz.441.1613664...RT.307.9507622 
    ##                            0.957                            0.760 
    ##  mz.575.1195077...RT.146.0944394  mz.613.3595637...RT.337.9140933 
    ##                            0.918                            0.882 
    ##  mz.627.3761655...RT.345.5595674  MZ.139.0080106....RT.266.163112 
    ##                            0.856                            0.919 
    ## MZ.152.0623154....RT.92.08466258 MZ.155.0134266....RT.247.4243894 
    ##                            0.934                            0.945 
    ## MZ.160.1339746....RT.11.59533148 MZ.176.1056353....RT.24.45721495 
    ##                            0.872                            0.905 
    ##  MZ.193.098677....RT.26.13727811 MZ.201.0551545....RT.307.8090355 
    ##                            0.718                            0.878 
    ## MZ.204.1250884....RT.32.18089594 MZ.204.1345526....RT.48.36620651 
    ##                            0.798                            0.872 
    ##   MZ.211.06049....RT.261.5162407 MZ.211.0714984....RT.89.67769158 
    ##                            0.882                            0.695 
    ## MZ.211.0948965....RT.264.4617553  MZ.212.033164....RT.141.8011052 
    ##                            0.889                            0.920 
    ## MZ.215.0172169....RT.32.24189695 MZ.215.1699905....RT.32.13363533 
    ##                            0.975                            0.972 
    ##  MZ.215.385388....RT.33.51105013 MZ.216.0207763....RT.32.54068371 
    ##                            0.967                            0.980 
    ## MZ.219.1746802....RT.13.77330432  MZ.227.0593737....RT.265.090873 
    ##                            0.812                            0.921 
    ## MZ.230.9897971....RT.32.21604228  MZ.230.9942814....RT.89.5136752 
    ##                            0.970                            0.843 
    ##  MZ.236.9996167....RT.33.1702461 MZ.237.0180091....RT.57.44198376 
    ##                            0.968                            0.933 
    ##  MZ.238.003293....RT.33.23679884 MZ.242.0470363....RT.65.85015765 
    ##                            0.971                            0.784 
    ## MZ.242.2633643....RT.106.4459779 MZ.243.1004849....RT.3.745031639 
    ##                            0.784                            0.931 
    ## MZ.244.1552798....RT.56.83217341 MZ.247.0970455....RT.23.66165344 
    ##                            0.806                            0.812 
    ## MZ.247.1325996....RT.241.3775945 MZ.247.1384435....RT.26.80803383 
    ##                            0.884                            0.804 
    ## MZ.251.1374383....RT.66.39697838 MZ.252.9745196....RT.32.25163064 
    ##                            0.947                            0.952 
    ## MZ.252.9764797....RT.91.85903111 MZ.253.1093357....RT.339.0808336 
    ##                            0.869                            0.938 
    ##  MZ.253.1446207....RT.291.876058  MZ.261.1217491....RT.259.733852 
    ##                            0.737                            0.724 
    ## MZ.264.1215224....RT.23.30921876  MZ.267.0752208....RT.47.5987377 
    ##                            0.902                            0.914 
    ##  MZ.268.0786356....RT.247.286278   MZ.268.971221....RT.39.4986666 
    ##                            0.928                            0.838 
    ## MZ.269.0906438....RT.246.6105621 MZ.269.1280232....RT.25.77123022 
    ##                            0.838                            0.912 
    ## MZ.273.9766249....RT.3.580064633 MZ.274.0931388....RT.132.7231618 
    ##                            0.887                            0.882 
    ##  MZ.275.1110199....RT.76.8900988  MZ.277.1417306....RT.324.896284 
    ##                            0.839                            0.914 
    ##  MZ.281.118992....RT.286.8747521  MZ.284.1869613....RT.241.953849 
    ##                            0.730                            0.926 
    ## MZ.286.0941817....RT.24.87727498 MZ.286.1439715....RT.6.463579088 
    ##                            0.705                            0.694 
    ##   MZ.286.201986....RT.284.740301 MZ.287.2051064....RT.284.9011055 
    ##                            0.839                            0.839 
    ## MZ.287.9919291....RT.117.5061824 MZ.289.0575398....RT.47.34458071 
    ##                            0.810                            0.901 
    ## MZ.292.1031298....RT.33.20948049 MZ.292.1586582....RT.193.7893427 
    ##                            0.803                            0.825 
    ## MZ.292.9720498....RT.112.3408676 MZ.294.1695373....RT.22.51003296 
    ##                            0.868                            0.796 
    ## MZ.298.1858108....RT.335.4980608 MZ.300.8461238....RT.28.21429718 
    ##                            0.879                            0.918 
    ## MZ.310.2017751....RT.301.5268685  MZ.311.996924....RT.141.4667774 
    ##                            0.895                            0.918 
    ## MZ.313.0857312....RT.146.5674858 MZ.314.2333368....RT.24.94394447 
    ##                            0.949                            0.911 
    ## MZ.319.1658447....RT.206.4334287 MZ.320.0974838....RT.89.10383275 
    ##                            0.905                            0.946 
    ## MZ.326.1979809....RT.265.0908809 MZ.332.0963401....RT.23.39433456 
    ##                            0.866                            0.931 
    ##  MZ.335.0677266....RT.146.150292  MZ.340.1192551....RT.319.580832 
    ##                            0.935                            0.895 
    ##  MZ.354.103273....RT.310.0341178 MZ.362.0851112....RT.139.0204419 
    ##                            0.854                            0.887 
    ## MZ.365.0463959....RT.63.61843229 MZ.365.1749732....RT.57.50138287 
    ##                            0.873                            0.787 
    ## MZ.367.1509493....RT.88.75312975 MZ.369.1565491....RT.89.08054371 
    ##                            0.919                            0.907 
    ## MZ.370.0525988....RT.21.63438275 MZ.370.2620353....RT.306.6197854 
    ##                            0.819                            0.861 
    ## MZ.372.9232556....RT.27.98128537 MZ.375.1153219....RT.2.261121688 
    ##                            0.968                            0.847 
    ## MZ.376.0842185....RT.309.9276143 MZ.380.1177542....RT.332.9323936 
    ##                            0.851                            0.861 
    ## MZ.381.0993332....RT.303.4481684  MZ.389.121976....RT.313.9559327 
    ##                            0.791                            0.433 
    ## MZ.396.1502435....RT.96.88806362 MZ.398.0670488....RT.309.3474288 
    ##                            0.847                            0.907 
    ##  MZ.407.0344674....RT.3.72067325 MZ.409.0248978....RT.34.83171973 
    ##                            0.961                            0.934 
    ## MZ.411.1649188....RT.0.336878992 MZ.412.1682325....RT.303.8321051 
    ##                            0.765                            0.657 
    ## MZ.415.9645978....RT.34.51465053 MZ.417.3378443....RT.335.7191899 
    ##                            0.966                            0.846 
    ## MZ.422.2140079....RT.131.2977363 MZ.423.0084949....RT.33.91615407 
    ##                            0.790                            0.950 
    ## MZ.424.0127573....RT.34.37186156 MZ.429.0226755....RT.33.89550352 
    ##                            0.957                            0.954 
    ## MZ.435.1999616....RT.313.5644411 MZ.437.9729528....RT.34.02426469 
    ##                            0.813                            0.959 
    ## MZ.438.9799872....RT.34.38718636  MZ.444.991422....RT.33.78759494 
    ##                            0.956                            0.965 
    ##   MZ.447.10803....RT.68.33300748 MZ.451.0076189....RT.33.65433463 
    ##                            0.947                            0.947 
    ##  MZ.453.1174655....RT.6.51171488 MZ.454.1827643....RT.84.47234231 
    ##                            0.873                            0.912 
    ## MZ.486.2571336....RT.125.0351305 MZ.497.1443441....RT.286.9818977 
    ##                            0.889                            0.775 
    ## MZ.520.9790465....RT.34.15643017 MZ.542.9586941....RT.34.39494134 
    ##                            0.943                            0.941 
    ## MZ.561.3432022....RT.328.9471712 MZ.584.2670695....RT.131.4255402 
    ##                            0.841                            0.783 
    ## MZ.595.3515218....RT.301.1663045 MZ.597.3652419....RT.330.2975575 
    ##                            0.704                            0.434 
    ## MZ.615.0353192....RT.34.28684079 MZ.631.3465895....RT.343.0206379 
    ##                            0.941                            0.905 
    ## MZ.638.3608799....RT.328.1278024  MZ.656.2017529....RT.3.99017038 
    ##                            0.863                            0.877

### Calculation of Bartlett’s test of Spericity

We examine the hypothesis that variables are uncorrelated in the
population. If variables are uncorrelated, this means FA is not
appropriate (will not make sense). Expected value of r should be smaller
than 0.05: means data does not produce an identity matrix.

``` r
EFAtools::BARTLETT(sub_feat_lung_data[,-1])
```

    ## ℹ 'x' was not a correlation matrix. Correlations are found from entered raw data.

    ## 
    ## ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
    ##   These data are probably suitable for factor analysis.
    ## 
    ##   𝜒²(20301) = 218007.6, p < .001

**Conclusion** : The data are suitable for FA.

## Step 2: Decision according to the number of factors

Expected number of factors: Cumulative factors explaining 60% - 70% of
variance should be retained in the model.

``` r
N_FACTORS(sub_feat_lung_data[,-1], method = "ULS")
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
    ##   𝜒²(20301) = 218007.6, p < .001
    ## 
    ## Kaiser-Meyer-Olkin criterion (KMO)
    ## 
    ## ✔ The overall KMO value for your data is marvellous with 0.903.
    ##   These data are probably suitable for factor analysis.
    ## 
    ## ── Number of factors suggested by the different factor retention criteria ──────
    ## 
    ## ◌ Comparison data: 21
    ## ◌ Empirical Kaiser criterion: 21
    ## ◌ Hull method with CAF: 1
    ## ◌ Hull method with CFI: 17
    ## ◌ Hull method with RMSEA: 17
    ## ◌ Kaiser-Guttman criterion with PCA: 44
    ## ◌ Kaiser-Guttman criterion with SMC: 31
    ## ◌ Kaiser-Guttman criterion with EFA: 20
    ## ◌ Parallel analysis with PCA: 26
    ## ◌ Parallel analysis with SMC: 44
    ## ◌ Parallel analysis with EFA: 27
    ## ◌ Sequential 𝜒² model tests: NA
    ## ◌ Lower bound of RMSEA 90% confidence interval: 24
    ## ◌ Akaike Information Criterion: 85

![](Lung-Cancer_EFA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](Lung-Cancer_EFA_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](Lung-Cancer_EFA_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

The most frequently recurring numbers of factors are: 17 and 44.
