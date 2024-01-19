BREAST Cancer: Factor Analysis with 10 Factors
================

### Load libraries

``` r
library(EFAtools)
```

### Load input data

``` r
sub_feat_breast_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_T.Feat.Select.csv")
```

# Compute factors

We use the varimax method which is the most common method.

``` r
res_efa = EFA(sub_feat_breast_data[,-1], n_factors = 10, rotation = "varimax", type = "EFAtools")
```

    ## ℹ 'x' was not a correlation matrix. Correlations are found from entered raw data.

Get factor scores.

``` r
fac_scores <- FACTOR_SCORES(sub_feat_breast_data[,-1], f = res_efa)
data_pos_scores <- fac_scores$scores
#dim(data_pos_scores)
```

Add target variable and save the data.

``` r
target = "Label"
data_pos_fa_target = cbind(Label=sub_feat_breast_data[,1], data_pos_scores)
head(data_pos_fa_target)

#write.csv(data_pos_fa_target, "/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_FA_10Factors.csv", row.names = F)
```

    ##      Label         F2         F7         F8          F1         F5          F6
    ## [1,]     0  0.7254082  0.8603587 -0.8057560  1.23624624 -0.9744746 -0.42942617
    ## [2,]     0  1.3137312 -0.4451486 -0.7740610 -0.62764495 -0.3360028  0.03290066
    ## [3,]     0  0.5154008  1.6531212  0.4817968  0.02562832  0.4695392  1.01780082
    ## [4,]     0  1.5495201 -1.2373516 -0.3348148 -0.70259325 -1.6709145  1.96141846
    ## [5,]     0  0.4946478  1.4278745 -0.5783482  0.63650282  0.5719200 -0.81583364
    ## [6,]     0 -0.8876144  0.1002976 -0.6655145 -2.29436178 -0.3059402 -2.32241804
    ##               F4          F3         F10         F9
    ## [1,] -0.03611056  1.68385557  1.65751186  0.6157526
    ## [2,] -1.27064263  0.63267813 -1.30000744  0.7485734
    ## [3,] -0.03605994 -0.03863892 -0.02795338  0.8092779
    ## [4,] -1.25783749  1.71345409 -1.01617004  0.3224478
    ## [5,] -0.05113923  0.49058873  0.87345601 -1.0438103
    ## [6,] -1.53581656  0.65241009  0.65932164  0.6142187
