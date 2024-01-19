Breast cancer: Model Binary Classification - Xgboost classifier - 7
Factors
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
res_efa = EFA(sub_feat_breast_data[,-1], n_factors = 7, rotation = "varimax", type = "EFAtools")
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

#write.csv(data_pos_fa_target, "/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/", row.names = F)
```

    ##      Label         F2         F7          F3         F1         F6         F5
    ## [1,]     0  0.6516441 -0.9734636  1.45333134  1.2390689  0.0733995 -0.7732051
    ## [2,]     0  1.5623772 -0.1993146 -0.59190823 -0.7697203  0.3555366 -0.7736563
    ## [3,]     0  0.5207677  1.1235617  1.02598556  0.2479018  0.9043604  0.4351496
    ## [4,]     0  1.6660238 -0.2915217 -0.95669899 -0.9751198  2.4033121 -2.0086394
    ## [5,]     0  0.5048496 -0.1214078  1.39902891  0.7419197 -0.4679644  0.5737471
    ## [6,]     0 -0.8538902 -0.1278801  0.09675444 -2.2558559 -1.6679258 -0.5349055
    ##              F4
    ## [1,]  0.5324827
    ## [2,] -1.7554015
    ## [3,] -0.2029116
    ## [4,] -1.0451945
    ## [5,]  0.6413497
    ## [6,] -0.4972138
