LUNG Cancer: Factor Analysis with 17 Factors
================

### Load libraries

``` r
library(EFAtools)
```

### Load input data

``` r
sub_feat_lung_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG_T.Feat.Select.csv")
```

# Compute factors

We use the varimax method which is the most common method.

``` r
res_efa_17 = EFA(sub_feat_lung_data[,-1], n_factors = 17, rotation = "varimax", type = "EFAtools")
```

    ## ℹ 'x' was not a correlation matrix. Correlations are found from entered raw data.

Get factor scores.

``` r
fac_scores <- FACTOR_SCORES(sub_feat_lung_data[,-1], f = res_efa_17)
data_pos_scores <- fac_scores$scores
#dim(data_pos_scores)
```

Add target variable and save the data.

``` r
target = "Label"
data_pos_fa_target = cbind(Label=sub_feat_lung_data[,1], data_pos_scores)
dim(data_pos_fa_target)

# write.csv(data_pos_fa_target, "/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG.FA.17F.csv", row.names = F)
```

    ## [1] 1005   18
