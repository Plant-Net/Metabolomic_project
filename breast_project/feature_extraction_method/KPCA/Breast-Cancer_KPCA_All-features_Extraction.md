Breast cancer: Run Kernel principal components analysis (KPCA) on all
features
================

### Load Libraries

``` r
library(kernlab)
```

### Load data

``` r
breast_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_T.ALL.Feat.csv")
```

# Run KPCA

### Compute KPCA on all features

``` r
data_pos.kpca = kpca(x=as.matrix(breast_data[, -c(1)]),  kernel = "rbfdot", features=0, kpar=list(sigma=0.2), th=0.0001)
```

### Select Kernel Principal Components (KPCs) that explain 95% of the variance in the data

``` r
var_explained <- pcv(data_pos.kpca)
cum_var_explained <- cumsum(var_explained)
best_k <- min(which(cum_var_explained >= 0.95))  # Choose the number of components to preserve at least 95% of variance
```

### Save the results of KPCA

``` r
target = "Label"
data_pos_kpca_target = cbind(Label=breast_data[,1], data_pos.kpca@rotated[,1:best_k])
#write.csv(data_pos_kpca_target,"/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_KPCA.ON.ALL.Feat.csv", row.names = F)
```
