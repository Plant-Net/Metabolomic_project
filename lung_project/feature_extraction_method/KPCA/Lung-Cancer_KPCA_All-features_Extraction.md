Lung cancer: Run Kernel principal components analysis (KPCA) on all
features
================

### Load Libraries

``` r
library(kernlab)
```

### Load data

``` r
lung_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG.ALL.FEAT.csv")
```

# Run KPCA

### Compute KPCA on all features

``` r
data_pos.kpca = kpca(x=as.matrix(lung_data[, -c(1)]),  kernel = "polydot", features=200, kpar=list(degree=3, scale=0.1, offset=1), th=0.0001)
```

### Save the results of KPCA

``` r
target = "Label"
data_pos_kpca_target = cbind(Label=lung_data[,1], data_pos.kpca@rotated)
#dim(data_pos_kpca_target)

data_pos_kpca_target <- as.data.frame(data_pos_kpca_target)
names(data_pos_kpca_target) <- c("Label", paste0("KPC",seq(1,200)))

#write.csv(data_pos_kpca_target, "/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG.KPCA.ALL.FEAT.csv", row.names = F)
```
