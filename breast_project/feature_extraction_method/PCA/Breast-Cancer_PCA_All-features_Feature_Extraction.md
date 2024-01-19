Brain cancer: Run Principal Component Analysis (PCA) on all features
================

### Load data

``` r
breast_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_T.ALL.Feat.csv")
```

# Run PCA

### Compute PCA on BREAST data

``` r
data_pos.pca = prcomp(breast_data[, -c(1)], center = TRUE, scale. = TRUE)
names(data_pos.pca)
```

    ## [1] "sdev"     "rotation" "center"   "scale"    "x"

### Look at PCA components

``` r
pos_pca.var <- data_pos.pca$sdev ^ 2
pos_pca.pvar <- pos_pca.var/sum(pos_pca.var)
paste0("Min & Max proportions of variance: ", min(pos_pca.pvar), " & ", max(pos_pca.pvar))
```

    ## [1] "Min & Max proportions of variance: 2.42765835211509e-05 & 0.172079220531903"

### Save loading matrix

``` r
loadings_matrix <- as.data.frame(data_pos.pca$rotation)
dim(loadings_matrix)
#head(loadings_matrix)
```

    ## [1] 162 162

### Select the principal components (PCs) that explain 90% of the variance in the data.

``` r
best_k <- sum(cumsum(pos_pca.pvar) <= 0.90)
best_k
```

    ## [1] 62

Plot PCA

``` r
par(mfrow=c(2,2))
plot(pos_pca.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
plot(cumsum(pos_pca.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
abline(h=0.9, col="blue")
abline(v=best_k, col="red")
screeplot(data_pos.pca)
screeplot(data_pos.pca,type="l")
par(mfrow=c(1,1))
```

<img src="Breast-Cancer_PCA_All-features_Feature_Extraction_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

### Save the results of PCA

``` r
data_pos_pca_target = cbind(Label=breast_data[,1], data_pos.pca$x[, 1:best_k])
#write.csv(data_pos_pca_target,"/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_PCA.ON.ALL.Feat.csv", row.names = F)
```
