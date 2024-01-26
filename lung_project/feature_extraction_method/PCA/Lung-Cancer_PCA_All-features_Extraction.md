Lung cancer: Run Principal Component Analysis (PCA) on all featuresr
================

### Load data

``` r
lung_data <- read.csv("/Users/justine_labory/Desktop/These_ANR/Latent_space_metabolomics/lung_project/LUNG.ALL.FEAT.csv")
```

# Run PCA

### Compute PCA on LUNG data

``` r
data_pos.pca = prcomp(lung_data[, -c(1)], center = TRUE, scale. = TRUE)
names(data_pos.pca)
```

    ## [1] "sdev"     "rotation" "center"   "scale"    "x"

### Look at PCA components

``` r
pos_pca.var <- data_pos.pca$sdev ^ 2
pos_pca.pvar <- pos_pca.var/sum(pos_pca.var)
paste0("Min & Max proportions of variance: ", min(pos_pca.pvar), " & ", max(pos_pca.pvar))
```

    ## [1] "Min & Max proportions of variance: 5.71622912816905e-31 & 0.0691084936986534"

### Save loading matrix

``` r
loadings_matrix <- as.data.frame(data_pos.pca$rotation)
dim(loadings_matrix)
#head(loadings_matrix)
```

    ## [1] 2944 1005

### Select the principal components (PCs) that explain 90% of the variance in the data.

``` r
best_k <- sum(cumsum(pos_pca.pvar) <= 0.90)
best_k
```

    ## [1] 399

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

<img src="Lung-Cancer_PCA_All-features_Extraction_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

### Save the results of PCA

``` r
data_pos_pca_target = cbind(Label=lung_data[,1], data_pos.pca$x[, 1:400])
# write.csv(data_pos_pca_target, "/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG.PCA.csv", row.names = F)
```
