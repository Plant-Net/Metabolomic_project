Pre-processing & Feature Selection on BRAIN Data
================

### Load libraries

``` r
library(stringr)
library(dplyr)
library(ggplot2)
library(GGally)
library(dgof)#install.packages("dgof")
```

### Load data

We load raw data from the article of Chardin *et al*
(<https://doi.org/10.1186/s12859-022-04900-x>).

``` r
data = read.csv("Metabolomic_project/brain_project/data/BRAIN_raw.csv", sep = ";")
t_data = t(data)
colnames(t_data) <- data[, 1]
brain_data = t_data[-c(1), ]
# Remove row names
rownames(brain_data) <- NULL
dim(brain_data)
```

    ## [1]   88 7018

The data concern **brain cancer** and are composed of **88 patients**
and **7’017 metabolites**.

# Modifications of data

### 1. Convert to numeric

``` r
brain_data <- as.data.frame(apply(brain_data, 2, function(x) as.numeric(x)))
```

### 2. Have a look on data

``` r
brain_data %>%
  mutate(Label = factor(Label)) %>%
  ggplot(aes(x=Label, fill=Label)) +
  geom_bar(stat="count", col = "black", show.legend = F) +
  scale_fill_manual(values = c("dodgerblue", "tomato3")) +
  labs(y = "Number of samples", title = "Number of samples per Label") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```

![](Brain-Cancer_Feature_Selection_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

**Labels 1** correspond to samples with **WT tumors** (n=38).<br>
**Labels 2** correspond to samples with **mutant tumors** (n=50).

Given that we have a very small number of patients, we can consider the
data to be unbalanced.

``` r
ggpairs(brain_data, columns = 12:18, ggplot2::aes(colour=as.character(brain_data[, 1]))) + theme_bw()
```

![](Brain-Cancer_Feature_Selection_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### 3. Transform target from 1/2 to 0/1

``` r
# Case: 2 >> 1 | Control: 1 >> 0
if (sum(brain_data[, 1] == 0) == 0) brain_data[, 1] = brain_data[, 1] - 1
# Look at first 5 rows
brain_data[1:5, 1]
#write.csv(sub_feat_brain_data, "Metabolomic_project/brain_project/data/BRAIN.ALL.FEAT.csv", row.names = F)
```

    ## [1] 1 1 1 1 0

Now, **labels 0** correspond to samples with **WT tumors** (n=38) and
**labels 1** correspond to samples with **mutant tumors** (n=50).

# Feature selection

Next, we want to reduce the number of features to see whether or not
this improves sample classification.<br> To do this, we chose the
Kolmogorov-Smirnov (KS) test to compare distribution of both populations
for each feature (metabolite).

### 1. Run KS test

``` r
Nb_features = NCOL(brain_data)
selected_features = c()
meta_pValues = rep(NA, 2)
#
target = "Label"
#
for (meta in 2:Nb_features){
  metabolite_data = brain_data[, meta]
  case_data = metabolite_data[brain_data[,target] == 1]
  control_data = metabolite_data[brain_data[,target] == 0]
  #Run KS test
  ks_result = ks.test(case_data, control_data, exact=FALSE)
  if (ks_result$p.value <= 0.05){ # keep the 200 significant p-value
    # Sufficient evidence to say that the two sample datasets do not come from the same distribution: so we keep the metabolite
    selected_features = c(selected_features, meta)
  }
  # Store the p.values
  meta_pValues[meta-1] = ks_result$p.value
}
```

### 2. Select all significative different features

We select only features with a pvalue \<= 0.05.

``` r
length(selected_features)
kept_feat <- selected_features
```

    ## [1] 269

In total, we have selected **269** of the 7’017 original features.

### 3. Have a look on selected features

``` r
set.seed(1)
meta = sample(kept_feat, 1) # Pick a feature
metabolite_data = brain_data[, meta]
case_data = metabolite_data[brain_data[,target] == 1]
control_data = metabolite_data[brain_data[,target] == 0]
plot(density(case_data), col="red") 
lines(density(control_data), col="blue")
```

![](Brain-Cancer_Feature_Selection_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
metabolite_data = brain_data[, 5]
case_data = metabolite_data[brain_data[,target] == 1]
control_data = metabolite_data[brain_data[,target] == 0]
plot(density(case_data), col="red") 
lines(density(control_data), col="blue")
```

![](Brain-Cancer_Feature_Selection_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### 4. Save the results

``` r
sub_feat_brain_data = cbind(Label=brain_data[,1], brain_data[, kept_feat])
dim(sub_feat_brain_data)
#write.csv(sub_feat_brain_data, "Metabolomic_project/brain_project/data/BRAIN_T.269.Feat.Select.csv", row.names = F)
```

    ## [1]  88 270
