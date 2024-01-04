Breast Cancer: Run Mixture of PPCA (MPPCA) on selected features
================
Justine LABORY
2023-09-15

### Load libraries

``` r
library(HDclassif)
library(pROC)
library(data.table)
```

### Load data

``` r
sub_feat_breast_data <- read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/breast_project/data/BREAST_T.Feat.Select.csv")
```

# Split data in train and test

``` r
set.seed(1234)
# Split original dataset in values (X) and labels (Y)
X <- sub_feat_breast_data[,-c(1)]
y <- factor(sub_feat_breast_data$Label)

# Define the number of folds for CV
k_folds <- 4       

# Add column for fold identifier
data_ex <- data.table(y, X, fold = 0)                              

# Create stratified folds
for (y_i in levels(y) ) {                                          
  nrow_i       <- nrow(data_ex[y == y_i,])
  n_per_fold_i <- ceiling(nrow_i / k_folds)
  data_ex[y == y_i, fold := sample(rep(1:k_folds, n_per_fold_i), nrow_i, replace = FALSE)]
}

# Number of classes per fold
data_ex[, table(y, fold)] 
```

# 2. Run MPPCA

``` r
accuracy = 0
balanced_accuracy = 0
precision = 0
recall = 0
specificity = 0
f1_score = 0
roc_auc = 0

# 5Fold CV
for (i in 1:k_folds) {
  
  ind_train <- data_ex$fold != i
  train_data <- cbind(y, X)[ind_train, ]
  test_data <- cbind(y, X)[!ind_train, ]

  # Fit the model
  prms <- hdda(train_data[,-1], model="AkjBkQkD", train_data[,1], scaling=TRUE)
  
  #the results on the test dataset
  pred_test <- predict(prms, test_data[,-1], test_data[,1])
  
  # Confusion matrix
  confusion_mat <- pred_test$confusion
  
  TP <- confusion_mat[2, 2]
  TN <- confusion_mat[1, 1]
  FP <- confusion_mat[2, 1]
  FN <- confusion_mat[1, 2]
  
  # Compute metrics
  
  accuracy[i] <- mean(pred_test$class == test_data[, 1])
  precision[i] <- TP / (TP + FP)
  ## Recall or sensitivity (True Positive Rate)
  recall[i] <- TP / (TP + FN)
  ## Specificity (True Negative Rate)
  specificity[i] <- TN / (TN + FP)
  
  balanced_accuracy[i] <- (recall[i] + specificity[i]) / 2
  
  f1_score[i] <- 2 * (precision[i] * recall[i]) / (precision[i] + recall[i])

  posterior = pred_test$posterior
  res.roc <- roc(test_data[, 1], posterior[,2])
  # ROC AUC
  roc_auc[i] <- auc(res.roc)[1]
}
```

# 3. Compute the mean of metrics

``` r
paste0("Mean Accuracy : ", round(mean(accuracy)*100,1), " (±", round(sd(accuracy)*100,1),")")
paste0("Mean Balanced accuracy : ", round(mean(balanced_accuracy)*100,1), " (", round(sd(balanced_accuracy)*100,1),")")
paste0("Mean Precision : ", round(mean(precision)*100,1), " (±", round(sd(precision)*100,1),")")
paste0("Mean Recall : ", round(mean(recall)*100,1), " (±", round(sd(recall)*100,1),")")
paste0("Mean F1 score : ", round(mean(f1_score)*100,1), " (±", round(sd(f1_score)*100,1),")")
paste0("Mean ROC AUC : ", round(mean(roc_auc)*100,1), " (±", round(sd(roc_auc)*100,1),")")
paste0("Mean Specificity : ", round(mean(specificity)*100,1), " (±", round(sd(specificity)*100,1),")")
```

    ## [1] "Mean Accuracy : 87.5 (±3.5)"
    ## [1] "Mean Balanced accuracy : 84.2 (4.5)"
    ## [1] "Mean Precision : 92.6 (±2.9)"
    ## [1] "Mean Recall : 90.7 (±4.9)"
    ## [1] "Mean F1 score : 91.6 (±2.5)"
    ## [1] "Mean ROC AUC : 89.9 (±4.1)"
    ## [1] "Mean Specificity : 77.8 (±9.8)"

# 4. Compute the 95% Confidence Interval

``` r
paste0("Accuracy 95% CI: [", round(quantile(accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Balanced accuracy 95% CI: [", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Precision 95% CI: [", round(quantile(precision, probs = c(0.025, 0.975))[1]*100, 1), " ; ",round(quantile(precision, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Recall 95% CI: [", round(quantile(recall, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(recall, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("F1 score 95% CI: [", round(quantile(f1_score, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(f1_score, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("ROC AUC 95% CI: [", round(quantile(roc_auc, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(roc_auc, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Specificity 95% CI: [", round(quantile(specificity, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(specificity, probs = c(0.025, 0.975))[2]*100, 1),"]")
```

    ## [1] "Accuracy 95% CI: [82.8 ; 89.7]"
    ## [1] "Balanced accuracy 95% CI: [80.4 ; 88.7]"
    ## [1] "Precision 95% CI: [89.3 ; 95.7]"
    ## [1] "Recall 95% CI: [84.8 ; 95.8]"
    ## [1] "F1 score 95% CI: [88.1 ; 93.1]"
    ## [1] "ROC AUC 95% CI: [84.8 ; 93.9]"
    ## [1] "Specificity 95% CI: [65.6 ; 87.1]"
