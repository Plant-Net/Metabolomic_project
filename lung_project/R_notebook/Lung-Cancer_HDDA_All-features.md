Lung Cancer: Run High dimensional discriminant analysis (HDDA) on all
features
================

### Load libraries

``` r
library(HDclassif)
library(pROC)
```

### Load data

``` r
lung_data <- as.matrix(read.csv("/Users/justine_labory/Desktop/github/plantnet/Metabolomic_project/lung_project/data/LUNG.ALL.FEAT.csv"))
```

# 1. Run HDDA

``` r
#Set the same seed
set.seed(2023)
#
B.iter = 20
accuracy = 0
balanced_accuracy = 0
precision = 0
recall = 0
specificity = 0
f1_score = 0
roc_auc = 0

#
N = nrow(lung_data)

# Run bootstrap
for (i in 1:B.iter){
  
  # Split Train/Test
  train_ratio = 0.8
  indices= sample(1:N, train_ratio*N, replace = TRUE)
  train_data= lung_data[indices,]
  test_data= lung_data[-indices,]
  
  # Train model and compute performance metrics
  prms <- hdda(train_data[,-1], train_data[,1], scaling=TRUE)
  
  # Predict on test set
  pred_test <- predict(prms, test_data[,-1], test_data[,1])
  posterior = pred_test$posterior
  
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

  res.roc <- roc(test_data[, 1], posterior[,2])
  # ROC AUC
  roc_auc[i] <- auc(res.roc)[1]
  
}
```

# 2. Compute the mean of metrics

``` r
paste0("Mean Accuracy : ", round(mean(accuracy)*100,1), " (±", round(sd(accuracy)*100,1),")")
paste0("Mean Balanced accuracy : ", round(mean(balanced_accuracy)*100,1), " (±", round(sd(balanced_accuracy)*100,1),")")
paste0("Mean Precision : ", round(mean(precision)*100,1), " (±", round(sd(precision)*100,1),")")
paste0("Mean Recall : ", round(mean(recall)*100,1), " (±", round(sd(recall)*100,1),")")
paste0("Mean F1 score : ", round(mean(f1_score)*100,1), " (±", round(sd(f1_score)*100,1),")")
paste0("Mean ROC AUC : ", round(mean(roc_auc)*100,1), " (±", round(sd(roc_auc)*100,1),")")
paste0("Mean Specificity : ", round(mean(specificity)*100,1), " (±", round(sd(specificity)*100,1),")")
```

    ## [1] "Mean Accuracy : 63.3 (±2.7)"
    ## [1] "Mean Balanced accuracy : 63.6 (±2.3)"
    ## [1] "Mean Precision : 59.4 (±3.8)"
    ## [1] "Mean Recall : 69.7 (±13.5)"
    ## [1] "Mean F1 score : 63.3 (±5)"
    ## [1] "Mean ROC AUC : 67.2 (±2.2)"
    ## [1] "Mean Specificity : 57.5 (±13.9)"

# 3. Compute the 95% Confidence Interval

``` r
paste0("Accuracy 95% CI: [", round(quantile(accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Balanced accuracy 95% CI: [", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Precision 95% CI: [", round(quantile(precision, probs = c(0.025, 0.975))[1]*100, 1), " ; ",round(quantile(precision, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Recall 95% CI: [", round(quantile(recall, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(recall, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("F1 score 95% CI: [", round(quantile(f1_score, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(f1_score, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("ROC AUC 95% CI: [", round(quantile(roc_auc, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(roc_auc, probs = c(0.025, 0.975))[2]*100, 1),"]")
paste0("Specificity 95% CI: [", round(quantile(specificity, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(specificity, probs = c(0.025, 0.975))[2]*100, 1),"]")
```

    ## [1] "Accuracy 95% CI: [58.3 ; 68]"
    ## [1] "Balanced accuracy 95% CI: [60.3 ; 68]"
    ## [1] "Precision 95% CI: [52.4 ; 65]"
    ## [1] "Recall 95% CI: [47.9 ; 90.8]"
    ## [1] "F1 score 95% CI: [53.6 ; 71.5]"
    ## [1] "ROC AUC 95% CI: [63.2 ; 71.6]"
    ## [1] "Specificity 95% CI: [31.4 ; 75.1]"
