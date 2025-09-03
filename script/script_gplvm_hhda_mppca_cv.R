#!/usr/bin/Rscript

library(optparse)
library(kernlab)#install.packages("kernlab")
library(HDclassif)
library(pROC)
library(data.table)
library(crayon)
library(ROCR)

############## Function


define_stratified_kfold <- function(X,y, nb_fold = 4) {
  ### """
  # The number of fold is set at 4 by default.
  # """
  
  # Define the number of folds for CV
  k_folds <- nb_fold       
  
  # Add column for fold identifier
  data_ex <- data.table(y, X, fold = 0)                              
  
  # Create stratified folds
  for (y_i in levels(y) ) {                                          
    nrow_i       <- nrow(data_ex[y == y_i,])
    n_per_fold_i <- ceiling(nrow_i / k_folds)
    data_ex[y == y_i, fold := sample(rep(1:k_folds, n_per_fold_i), nrow_i, replace = FALSE)]
  }
  
  # Number of classes per fold
  # data_ex[, table(y, fold)] 
  
  return(data_ex)
  
}

gplvm_model <- function(train_data, test_data) {
  
  # Fit the model
  gplvm <- gausspr(x=train_data[,-1], y=factor(train_data[,1]), type="classification", scaled=TRUE, kernel="rbfdot", kpar="automatic")
  
  #the results on the test dataset
  y_pred <- predict(gplvm, test_data[,-1])
  
  y_pred_proba_all <- predict(gplvm, test_data[,-1], type="probabilities")
  
  y_pred_proba <- y_pred_proba_all[,2]
  
  return(list(model = gplvm,
              y_pred = y_pred,
              y_pred_proba = y_pred_proba))
  
}

hdda_model <- function(train_data, test_data) {
  
  prms <- hdda(train_data[,-1], train_data[,1], scaling=TRUE)
  
  #the results on the test dataset
  res <- predict(prms, test_data[,-1], test_data[,1])
  
  y_pred = res$class
  
  y_pred_proba <- res$posterior[,2]
  
  
  return(list(model = prms,
              res = res,
              y_pred = y_pred,
              y_pred_proba = y_pred_proba))
  
}

mppca_model <- function(train_data, test_data) {
  
  # Fit the model
  prms <- hdda(train_data[,-1], model="AkjBkQkD", train_data[,1], scaling=TRUE)
  
  #the results on the test dataset
  res <- predict(prms, test_data[,-1], test_data[,1])
  
  y_pred = res$class
  
  y_pred_proba <- res$posterior[,2]
  
  
  return(list(model = prms,
              res = res,
              y_pred = y_pred,
              y_pred_proba = y_pred_proba))
  
}

confusion_matrix <- function(test_data, y_pred) {
  
  if (analysis_name == "GPLVM") {
    
    # Confusion matrix
    confusion_mat = as.matrix(table(Predicted_Values = y_pred, Actual_Values = test_data[, 1]))
    
  } else {
    
    confusion_mat <- res$confusion
    
  }
  
  return(confusion_mat)
  
  
}

get_tpr <- function(test_data, y_pred_proba) {
  
  roc_obj <- prediction(y_pred_proba, test_data[, 1])
  perf <- performance(roc_obj, "tpr", "fpr")
  
  fpr <- unlist(perf@x.values[[1]])
  tpr <- unlist(perf@y.values[[1]])
  auc <- round(performance(roc_obj, "auc")@y.values[[1]], 3)
  
  aucs <- c(aucs, auc)
  tpr <- approx(fpr, tpr, xout = base_fpr)$y
  tpr[1] <- 0
  
  return(tpr)

  
}

################# MAIN #######################@

# create an option parser
option_list <- list(
  make_option(c("-i", "--input"), dest="input_file",
              help="path to input file"),
  make_option(c("-o", "--output_path"), dest="output_path",
              help="path to output directory"),
  make_option(c("-n", "--number_fold"), dest="nb_fold",
              default = 4, help="number of folds for cross-validation"),
  make_option(c("-r", "--number_repetitions"), dest="nb_repet",
              default = 5, help="number of repetitions for cross-validation"),
  make_option(c("-a", "--analysis_name"), dest="analysis_name",
              help="Name of the model to use for feature extraction. It can be GPLVM, HDDA, MPPCA"),
  make_option(c("-t", "--cancer_tissue"), dest="cancer_tissue",
              help="Name of the tissue of the analysed cancer"),
  make_option(c("-m", "--omics_type"), dest="omics_type",
              help="Name of the tissue of the analysed cancer")
)



parser <- OptionParser(option_list=option_list)

# parse the command line arguments
args <- parse_args(parser)

# use the parsed arguments
input_file <- args$input_file
output_path <- args$output_path
nb_fold <- as.integer(args$nb_fold)
nb_repet <- as.integer(args$nb_repet)
analysis_name <- toupper(args$analysis_name)
cancer_tissue <- args$cancer_tissue
omics_type <- args$omics_type

cat(bgMagenta(paste0(Sys.time(),": Load data")))
cat("\n")

cancer_data <- read.csv(input_file)

# if (omics_type == "metabolomics") {
#   cancer_data <- read.csv(input_file)
# } else {
#   cancer_data <- read.csv(input_file, row.names = 1)
# }

# Split original dataset in values (X) and labels (Y)
X <- cancer_data[,-c(1)]
y <- factor(cancer_data$Label)

all_seed <- c()
accuracy <- c()
balanced_accuracy <- c()
precision <- c()
recall <- c()
specificity <- c()
f1_score <- c()
roc_auc <- c()

cat(bgMagenta(paste0(Sys.time(),": Beginning of the analysis")))
cat("\n")
start_time = Sys.time()

set.seed(2024)

for (j in 1:nb_repet) {
  
  data_ex <- define_stratified_kfold(X,y, nb_fold = nb_fold)
  
  fprs <- list()
  tprs <- list()
  base_fpr <- seq(0, 1, length.out = 1001)
  aucs <- numeric()
  
  for (i in 1:nb_fold) {
    
    ind_train <- data_ex$fold != i
    train_data <- cbind(y, X)[ind_train, ]
    test_data <- cbind(y, X)[!ind_train, ]
    
    if (analysis_name == "GPLVM") {
      
      print("GPLVM")
      
      output <- gplvm_model(train_data, test_data)
      y_pred = output$y_pred
      y_pred_proba <- output$y_pred_proba
      
    } else if (analysis_name == "HDDA") {
      
      print("HDDA")
      
      output <- hdda_model(train_data, test_data)
      res = output$res
      y_pred = output$y_pred
      y_pred_proba <- output$y_pred_proba
      
    } else if (analysis_name == "MPPCA") {
      
      print("MPPCA")
      
      output <- mppca_model(train_data, test_data)
      res = output$res
      y_pred = output$y_pred
      y_pred_proba <- output$y_pred_proba
    }
    
    
    # Confusion matrix
    confusion_mat = confusion_matrix(test_data, y_pred)
    print(confusion_mat)
    
    TP <- confusion_mat[2, 2]
    TN <- confusion_mat[1, 1]
    FP <- confusion_mat[2, 1]
    FN <- confusion_mat[1, 2]
    
    # Compute metrics
    
    accuracy <- c(accuracy, mean(y_pred == test_data[, 1]))
    
    precision_ <- TP / (TP + FP)
    precision <- c(precision, precision_)
    ## Recall or sensitivity (True Positive Rate)
    recall_ <- TP / (TP + FN)
    recall <- c(recall, recall_)
    ## Specificity (True Negative Rate)
    specificity_ <- TN / (TN + FP)
    specificity <- c(specificity, specificity_)
    balanced_accuracy <- c(balanced_accuracy, (recall_ + specificity_) / 2 )
    f1_score <- c(f1_score, 2 * (precision_ * recall_) / (precision_ + recall_))
    res.roc <- roc(test_data[, 1], y_pred_proba)
    # ROC AUC
    roc_auc <- c(roc_auc, auc(res.roc)[1])
    
    tpr <- get_tpr(test_data, y_pred_proba)
    
    # Calculate the index for DataFrame update
    index = (j - 1) * nb_fold * nb_repet + i
    
    tprs[[index]] <- tpr
    
  }
  
}

tprs <- do.call(rbind, tprs)
mean_tpr <- colMeans(tprs, na.rm = TRUE)
std <- apply(tprs, 2, sd, na.rm = TRUE)
mean_auc <- round(mean(aucs), 3)

# Save tpr
write.table(mean_tpr, paste0(output_path, "/",cancer_tissue, "_",omics_type,"_",tolower(analysis_name), "_tpr.csv"), sep = ",",row.names = F, col.names = F)

# Save df of metrics
df_metrics <- data.frame(accuracy, balanced_accuracy, precision, recall, f1_score, roc_auc,specificity)
write.csv(df_metrics, paste0(output_path, "/",cancer_tissue, "_",omics_type,"_",tolower(analysis_name), "_metrics_table.csv"), row.names = F)

df_mean_sd = as.data.frame(matrix(nrow = 2, ncol = ncol(df_metrics)))
rownames(df_mean_sd) <- c("mean", "sd")
names(df_mean_sd) <- names(df_metrics)
df_mean_sd[1,] <- c(paste0(round(mean(accuracy)*100,1), " (±", round(sd(accuracy)*100,1),")"),
                   paste0(round(mean(balanced_accuracy)*100,1), " (±", round(sd(balanced_accuracy)*100,1),")"),
                   paste0(round(mean(precision)*100,1), " (±", round(sd(precision)*100,1),")"),
                   paste0(round(mean(recall)*100,1), " (±", round(sd(recall)*100,1),")"),
                   paste0(round(mean(f1_score)*100,1), " (±", round(sd(f1_score)*100,1),")"),
                   paste0(round(mean(roc_auc)*100,1), " (±", round(sd(roc_auc)*100,1),")"),
                   paste0(round(mean(specificity)*100,1), " (±", round(sd(specificity)*100,1),")")
)

df_mean_sd[2,] <- c(paste0("[", round(quantile(accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]"),
                   paste0("[", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]"),
                   paste0("[", round(quantile(precision, probs = c(0.025, 0.975))[1]*100, 1), " ; ",round(quantile(precision, probs = c(0.025, 0.975))[2]*100, 1),"]"),
                   paste0("[", round(quantile(recall, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(recall, probs = c(0.025, 0.975))[2]*100, 1),"]"),
                   paste0("[", round(quantile(f1_score, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(f1_score, probs = c(0.025, 0.975))[2]*100, 1),"]"),
                   paste0("[", round(quantile(roc_auc, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(roc_auc, probs = c(0.025, 0.975))[2]*100, 1),"]"),
                   paste0("[", round(quantile(specificity, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(specificity, probs = c(0.025, 0.975))[2]*100, 1),"]")
)

df_mean_sd$feature_extraction_method <- analysis_name

write.csv(df_mean_sd, paste0(output_path,"/", cancer_tissue, "_",omics_type,"_",tolower(analysis_name),"_satistics_table.csv"))

cat("\n")
end_time = Sys.time()

time_diff_str <- paste(as.character(round(end_time - start_time,2)), attr(end_time - start_time, "units"))
cat(bgRed(paste0("Execution time : ", time_diff_str)))
cat("\n")

# Save info file

file_name <- paste0(output_path,"/", cancer_tissue, "_",omics_type,"_",tolower(analysis_name),"_summary_file.txt")

writeLines("##SUMMARY FILE OF THE ANALYSIS##\n", file_name)

info_line <- c("\n1- Information about the analysis:",
               paste0("Cancer type: ", cancer_tissue),
               paste0("Omic analyzed: ", omics_type),
               paste0("Number of samples: ", nrow(X)),
               paste0("Number of features: ", ncol(X)),
               paste0("Execution time: ", time_diff_str))

write(info_line, file_name, append = T)

write("\n2- Mean and Standard deviation\n", file_name, append = T)
lines <- c(paste0("Mean Accuracy : ", round(mean(accuracy)*100,1), " (±", round(sd(accuracy)*100,1),")"),
           paste0("Mean Balanced accuracy : ", round(mean(balanced_accuracy)*100,1), " (±", round(sd(balanced_accuracy)*100,1),")"),
           paste0("Mean Precision : ", round(mean(precision)*100,1), " (±", round(sd(precision)*100,1),")"),
           paste0("Mean Recall : ", round(mean(recall)*100,1), " (±", round(sd(recall)*100,1),")"),
           paste0("Mean F1 score : ", round(mean(f1_score)*100,1), " (±", round(sd(f1_score)*100,1),")"),
           paste0("Mean ROC AUC : ", round(mean(roc_auc)*100,1), " (±", round(sd(roc_auc)*100,1),")"),
           paste0("Mean Specificity : ", round(mean(specificity)*100,1), " (±", round(sd(specificity)*100,1),")")
)
write(lines, file_name, append = T)

write("\n3- Confidence interval:\n", file_name, append = T)
lines_2 <- c(paste0("Accuracy 95% CI: [", round(quantile(accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]"),
             paste0("Balanced accuracy 95% CI: [", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(balanced_accuracy, probs = c(0.025, 0.975))[2]*100, 1),"]"),
             paste0("Precision 95% CI: [", round(quantile(precision, probs = c(0.025, 0.975))[1]*100, 1), " ; ",round(quantile(precision, probs = c(0.025, 0.975))[2]*100, 1),"]"),
             paste0("Recall 95% CI: [", round(quantile(recall, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(recall, probs = c(0.025, 0.975))[2]*100, 1),"]"),
             paste0("F1 score 95% CI: [", round(quantile(f1_score, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(f1_score, probs = c(0.025, 0.975))[2]*100, 1),"]"),
             paste0("ROC AUC 95% CI: [", round(quantile(roc_auc, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(roc_auc, probs = c(0.025, 0.975))[2]*100, 1),"]"),
             paste0("Specificity 95% CI: [", round(quantile(specificity, probs = c(0.025, 0.975))[1]*100, 1), " ; ", round(quantile(specificity, probs = c(0.025, 0.975))[2]*100, 1),"]")
)

write(lines_2, file_name, append = T)

cat(bgMagenta(paste0(Sys.time(),": End of the analysis")))
cat("\n")



