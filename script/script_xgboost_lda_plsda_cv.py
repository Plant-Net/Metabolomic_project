#Import packages
from email.policy import default
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#import xgboost
import xgboost as xgb
from xgboost import XGBClassifier, plot_importance

#import sklearn
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import RepeatedStratifiedKFold
#from sklearn.metrics import accuracy_score
from sklearn import metrics
from sklearn.metrics import *

import shap

import argparse
import os
import glob
import time
from datetime import datetime

import  warnings
warnings.filterwarnings('ignore')


#start_time = datetime.now()


def timer(start_time=None):
    if not start_time:
        start_time = datetime.now()
        return start_time
    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\n Time taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))


def xgboost_model(X_train, y_train, X_test):

    xgboost = XGBClassifier(seed=1234)

    xgboost.fit(X_train,y_train)
    # # Get predicted probability
    y_pred_proba = xgboost.predict_proba(X_test)[:,1]
    y_pred = xgboost.predict(X_test)

    return xgboost, y_pred, y_pred_proba


def pls_da(X_train,y_train, X_test):
    
    # Define the PLS object for binary classification
    plsda = PLSRegression(n_components=2)
    
    # Fit the training set
    plsda.fit(X_train, y_train)

    y_pred_proba = plsda.predict(X_test)[:,0]
    
    # Binary prediction on the test set, done with thresholding
    y_pred = (plsda.predict(X_test)[:,0] > 0.5).astype('uint8')
    
    return plsda, y_pred, y_pred_proba

def lda(X_train,y_train, X_test):
    
    # Define the model
    lda_model = LinearDiscriminantAnalysis()
    
    # Fit the training set
    lda_model.fit(X_train, y_train)

    y_pred_proba = lda_model.predict_proba(X_test)[:,1]
    
    # Binary prediction
    y_pred = lda_model.predict(X_test)
    
    return lda_model, y_pred, y_pred_proba



def shap_f(feat_data, fit_model):
    shap.initjs()
    
    explainer = shap.TreeExplainer(fit_model)
    shap_values = explainer.shap_values(feat_data) # this is a vector not a matrix
        
    figure = plt.figure()
    
    return shap.summary_plot(shap_values, feat_data, max_display=20, show = False, plot_size=[15,5])


if __name__ == '__main__':
    # Parse the command-line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='Path to the input file')
    
    parser.add_argument('-o', '--output_folder_path', type=str, required=True,
                        help='Path to the output folder containing the results')
    
    parser.add_argument('-a', '--analysis_name', type=str, required=True,
                        help='Name to the analysis. It can be extraction method name like : FA, GPLVM, HDDA, KPCA, MPPCA, PCA or NFE for No Feature Extraction.')
    
    parser.add_argument('-t', '--tissue_type', type=str, required=True,
                        help='Type of tissue analysed')
    
    parser.add_argument('-m', '--omics_type', type=str, required=True,
                        help='Type of omics')
    
    parser.add_argument('-n', '--fold_number', type=int, default = 4,
                        help='Number of fold for cross-validation.')
    
    parser.add_argument('-r', '--repetition_number', type=int, default=5,
                        help='Number of repetitions for cross validation.')
    
    
    args = parser.parse_args()

    # Get content of passed arguments
    outputpath = args.output_folder_path
    analysisname = args.analysis_name
    omics_type = args.omics_type
    n_splits = int(args.fold_number)
    n_repeats = int(args.repetition_number)
    tissue_type = args.tissue_type

    start_time = timer(None)

    # Load the input file
    if omics_type == 'metabolomics':
        df = pd.read_csv(args.input_file)
    
    else:
        df = pd.read_csv(args.input_file, index_col=0)

    #if analysisname.upper() == 'KPCA':
        #df['KPC3'] = df['KPC3'].astype(float)
        #df['KPC4'] = df['KPC3'].astype(float)

    X=df.drop("Label", axis=1)
    y=df["Label"]


    rkf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=0)

    # Create a dataframe in which all metrics of all iterations will be stored.
    df_metrics = pd.DataFrame(columns=['accuracy', 'balanced_accuracy', 'precision', 'recall', 'f1score',
                                    'roc_auc', 'specificity'],
                            index = ['fold' + str(x) for x in range(1,(n_splits*n_repeats +1))])


    for i, (train_index, test_index) in enumerate(rkf.split(X,y)):
        
        X_train=X.iloc[train_index]
        y_train=y.iloc[train_index]
        X_test=X.iloc[test_index] 
        y_test=y.iloc[test_index]

        if analysisname.upper() == "LDA":
            model, y_pred, y_pred_proba = lda(X_train, y_train, X_test)
            print("LDA")
        
        elif analysisname.upper() == "PLSDA":
            model, y_pred, y_pred_proba = pls_da(X_train, y_train, X_test)
            print("PLSDA")
        
        else:
            xgboost, y_pred, y_pred_proba = xgboost_model(X_train, y_train, X_test)
            print("XGBoost")


        cm = confusion_matrix(y_test, y_pred)
        print(cm)

        tn = cm[0, 0]
        fp = cm[0, 1]
        fn = cm[1, 0]
        tp = cm[1, 1]
        
        acc = accuracy_score(y_test, y_pred)
        bal_acc = balanced_accuracy_score(y_test, y_pred)
        f1score = metrics.f1_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        roc_auc = roc_auc_score(y_test, y_pred_proba)
        spe= tn / (tn + fp)
        
        # Fill in the dataframe with metrics
        df_metrics.iloc[i] = [acc, bal_acc, f1score, precision, recall, roc_auc, spe]


    # Save the dataframe of metrics
    df_metrics.to_csv(outputpath+'/'+tissue_type+'_'+omics_type+'_'+analysisname+"_metrics_table.csv")

    # MEAN AND SD
    df_stats = pd.DataFrame(columns=['accuracy', 'balanced_accuracy', 'precision', 'recall', 'f1score','roc_auc', 'specificity'],
                            index=['mean_sd', 'confidence_interval'])
    
    list_row = [f'{np.mean(df_metrics["accuracy"])*100:.1f} (±{np.std(df_metrics["accuracy"])*100:.1f})',
                f'{np.mean(df_metrics["balanced_accuracy"])*100:.1f} (±{np.std(df_metrics["balanced_accuracy"])*100:.1f})',
                f'{np.mean(df_metrics["precision"])*100:.1f} (±{np.std(df_metrics["precision"])*100:.1f})',
                f'{np.mean(df_metrics["recall"])*100:.1f} (±{np.std(df_metrics["recall"])*100:.1f})',
                f'{np.mean(df_metrics["f1score"])*100:.1f} (±{np.std(df_metrics["f1score"])*100:.1f})',
                f'{np.mean(df_metrics["roc_auc"])*100:.1f} (±{np.std(df_metrics["roc_auc"])*100:.1f})',
                f'{np.mean(df_metrics["specificity"])*100:.1f} (±{np.std(df_metrics["specificity"])*100:.1f})']
    
    df_stats.iloc[0]=list_row
    
    # CONFIDENCE INTERVAL
    acc_ci= np.percentile(df_metrics["accuracy"], [2.5,97.5])
    bal_acc_ci= np.percentile(df_metrics["balanced_accuracy"], [2.5,97.5])
    precision_ci = np.percentile(df_metrics["precision"], [2.5,97.5])
    recall_ci= np.percentile(df_metrics["recall"], [2.5,97.5])
    f1_ci = np.percentile(df_metrics["f1score"], [2.5,97.5])
    roc_auc_ci = np.percentile(df_metrics["roc_auc"], [2.5,97.5])
    specificity_ci = np.percentile(df_metrics['specificity'], [2.5,97.5])

    list_row = [f'[{acc_ci[0]:.3f} ; {acc_ci[1]:.3f}]',
                f'[{bal_acc_ci[0]:.3f} ; {bal_acc_ci[1]:.3f}]',
                f'[{precision_ci[0]:.3f} ; {precision_ci[1]:.3f}]',
                f'[{recall_ci[0]:.3f} ; {recall_ci[1]:.3f}]',
                f'[{f1_ci[0]:.3f} ; {f1_ci[1]:.3f}]',
                f'[{roc_auc_ci[0]:.3f} ; {roc_auc_ci[1]:.3f}]',
                f'[{specificity_ci[0]:.3f} ; {specificity_ci[1]:.3f}]'
                ]
    
    df_stats.iloc[1]=list_row

    df_stats.to_csv(outputpath+'/'+tissue_type+'_'+omics_type+'_'+analysisname+'_statistics_table.csv')

    end_time = datetime.now()
    print('Program executed in: {}'.format(end_time - start_time))

    # Create a file summary
    df=open(outputpath+'/'+tissue_type+'_'+omics_type+'_'+analysisname+'_summary_file.txt','w')
    df.write('1- Information about the analysis \n')
    df.write('Cancer type: '+tissue_type+'\n')
    df.write('Omic analyzed: '+omics_type+'\n')
    df.write('Number of samples: '+str(len(X.axes[0]))+'\n')
    df.write('Number of features: '+str(len(X.axes[1]))+'\n')
    df.write('Execution time: {} \n'.format(end_time - start_time))


    df.write("\n2- Mean and Standard deviation: \n")
    df.write(f'The mean accuracy is: {np.mean(df_metrics["accuracy"])*100:.1f} (±{np.std(df_metrics["accuracy"])*100:.1f}) \n')
    df.write(f'The mean balanced accuracy is: {np.mean(df_metrics["balanced_accuracy"])*100:.1f} (±{np.std(df_metrics["balanced_accuracy"])*100:.1f}) \n')
    df.write(f'The mean precision is: {np.mean(df_metrics["precision"])*100:.1f} (±{np.std(df_metrics["precision"])*100:.1f}) \n')
    df.write(f'The mean recall is: {np.mean(df_metrics["recall"])*100:.1f} (±{np.std(df_metrics["recall"])*100:.1f}) \n')
    df.write(f'The mean F1 score is: {np.mean(df_metrics["f1score"])*100:.1f} (±{np.std(df_metrics["f1score"])*100:.1f}) \n')
    df.write(f'The mean ROC AUC is: {np.mean(df_metrics["roc_auc"])*100:.1f} (±{np.std(df_metrics["roc_auc"])*100:.1f}) \n')
    df.write(f'The mean specificity is: {np.mean(df_metrics["specificity"])*100:.1f} (±{np.std(df_metrics["specificity"])*100:.1f}) \n')

    

    df.write("\n3- Confidence interval: \n")
    df.write(f'The confidence interval for accuracy is [{acc_ci[0]:.3f} ; {acc_ci[1]:.3f}]\n')
    df.write(f'The confidence interval for balanced accucary is [{bal_acc_ci[0]:.3f} ; {bal_acc_ci[1]:.3f}]\n')
    df.write(f'The confidence interval for precision is [{precision_ci[0]:.3f} ; {precision_ci[1]:.3f}]\n')
    df.write(f'The confidence interval for recall is [{recall_ci[0]:.3f} ; {recall_ci[1]:.3f}]\n')
    df.write(f'The confidence interval for F1 score is [{f1_ci[0]:.3f} ; {f1_ci[1]:.3f}]\n')
    df.write(f'The confidence interval for ROC AUC score is [{roc_auc_ci[0]:.3f} ; {roc_auc_ci[1]:.3f}]\n')
    df.write(f'The confidence interval for specificity score is [{specificity_ci[0]:.3f} ; {specificity_ci[1]:.3f}]\n')


    df.close()

    



    