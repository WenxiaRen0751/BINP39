# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia
Date: 2025_10-01
Description:
This script is to train an XGBoost model for PBMC and non-PBMC classification using the training dataset generated from PBMC_traindata_preprocess.py.
The trained model and the best threshold based on F1 score will be saved for future prediction.

"""

# Import packages
import numpy as np
import pandas as pd
import joblib,torch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (classification_report, confusion_matrix, f1_score, roc_auc_score,precision_recall_curve, auc)
from xgboost import XGBClassifier

# define feature names
RAW_FEATURE_NAMES = ["SSC 488/10-A","FSC 488/10-A","FSC 488/10-H","FSC_A_Normal", "SSC_A_Normal", "FSC_H_Normal"]
MODEL_FEATURES =["SSC 488/10-A","FSC 488/10-A","FSC 488/10-H","FSC_A_Normal", "SSC_A_Normal", "FSC_H_Normal","FSC_FSC"]
# feature transformation, return processed training data with selected features
def build_features(X_raw):
    df = pd.DataFrame(X_raw, columns=RAW_FEATURE_NAMES)
    df["FSC_FSC"] = df["FSC 488/10-A"] * df["FSC 488/10-H"]
    X_model = df[MODEL_FEATURES].values
    return X_model

# find the best threshold based on F1 score
def find_best_threshold(y_true, probs):
    best_f1 = 0
    best_thresh = 0.5
    for thresh in np.linspace(0.01, 0.99, 99):
        pred = (probs >= thresh).astype(int)
        f1 = f1_score(y_true, pred)
        if f1 > best_f1:
            best_f1 = f1
            best_thresh = thresh
    return best_thresh, best_f1

def evaluate_with_kfold(X, y, n_splits=5, random_state=11):

    """
    Evaluate XGBoost model with K-Fold cross-validation
    inputs:
        X: feature matrix
        y: binary labels
        n_splits : number of folds for cross-validation
        random_state: Random seed for reproducibility
    Returns:
        k_fold_records: List containing performance metrics for each fold       
    """

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    k_fold_records = []
    # 5-Fold cross-validation
    for fold, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
        print(f"\nFold {fold}/{n_splits}")
        X_train_fold, X_val_fold = X[train_idx], X[val_idx]
        y_train_fold, y_val_fold = y[train_idx], y[val_idx]
        # train XGBoost model
        clf = XGBClassifier(
            n_estimators=200,
            learning_rate=0.1,
            max_depth=4,
            subsample=0.8,
            colsample_bytree=1.0,
            tree_method="hist",
            reg_alpha=0.1,
            reg_lambda=1.0,
            random_state=random_state,
            use_label_encoder=False,
            eval_metric="logloss",
            n_jobs=-1
        )
        # fit model
        clf.fit(X_train_fold, y_train_fold)
        # validate model
        probs = clf.predict_proba(X_val_fold)[:, 1]
        # determine best threshold
        best_thresh, best_f1 = find_best_threshold(y_val_fold, probs)
        # make predictions
        preds = (probs >= best_thresh).astype(int)
        # calculate metrics
        roc_auc = roc_auc_score(y_val_fold, probs)
        precision, recall, _ = precision_recall_curve(y_val_fold, probs)
        pr_auc = auc(recall, precision)
        f1 = f1_score(y_val_fold, preds)
        # record results
        k_fold_records.append({
            "fold": fold,
            "roc_auc": roc_auc,
            "pr_auc": pr_auc,
            "f1": f1
        })
        print(f"Fold {fold} - ROC AUC: {roc_auc:.4f}, PR AUC: {pr_auc:.4f}, F1: {f1:.4f}, Best Threshold: {best_thresh:.3f}")

    return k_fold_records


def model_training():

    """
    Train XGBoost model for PBMC vs non-PBMC classification on the full training dataset
    Returns: model_package: trained model, best threshold, and feature names 

    """
    # Load full training data generated from PBMC_traindata_preprocess.py
    train_data_path = "/home/wenxiaren/Part_3_ML/Dataset/PBMC_dataset/train_data/train_dataset.pt"
    data = torch.load(train_data_path, map_location="cpu")
    X = data["X"].cpu().numpy()   
    y = data["y"].cpu().numpy()  
    # binary labels
    y_binary = y.astype(int)
    # Build features
    X_processed= build_features(X)
    # Evaluate with 5-Fold cross-validation 
    kfold_results = evaluate_with_kfold(
        X_processed, y_binary, n_splits=5, random_state=11
    )
    # Train final model on the full train dataset
    final_model = XGBClassifier(
        n_estimators=200,
        learning_rate=0.1,
        max_depth=4,
        subsample=0.8,
        colsample_bytree=1.0,
        tree_method="hist",
        reg_alpha=0.1,
        reg_lambda=1.0,
        random_state=11,
        use_label_encoder=False,
        eval_metric="logloss",
        n_jobs=-1
    )
    # Fit model
    final_model.fit(X_processed, y_binary)
    # Determine best threshold on full data
    final_probs = final_model.predict_proba(X_processed)[:, 1]
    best_thresh, _ = find_best_threshold(y_binary, final_probs)
    final_preds = (final_probs >= best_thresh).astype(int)
    print(f"Best threshold on full data: {best_thresh:.3f}")
    # Package model
    model_package = {
        "model": final_model,
        "threshold": best_thresh,
        "feature_names": MODEL_FEATURES
    }
    # save model 
    joblib.dump(model_package, "PBMC_xgboost_model.joblib")
    print(f"\n=== Final Model Performance (threshold={best_thresh:.3f}) ===")
    print("Classification Report:")
    print(classification_report(y_binary, final_preds, digits=4))
    print("\nConfusion Matrix:")
    print(confusion_matrix(y_binary, final_preds))
    # return model 
    return model_package

if __name__ == "__main__":
    model_package = model_training()