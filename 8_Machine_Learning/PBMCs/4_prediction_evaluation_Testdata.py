# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia 
Date: 2025-10-01
Description:
This script is to evaluate the performance of the trained XGBoost model for PBMC and non-PBMC classification.
It loads the trained model, makes predictions on the test dataset, return the gating plots for each sample and computes performance metrics.

This needs AliGater package for gating plot generation.

"""

import aligater as ag
from sklearn.preprocessing import StandardScaler
import glob,joblib,os,torch
import numpy as np
from sklearn.metrics import (accuracy_score, f1_score, precision_score, recall_score, confusion_matrix, roc_auc_score, average_precision_score)


################################################# gating PBMC based on the trained model, adn return gate plot for each sample ##########################################################
# Load the trained XGBoost model from the output of PBMC_xgboost_training.py
model_path = "/home/wenxiaren/Part_3_ML/Training/PBMC/PBMC_xgboost_model.joblib"
clf = joblib.load(model_path)
xgb_model = clf["model"]
threshold = clf["threshold"]
# load tesdata path from fcs_path_train_test_split.py
test_data_path = "/home/wenxiaren/Part_3_ML/Data_Processing/test_path.cs"
# predict sample by sample
with open(test_data_path, "r") as f:
    for i in f:
        i = i.strip()
        fcs = ag.loadFCS(path=i, return_type="agsample",compensate=True, flourochrome_area_filter=True)
        if fcs is None:
            print("Skipping sample: fcs is None (filtered out or failed to load)")
            continue  
        base_path = "/home/wenxiaren/cbio4/data/Aitzkoa/ImmSight/"
        root_path = os.path.relpath(fcs.filePath,base_path)
        Image_ID_prefix_1 = root_path.replace("/","+")
        # Image_ID_prefix will be part of the name of the images
        Image_ID_prefix =Image_ID_prefix_1.replace(".fcs", "")
        metadata = fcs()
        features = ["FSC 488/10-A","SSC 488/10-A","FSC 488/10-H"]
        out =  metadata[features].copy()
        scaler = StandardScaler()
        scaler.fit(out[["FSC 488/10-A", "SSC 488/10-A", "FSC 488/10-H"]])
        out[["FSC_A_Normal", "SSC_A_Normal", "FSC_H_Normal"]] = scaler.transform(out[["FSC 488/10-A", "SSC 488/10-A", "FSC 488/10-H"]])
        out["FSC_FSC"] = out["FSC_A_Normal"] *out["FSC_H_Normal"]
        X_processed = np.column_stack([
                out["SSC 488/10-A"],
                out["FSC 488/10-A"],
                out["FSC 488/10-H"],
                out["FSC_A_Normal"],   # FSC_q 
                out["SSC_A_Normal"],   # SSC_q
                out["FSC_H_Normal"],   # FSC_H
                out["FSC_FSC"]
            ])  
        if calibrator is not None:
            y_probs = calibrator.predict_proba(X_processed)[:, 1]
        else:
            y_probs = xgb_model.predict_proba(X_processed)[:, 1]
    
        preds= (y_probs >= threshold).astype(int)
        out["pred_is_pbmc"] = preds
        pbmc_indices = out.index[out["pred_is_pbmc"] == 1].tolist()
        PBMC = ag.AGgate(pbmc_indices,None, "FSC 488/10-A","SSC 488/10-A", "PBMC")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/PBMC/tmp/"+Image_ID_prefix+"-PBMCs.png"
        # return PBMC gate plot
        ag.backGate(fcs,population=PBMC,background_population=None,xCol="FSC 488/10-A", yCol="SSC 488/10-A",markersize=0.1, filePlot=fileName)

################################################### evaluate the XGBoost model on test data ############################################################################3
def evaluate_xgboost(xgb_model, pt_files, threshold, calibrator=None):
    """
    Evaluate the XGBoost model on test .pt files and compute performance metrics.
    inputs:
        xgb_model: Trained XGBoost model.
        pt_files: List of paths to .pt files containing test data.
        threshold: Probability threshold for classification.
        calibrator: Optional probability calibrator.
    Returns:
        performance metrics.
        Confusion matrix.
    """
    
    y_true_all, y_pred_all, y_probs_all = [], [], []
    # make predictions for each .pt file or each sample
    for pt_file in pt_files:
        try:
            data = torch.load(pt_file, map_location="cpu")
            X, y_true, ID = data["X"], data["y"], data["ID"]
            if len(X) == 0:
                print(f" file {ID} is not valid,jumping...")
                continue
            if isinstance(X, torch.Tensor):
                X = X.cpu().numpy()
            if isinstance(y_true, torch.Tensor):
                y_true = y_true.cpu().numpy()

            # feature processing for the test data to fit the model
            keep_features =  ["SSC 488/10-A","FSC 488/10-A","FSC 488/10-H"]
            all_features =["SSC 488/10-A","FSC 488/10-A","FSC 488/10-H"]
            cols = [all_features.index(f) for f in keep_features]
            X_raw = X[:, cols] 
            scaler = StandardScaler()
            scaler.fit(X_raw)       
            X_scaled = scaler.transform(X_raw)
            FSC_FSC = X_raw[:, 0] * X_raw[:, 2]  # FSC-A * FSC-H
            X_processed = np.column_stack([
                X_raw[:, 0],
                X_raw[:, 1],
                X_raw[:, 2],
                X_scaled[:, 1],   # FSC_A_s
                X_scaled[:, 0],   # SSC_A_s
                X_scaled[:, 2],   # FSC_H_s"
                FSC_FSC
            ]) 
            y_binary =  y_true
            # make predictions
            if calibrator is not None:
                y_probs = calibrator.predict_proba(X_processed)[:, 1]
            else:
                y_probs = xgb_model.predict_proba(X_processed)[:, 1]
        
            y_pred = (y_probs >= threshold).astype(int)
            
            # save predictions
            y_true_all.extend(y_binary)
            y_pred_all.extend(y_pred)
            y_probs_all.extend(y_probs)
            
        except Exception as e:
            print(f"dealing with file {pt_file} error: {e}")
            import traceback
            traceback.print_exc()
            continue

    # calculate performance metrics across all test samples
    # convert to numpy arrays
    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)
    y_probs_all = np.array(y_probs_all)
    # calculate performance metrics
    acc = accuracy_score(y_true_all, y_pred_all)
    precision = precision_score(y_true_all, y_pred_all)
    recall = recall_score(y_true_all, y_pred_all)
    f1 = f1_score(y_true_all, y_pred_all)
    roc_auc = roc_auc_score(y_true_all, y_probs_all)
    pr_auc = average_precision_score(y_true_all, y_probs_all)
    # the confusion matrix
    cm = confusion_matrix(y_true_all, y_pred_all)
    tn, fp, fn, tp = cm.ravel()
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    # compile performance metrics
    performance = {
        "Accuracy": acc,
        "Precision": precision,
        "Recall": recall,
        "F1_Score": f1,
        "Specificity": specificity,
        "ROC_AUC": roc_auc,
        "PR_AUC": pr_auc,
        "Threshold": threshold,
        "Total_Samples": len(y_true_all),
        "Class_0_Count": int(np.sum(y_true_all == 0)),
        "Class_1_Count": int(np.sum(y_true_all == 1)),
        "TP": int(tp),
        "TN": int(tn),
        "FP": int(fp),
        "FN": int(fn)
    }
    return performance, cm

if __name__ == "__main__":

    # Load the trained XGBoost model from the output of PBMC_xgboost_training.py
    model_path = "/home/wenxiaren/Part_3_ML/Training/PBMC/PBMC_xgboost_model.joblib"
    try:
        PBMC_xgboost_model = joblib.load(model_path)
        xgb_model = PBMC_xgboost_model["model"]
        threshold = PBMC_xgboost_model["threshold"] 
    except FileNotFoundError:
        print(f" {model_path}")
        exit(1)
    except Exception as e:
        print(f" {e}")
        exit(1)

    # Load test .pt files from the test_data_pre.py output
    test_dir = "/home/wenxiaren/Part_3_ML/Dataset/PBMC_dataset/test_data"
    pt_files = glob.glob(f"{test_dir}/*.pt")
    # Evaluate model
    performance, cm = evaluate_xgboost(
        xgb_model=xgb_model,
        pt_files=pt_files,
        threshold=threshold,
        calibrator=None 
    )
    print(performance)
    print(cm)        



