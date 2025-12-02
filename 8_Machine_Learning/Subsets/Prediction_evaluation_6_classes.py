# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Author: Wenxia 
Date: 2025-10-11
Description:
This script is to evaluate the performance of the trained MPL model for 6 subsets classification from singlets
It loads the trained model, makes predictions on the test dataset, and computes performance metrics.

"""

import os, glob
import torch
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import (accuracy_score, f1_score, precision_score, recall_score, confusion_matrix,roc_auc_score, average_precision_score)
import torch.nn as nn
import torch.nn.functional as F

# loading the model
class MLP(nn.Module):
    def __init__(self, in_dim, hidden_dim1, hidden_dim2, out_dim):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(in_dim, hidden_dim1)
        self.bn1 = nn.BatchNorm1d(hidden_dim1)
        self.drop1 = nn.Dropout(0.3)

        self.fc2 = nn.Linear(hidden_dim1, hidden_dim2)
        self.bn2 = nn.BatchNorm1d(hidden_dim2)
        self.drop2 = nn.Dropout(0.3)

        self.fc3 = nn.Linear(hidden_dim2, out_dim)

    def forward(self, x):
        x = F.relu(self.bn1(self.fc1(x)))
        x = self.drop1(x)
        x = F.relu(self.bn2(self.fc2(x)))
        x = self.drop2(x)
        return self.fc3(x)

# map leaves → parent nodes (for hierarchical analysis)
def build_parent_mapping():
    """
    Build mapping from leaf classes to parent nodes in the hierarchy.
    """
    leaf_to_node = {
        0: "singlet/CD45+CD3-/CD19+",
        1: "singlet/CD45+CD3-/CD19-/CD14+/CM",
        2: "singlet/CD45+CD3-/CD19-/CD14+/non_CM",
        3: "singlet/CD45+CD3-/CD19-/CD14-",
        4: "singlet/CD45+CD3+",
        5: "singlet/rest_of_singlets"
    }
    parent_map = {}
    for leaf, node in leaf_to_node.items():
        parents = []
        parts = node.split("/")
        for i in range(1, len(parts)+1):
            parent = "/".join(parts[:i])
            parents.append(parent)
        parent_map[leaf] = parents

    return leaf_to_node, parent_map

# Evaluation
def evaluate(model, pt_files, leaf_to_name, num_classes,device="cpu", batch_size=8192):
    
    """
    Evaluate the MLP model on test data and calculate performance metrics.
    inputs:
        model: Trained MLP model
        pt_files: List of paths to .pt files containing test data.
        leaf_to_name: Mapping from leaf class indices to class names.
        num_classes: Number of classes. 6 classes here.
    outputs:
        df_all: DataFrame with each class performance metrics.
        df_parent: DataFrame with parent node performance metrics.

    """
    y_true_all, y_pred_all = [], []
    leaf_probs_all = []
    # Build hierarchical mappings for parent-level analysis
    leaf_to_node, parent_map = build_parent_mapping()
    parent_true_all = {}
    with torch.no_grad():
        # Process each .pt file for each sample
        for pt_file in pt_files:
            data = torch.load(pt_file, map_location="cpu")
            X, y_true, ID = data["X"].to(device), data["y"].to(device), data["ID"]
            # rename classes from 1-6 to 0-5
            y_true[y_true==2]=1
            y_true[y_true==3]=2
            y_true[y_true==4]=3
            y_true[y_true==5]=4
            y_true[y_true==6]=5
            features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A","CD14","CD19","CD3","CD45","CD16",'BV786 CD14-A_v2',
                'BV605 CD19-A_v2',
                'Alexa 700 CD3-A_v2',
                'APC-H7 CD45-A_v2',
                'APC CD16-A_v2']
            expected_features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A",'BV786 CD14-A_v2',
                'BV605 CD19-A_v2',
                'Alexa 700 CD3-A_v2',
                'APC-H7 CD45-A_v2',
                'APC CD16-A_v2']
            feature_indices = [features.index(feat) for feat in expected_features]
            X = X[:, feature_indices] 
            cd14_idx = expected_features.index('BV786 CD14-A')
            cd16_idx = expected_features.index('APC CD16-A')
            # extract columns 
            cd14_col = X[:, cd14_idx]
            cd16_col = X[:, cd16_idx]
            # new features. the ratio of CD14/CD16
            CD14_by_16 = (cd16_col / (cd14_col + 1e-8)).reshape(-1, 1)  
            X = torch.cat((X,CD14_by_16),1)
            expected_features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A",'BV786 CD14-A_v2',
            'BV605 CD19-A_v2','Alexa 700 CD3-A_v2','APC-H7 CD45-A_v2','APC CD16-A_v2',"CD14_by_16"]
            if len(X) == 0:
                print(f"No valid samples in {ID}, skipping...")
                continue
            # prediction
            probs_all, pred_all = [], []
            # Process in batches
            for i in range(0, len(X), batch_size):
                xb = X[i:i+batch_size]
                logits = model(xb)
                probs = torch.softmax(logits, dim=1) 
                pred = torch.argmax(probs, dim=1)     
                probs_all.append(probs)
                pred_all.append(pred)
            # concatenate batch results
            probs = torch.cat(probs_all, dim=0)   
            pred = torch.cat(pred_all, dim=0)     
            y_true_np = y_true.cpu().numpy()
            y_pred_np = pred.cpu().numpy()
            probs_np = probs.cpu().numpy()
            # combine results
            y_true_all.extend(y_true_np)
            y_pred_all.extend(y_pred_np)
            leaf_probs_all.append(probs_np)
            # Parent-level accuracy analysis
            for t, p in zip(y_true_np, y_pred_np):
                if t in parent_map and p in parent_map:
                    true_parents = parent_map[t]
                    pred_parents = parent_map[p]
                    for parent in true_parents:
                        parent_true_all.setdefault(parent, []).append(1 if parent in pred_parents else 0)

    y_true_all = np.array(y_true_all, dtype=int)
    y_pred_all = np.array(y_pred_all, dtype=int)
    leaf_probs_all = np.concatenate(leaf_probs_all, axis=0) 

    # calculate overall metrics
    cm = confusion_matrix(y_true_all, y_pred_all, labels=list(range(num_classes)))
    tp = np.diag(cm).astype(float)
    support = cm.sum(1).astype(float)
    pred_total = cm.sum(0).astype(float)
    FN = support - tp
    FP = pred_total - tp
    TN = cm.sum() - (tp + FP + FN)
    precision = tp / np.maximum(pred_total, 1.0)
    recall = tp / np.maximum(support, 1.0)
    specificity = TN / np.maximum(TN + FP, 1.0)
    f1 = 2 * precision * recall / np.maximum(precision + recall, 1e-12)

    # AUC for each class
    perclass_auc = []
    for k in range(num_classes):
        y_k = (y_true_all == k).astype(int)
        try:
            auc_k = roc_auc_score(y_k, leaf_probs_all[:, k])
        except:
            auc_k = np.nan
        perclass_auc.append(auc_k)

    # save class-level metrics
    df_all = pd.DataFrame({
        "Class ID": np.arange(num_classes),
        "Precision": precision,
        "Recall": recall,
        "Specificity": specificity,
        "Accuracy": recall,   # same as class-level accuracy
        "F1 Score": f1,
        "AUC": perclass_auc
    })
    # calculate overall averages
    df_all.loc["Overall"] = [
        "", precision.mean(), recall.mean(),
        specificity.mean(), recall.mean(),
        f1.mean(), np.nanmean(perclass_auc)
    ]
    # save 
    df_all.to_csv("classes_metrics.csv")
 
    # Parent performance
    parent_metrics = []
    for parent, vals in parent_true_all.items():
        vals = np.array(vals)
        acc_parent = vals.mean() if len(vals) > 0 else np.nan
        parent_metrics.append({"Node": parent, "Support": int(len(vals)), "Accuracy": acc_parent})
    df_parent = pd.DataFrame(parent_metrics)
    # Parent accuracy heatmap 
    parent_name = df_parent["Node"].tolist()
    mat = np.zeros((len(parent_name), len(parent_name)))
    for i, row in df_parent.iterrows():
        if not np.isnan(row["Accuracy"]):
            mat[i, i] = row["Accuracy"] * 100
    plt.figure(figsize=(10, 8))
    sns.heatmap(mat, annot=True, fmt=".1f", cmap="YlGnBu",
                xticklabels=parent_name, yticklabels=parent_name,
                cbar_kws={'label': 'Accuracy %'})
    plt.title("Parent Node Accuracy")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig("./accuracy_all.png", dpi=600, bbox_inches='tight')
    plt.close()

    return df_all, df_parent

if __name__ == "__main__":
    device = "cuda" if torch.cuda.is_available() else "cpu"
    # Load the model for 6 classes
    model_path = "/home/wenxia/Part_3_ML/training/mlp_6classes_2.pt"  
    checkpoint = torch.load(model_path, map_location=device)
    features = checkpoint["features"]
    leaf_to_name = checkpoint["leaf_to_name"]
    num_classes = checkpoint["num_classes"]
    # Initialize model
    model = MLP(in_dim=len(features), hidden_dim1=128, hidden_dim2=64, out_dim=num_classes).to(device)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    # load test files
    test_dir = "/home/wenxia/Part_3_ML/Dataset/test_data"
    pt_files = glob.glob(f"{test_dir}/*.pt")
    # evaluation
    df_all, df_parent = evaluate(model, pt_files, leaf_to_name, num_classes, device=device)
    print("class-level metrics saved to: classes_metrics.csv")
    print("Parent node accuracy saved to: accuracy_all.png")





  
