# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia 
Date: 2025-10-11
Description:
This script is to train an MLP model for 6 subsets classification of flow cytometry data.

"""

import torch,os
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, random_split

# earlyStopping
class EarlyStopping:
    def __init__(self, patience=10, min_delta=1e-4):
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.best_loss = float("inf")
        self.early_stop = False
    def __call__(self, val_loss):
        if val_loss < self.best_loss - self.min_delta:
            self.best_loss = val_loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True

# dataset
class FlowCytometryDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)
    def __len__(self):
        return len(self.y)
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

# define MLP
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
    
# Train
def MLP_training(X, y, features,out_path="mlp6.pt", epochs=25):

    """
    Train a 3-layer MLP model for 6-class flow cytometry classification.
    """
    device = "cuda" if torch.cuda.is_available() else "cpu"
    leaf_to_name = {
        0: "singlet/CD45posCD3neg/CD19pos",
        1: "singlet/CD45posCD3neg/CD19neg/CD14pos/CM",
        2: "singlet/CD45posCD3neg/CD19neg/CD14pos/non_CM",
        3: "singlet/CD45posCD3neg/CD19neg/CD14neg",
        4: "singlet/CD45posCD3pos",
        5: "singlet/rest_of_singlets"
    }
    num_classes = len(leaf_to_name)  # 7 classes (0-6)
    # Calculate class weights for imbalanced dataset
    class_counts = np.bincount(y, minlength=num_classes)
    class_weights = class_counts.max() / (class_counts + 1e-6)
    class_weights = class_weights / class_weights.sum() * len(class_weights)
    class_weights = torch.tensor(class_weights, dtype=torch.float32).to(device)
    # Create dataset
    dataset = FlowCytometryDataset(X, y)
    n_val = int(0.1 * len(dataset))
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val])
    train_loader = DataLoader(train_set, batch_size=2048, shuffle=True, num_workers=16, pin_memory=True)
    val_loader = DataLoader(val_set, batch_size=2048, shuffle=False, num_workers=16, pin_memory=True)
    # Model, loss & optimizer
    model = MLP(in_dim=len(features), hidden_dim1=128, hidden_dim2=64, out_dim=num_classes).to(device)
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=15, gamma=0.5)
    early_stopper = EarlyStopping(patience=8, min_delta=1e-4)
    best_val_loss = float("inf")
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            loss = criterion(logits, yb)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        train_loss = total_loss / len(train_loader)
        # Validation
        model.eval()
        val_loss = 0
        correct = 0
        total = 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(device), yb.to(device)
                logits = model(xb)
                loss = criterion(logits, yb)
                val_loss += loss.item()
                # accuracy
                pred = torch.argmax(logits, dim=1)
                correct += (pred == yb).sum().item()
                total += yb.size(0)
        val_loss /= len(val_loader)
        val_acc = correct / total
        scheduler.step()
   
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save({
                "model_state": model.state_dict(),
                "leaf_to_name": leaf_to_name,
                "features": features,
                "num_classes": num_classes,
                "class_weights": class_weights.cpu()
            }, out_path)

        early_stopper(val_loss)
        if early_stopper.early_stop:
            print(f"Early stopping at epoch {epoch+1}")
            break
    return model, leaf_to_name

# training
features =["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A","CD14","CD19","CD3","CD45","CD16",'BV786 CD14-A_v2', 'BV605 CD19-A_v2','Alexa 700 CD3-A_v2', 'APC-H7 CD45-A_v2', 'APC CD16-A_v2']
expected_features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A",'BV786 CD14-A_v2', 'BV605 CD19-A_v2','Alexa 700 CD3-A_v2', 'APC-H7 CD45-A_v2', 'APC CD16-A_v2']
# Load data
train_data_path = "/home/wenxia/Part_3_ML/Dataset/train_dataset_6classes.pt"
data = torch.load(train_data_path, map_location="cpu")
X = data["X"].cpu().numpy()   # shape (N, 15)
y = data["y"].cpu().numpy()   # shape (N,)
feature_indices = [features.index(feat) for feat in expected_features]
X = X[:, feature_indices] 
# add new feature to seperate the non-CM
cd14_idx = expected_features.index('BV786 CD14-A')
cd16_idx = expected_features.index('APC CD16-A')
# Extract columns (these will be 1D)
cd14_col = X[:, cd14_idx]
cd16_col = X[:, cd16_idx]
# new feature
CD14_by_16 = (cd16_col / (cd14_col + 1e-8)).reshape(-1, 1)
expected_features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A",'BV786 CD14-A_v2','BV605 CD19-A_v2', 'Alexa 700 CD3-A_v2', 'APC-H7 CD45-A_v2', 'APC CD16-A_v2',"CD14_by_16"]
# Train 
MLP_training(X, y,expected_features, out_path="mlp_6classes.pt", epochs=25)
