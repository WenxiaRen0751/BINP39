# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import os
import random
import pandas as pd

"""
This script splits a list of flow cytometry data file paths into training and testing sets.

"""
# seperate the original fcs data pathinto train and test split
path = "/home/wenxiaren/Part_3_ML/Data_Processing/old_data_path_B_panel_final_1007.csv"   
with open(path) as f:
    all_files = [line.strip() for line in f if line.strip()]
random.shuffle(all_files)
random.seed(40)  
n = len(all_files)
train_n = int(0.9 * n)
train = all_files[:train_n]
test  = all_files[train_n:]
pd.Series(train).to_csv("/home/wenxiaren/Machine_Learning/Data_Processing/train_path.csv",
                        index=False, header=False)
pd.Series(test).to_csv("/home/wenxiaren/Machine_Learning/Data_Processing/test_path.csv",
                    index=False, header=False)
print(f"Total files: {n}, train = {len(train)}, test = {len(test)}")