# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia 
Date: 2025-10-11
Description:
This script is to evaluate the performance of the trained MPL model for 6 subsets classification from singlets
It loads the trained model, makes predictions on the test dataset, return the gating plots for each sample 
This needs AliGater package for gating plot generation.

"""
import aligater as ag
from math import inf
import numpy as np
import pandas as pd
import os, glob,yaml,torch
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler

# loading the model for 7 Classes classification
model_path = "/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/mlp_7classes.pt"  
device = "cuda" if torch.cuda.is_available() else "cpu"  
checkpoint = torch.load(model_path, map_location=device)
features = checkpoint["features"]
leaf_to_name = checkpoint["leaf_to_name"]
num_classes = checkpoint["num_classes"]
# define the MLP
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
model = MLP(in_dim=len(features), hidden_dim1=128, hidden_dim2=64, out_dim=num_classes).to(device)
model.load_state_dict(checkpoint["model_state"])
model.eval()

################################################# gating based on the trained model, adn return gate plot for each sample #########################################################
SSC="SSC 488/10-A"
FSC="FSC 488/10-A"
IgA="FITC IgA-A"
IgD="PerCP-Cy5.5 IgD-A"
CD141="PE-Cy7 CD141-A"
CD123="PE-Cy5 CD123-A"
CD27="PE (R-phycoerythrin) CD27-A"
CD34="PE-Dazzle 594 CD34-A"
HLADR="BV650 HLA-DR-A"
CD56="BV711 CD56-A"
CD14="BV786 CD14-A"
CD38="BV421 CD38-A"
CD19="BV605 CD19-A"
CD24="BV510 CD24-A"
CD3="Alexa 700 CD3-A"
CD45="APC-H7 CD45-A"
CD16="APC CD16-A"

# Read the corrections yaml file
corrections_yaml_file ="/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_B_Panel_1006/Scripts/correction.yaml"
with open(corrections_yaml_file, 'r') as stream:
    corrections_dict_raw = yaml.safe_load(stream)
    correction_1 =  corrections_dict_raw['PBMC Gating']['CUSTOM_100000']
    correction_2 = corrections_dict_raw['PBMC Gating']['CUSTOM_95000']
    correction_3 = corrections_dict_raw['PBMC Gating']['CUSTOM_90000']
    correction_4 = corrections_dict_raw['PBMC Gating']['Ylim_100K']
    correction_5 = corrections_dict_raw['PBMC Gating']['Ylim_110K']
    correction_6 = corrections_dict_raw['PBMC Gating']['Ylim_130K']
    correction_7 = corrections_dict_raw['PBMC Gating']['Ylim_90K']
    correction_8 = corrections_dict_raw['PBMC Gating']['xlim_65K']
    correction_9 = corrections_dict_raw['PBMC Gating']['xlim_70K']
    correction_10 = corrections_dict_raw['PBMC Gating']['xlim_80K']
    correction_11 = corrections_dict_raw['PBMC Gating']['No_clutter_240K']
    correction_12 = corrections_dict_raw['PBMC Gating']['Ylim_75K']
    correction_13 = corrections_dict_raw['PBMC Gating']['Dijkstra failure']
    correction_14 = corrections_dict_raw['Singlet Gating']['adjustAngle_10']
    correction_15 = corrections_dict_raw['Singlet Gating']['adjustAngle_5']
    correction_16 = corrections_dict_raw['Singlet Gating']['heightScale_3']
    correction_17 = corrections_dict_raw['Singlet Gating']['CUSTOM_2']
    correction_18 = corrections_dict_raw['Singlet Gating']['widthScale_4']
    correction_19 = corrections_dict_raw['Singlet Gating']['adjustAngle_14']
    correction_20 = corrections_dict_raw['Singlet Gating']['heightScale_1']
    correction_21 = corrections_dict_raw['Singlet Gating']['widthScale_6']
    correction_22 = corrections_dict_raw['CD45_CD19 Gating']['CUSTOM_x_1900']
    correction_23 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+100']
    correction_24 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+200']
    correction_25 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+500']
    correction_26 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+1000']
    correction_27 = corrections_dict_raw['CD45_CD19 Gating']['xlim_-100']
    correction_28 = corrections_dict_raw['CD45_CD19 Gating']['CUSTOM_y_to_2500']
    correction_29 = corrections_dict_raw['CD45_CD19 Gating']['ylim_-100']
    correction_30 = corrections_dict_raw['CD45_CD19 Gating']['ylim_-200']


test_data_path = "/home/wenxiaren/Part_3_ML/Data_Processing/test_path.csv"
# gating PBMC based on the trained model
with open(test_data_path, "r") as f:
    for i in f:
        i = i.strip()
        fcs = ag.loadFCS(path=i, return_type="agsample",compensate=True, flourochrome_area_filter=True)
        if fcs is None:
            print("Skipping sample: fcs is None (filtered out or failed to load)")
            continue  
        ag.agconf.ag_verbose=False
        date_plate = ag.getFileName(ag.getParent(fcs.filePath))
        date = date_plate.split("B panel")[0]
        plate = date_plate.split("B panel")[1]
        sampleName=ag.getFileName(fcs.filePath)
        # sample_ids are used to match wtih the id in the yaml file
        sample_ids = date+"-"+plate+"-"+sampleName
        if  sample_ids in correction_11:
            no_clutter_thresh = 240000
        else:
            no_clutter_thresh = 214000
        no_clutter1=ag.gateThreshold(fcs,"no_clutter", "FSC 488/10-A", "FSC 488/10-H",thresh= no_clutter_thresh, orientation='vertical',population="lower")
        no_clutter=ag.gateThreshold(fcs,"no_clutter","FSC 488/10-A", "FSC 488/10-H", parentGate=no_clutter1,thresh=214000, orientation='horizontal',population="lower")
        halfcut_middle = ag.gateThreshold(fcs, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=110000,
                                        parentGate=no_clutter, orientation='vertical', population='upper')
        ylim_bot = ag.densityDelimitation(fcs, xCol='SSC 488/10-A', parentGate=halfcut_middle, interval=[10000,20000], limit_threshold=0.05, direction='left',scale='linear')
        if ylim_bot == inf:
            ylim_bot = 10000

        halfcut_tail = ag.gateThreshold(fcs, name="right_tail",  xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=180000, parentGate=no_clutter, orientation='vertical', population='upper')
        ylim_top = ag.densityDelimitation(fcs, xCol='SSC 488/10-A', parentGate=halfcut_tail, interval=[50000, 125000], limit_threshold=0.2, direction='right',scale='linear')+25000
        if ylim_top == inf:
            ylim_top = 150000

        PBMC_step1 = ag.gateCorner(fcs, name="cut_corner", xCol='FSC 488/10-A', yCol='SSC 488/10-A', xThresh=65000, yThresh=ylim_bot+20000, 
                                xOrientation='lower', yOrientation='lower', Outer=True, parentGate=no_clutter)
        
        if  sample_ids in correction_13:
            startY_1= 30000
            xboundaries_1 = [70000,214000]
            direction_1 = 'up'
        else:
            startY_1= ylim_bot+10000
            xboundaries_1 = None
            direction_1 = 'both'

        PBMC_step2 = ag.horizontalPath(fcs, name="PBMC",
                        xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                        startY=startY_1, endY=ylim_top, xboundaries=xboundaries_1,
                        yboundaries=[ylim_bot-5000,ylim_top+5000],
                        leftRight=True , direction=direction_1,
                        maxStep=2, phi=0.1, bins=100, sigma=1,
                        scale='linear', parentGate=PBMC_step1)
        if  sample_ids in correction_8:
            cut_off_xlim = 65000
        elif sample_ids in correction_9:
            cut_off_xlim = 70000
        elif sample_ids in correction_10:
            cut_off_xlim = 80000
        elif sample_ids in correction_1:
            cut_off_xlim = 100000
        elif sample_ids in correction_2: 
            cut_off_xlim = 95000
        elif sample_ids in correction_3:
            cut_off_xlim = 90000
        else:
            cut_off_xlim = 75000
        outline = ag.gateThreshold(fcs, name="outline", xCol=FSC, yCol=SSC, thresh=cut_off_xlim,parentGate=PBMC_step2,
                                        orientation='vertical', population="upper",filePlot=None)
        if  sample_ids in correction_4:
            cut_off_ylim = 100000
        elif sample_ids in correction_5:
            cut_off_ylim = 110000
        elif sample_ids in correction_6:
            cut_off_ylim = 130000
        elif sample_ids in correction_7:
            cut_off_ylim = 90000
        elif sample_ids in correction_12:  
            cut_off_ylim = 75000
        else:
            cut_off_ylim = 120000
        outline_1 = ag.gateThreshold(fcs, name="outline_1", xCol=FSC, yCol=SSC, thresh=cut_off_ylim,parentGate=outline,
                                        orientation='horizontal', population="lower",filePlot=None)
        PBMC=ag.gatePC(fcs,name="PBMC",xCol=FSC,yCol=SSC,center='centroid', adjustAngle=2, widthScale=8, heightScale = 6, parentGate=outline_1,filePlot=None)
        ag.backGate(fcs,population=PBMC,background_population=None,xCol=FSC, yCol=SSC,markersize=0.1, filePlot=None)
        fcs.update(ag.AGgate(PBMC,None,FSC,SSC,"PBMC"), QC=True, xlim=[0,215000], ylim=[0,215000],MFI=True, MFI_type='all', extra_MFI=None)
        # Next: Singlet gating
        fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/singlets/"+date+"-"+plate+"-"+sampleName+"-singlets.png"
        # Setting the Angle
        if sample_ids in correction_14:
            adjustAngle_1 = 10
        elif sample_ids in correction_15:
            adjustAngle_1 = 5
        elif sample_ids in correction_17:
            adjustAngle_1 = 4
        elif sample_ids in correction_19:
            adjustAngle_1 = 14
        else:
            adjustAngle_1 = 1
        # Setting the height
        if sample_ids in correction_16 or sample_ids in correction_17:
            heightScale_1 = 3
        elif sample_ids in correction_20: 
            heightScale_1 = 0.5
        else:
            heightScale_1 = 4 
        # Setting width
        if sample_ids in correction_18:
            widthScale_1 = 4
        elif sample_ids in correction_21:
            widthScale_1 = 6
        else:
            widthScale_1 = 8
        singlets=ag.gatePC(fcs,xCol=FSC, yCol="FSC 488/10-H",name="singlets",center='density', adjustAngle= adjustAngle_1, widthScale=widthScale_1, 
                           heightScale=heightScale_1,parentGate=PBMC,filePlot=None)
        fcs.update(ag.AGgate(singlets,PBMC,FSC,"FSC 488/10-H","singlets"), QC=True, xlim=[25000,214000], ylim=[25000,214000],MFI=True, MFI_type='all', extra_MFI=None)
        
        ############################################################ # from here to do the prediction for 6 classes #########################################################
        metadata = fcs()
        features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A"]
        df_all =  metadata[features].copy()
        scaler = StandardScaler()
        scaler.fit(df_all[["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A"]])
        df_all[["CD14","CD19","CD3","CD45","CD16"]] = scaler.transform(df_all[["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A"]])
        features =["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A","CD14","CD19","CD3","CD45","CD16",]
        
        logtrans_dict = {"BV786 CD14-A":"BV786 CD14-A_v2","BV605 CD19-A":"BV605 CD19-A_v2","Alexa 700 CD3-A":"Alexa 700 CD3-A_v2","APC-H7 CD45-A":"APC-H7 CD45-A_v2","APC CD16-A":"APC CD16-A_v2"}
        for col, new_col in logtrans_dict.items():
            if col == "BV605 CD19-A":
                df_all[new_col] = [ag.transformWrapper(x, 1000, "bilog") for x in  df_all[col]]
            else:
                df_all[new_col] = [ag.transformWrapper(x, 500, "bilog") for x in  df_all[col]]
        df_all["CD14_by_16"] = df_all["APC CD16-A"]/df_all["BV786 CD14-A"]
        singlets_indices = singlets()
        expected_features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A","BV786 CD14-A_v2",
             "BV605 CD19-A_v2",
             "Alexa 700 CD3-A_v2",
             "APC-H7 CD45-A_v2",
             "APC CD16-A_v2","CD14_by_16"]
        out = df_all.loc[singlets_indices, expected_features].copy() 
        X_processed = out.to_numpy()
        X_input = torch.tensor(X_processed, dtype=torch.float32)
        # starting predicting new dataset
        logits = model(X_input)
        probs = torch.softmax(logits, dim=1)
        pred = torch.argmax(probs, dim=1)
        # assigning the prediction to the column "label"
        out["label"] = pred.numpy()
        # class mapping
        class_dict = {"CD19pos":0,"CM":1,"Inter_Mono":2,"non_CM":3,"CD14neg":4,"CD45posCD3pos":5,"rest_of_singlets":6}
        indices_dict = {}
        # extracting the class index
        for class_name, class_id in class_dict.items():
            indices = out.index[out["label"] == class_id].tolist()
            indices_dict[class_name] = indices
            # the output for this one is like: non_PBMCs = indices_dict["non_PBMCs"]
        total_event = len(metadata)
        total_indexs = np.arange(total_event)
        # plotting the CD45+CD3+
        CD45posCD3pos_indices = indices_dict["CD45posCD3pos"]           
        CD45posCD3pos = ag.AGgate(CD45posCD3pos_indices,singlets,"Alexa 700 CD3-A","APC-H7 CD45-A","CD45pos")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-CD45posCD3pos.png"
        ag.backGate(fcs,population=CD45posCD3pos,background_population=singlets,xCol="Alexa 700 CD3-A", yCol="APC-H7 CD45-A",markersize=0.1, scale='bilog', T=500,filePlot=fileName)
        fcs.update(ag.AGgate(CD45posCD3pos,singlets,CD3,CD45,"CD45pos"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
        # plotting the CD45+CD3- 
        CD45posCD3neg_indices =  list(np.concatenate([indices_dict["CD19pos"], indices_dict["CM"],indices_dict["Inter_Mono"],indices_dict["CD14neg"],indices_dict["non_CM"]]))
        CD45posCD3neg = ag.AGgate(CD45posCD3neg_indices ,singlets,"Alexa 700 CD3-A","APC-H7 CD45-A","CD45pos")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-CD45posCD3neg.png"
        ag.backGate(fcs,population=CD45posCD3neg,background_population=singlets,xCol="Alexa 700 CD3-A", yCol="APC-H7 CD45-A",markersize=0.1, scale='bilog', T=500,filePlot=fileName)
        fcs.update(ag.AGgate(CD45posCD3neg,singlets,CD3,CD45,"CD45posCD3pos"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
        # plotting the CD19+
        CD19pos = ag.AGgate(indices_dict["CD19pos"],CD45posCD3neg,"BV605 CD19-A","APC-H7 CD45-A","CD19pos")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-CD19pos.png"
        ag.backGate(fcs,population=CD19pos,background_population=CD45posCD3neg,xCol="BV605 CD19-A", yCol="APC-H7 CD45-A",markersize=0.2, scale='bilog', T=1000,filePlot=fileName)
        fcs.update(ag.AGgate(CD19pos,CD45posCD3neg,CD19,CD45,"CD19pos"), QC=True, xlim=[0,20000], ylim=[0,20000],scale='bilog', T=1000,MFI=True, MFI_type='all', extra_MFI=None)
        CD19neg_indices =  list(np.concatenate([indices_dict["CM"],indices_dict["Inter_Mono"],indices_dict["CD14neg"],indices_dict["non_CM"]]))
        # plotting the CD19-
        CD19neg = ag.AGgate(CD19neg_indices,CD45posCD3neg,"BV605 CD19-A","APC-H7 CD45-A","CD19pos")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-CD19neg.png"
        ag.backGate(fcs,population=CD19neg,background_population=CD45posCD3neg,xCol="BV605 CD19-A", yCol="APC-H7 CD45-A",markersize=0.1, scale='bilog', T=1000,filePlot=fileName)
        fcs.update(ag.AGgate(CD19neg,CD45posCD3neg,CD19,CD45,"CD19neg"), QC=True, xlim=[0,20000], ylim=[0,20000],scale='bilog', T=1000,MFI=True, MFI_type='all', extra_MFI=None)
        # plotting the CM
        CM = ag.AGgate(indices_dict["CM"],CD19neg,"BV786 CD14-A","APC CD16-A","CM")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-CM.png"
        ag.backGate(fcs,population=CM,background_population=CD19neg,xCol="BV786 CD14-A", yCol="APC CD16-A",markersize=0.1, scale='bilog', T=500,filePlot=fileName)
        # plotting the Inter_Mono
        Inter_Mono = ag.AGgate(indices_dict["Inter_Mono"],CD19neg,"BV786 CD14-A","APC CD16-A","CM")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-Inter_Mono.png"
        ag.backGate(fcs,population=Inter_Mono,background_population=CD19neg,xCol="BV786 CD14-A", yCol="APC CD16-A",markersize=0.2, scale='bilog', T=500,filePlot=fileName)
        # plotting the non_CM
        non_CM = ag.AGgate(indices_dict["non_CM"],CD19neg,"BV786 CD14-A","APC CD16-A","non_CM")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-non_CM.png"
        ag.backGate(fcs,population=non_CM,background_population=CD19neg,xCol="BV786 CD14-A", yCol="APC CD16-A",markersize=0.2, scale='bilog', T=500,filePlot=fileName)
        # plotting the CD14-
        CD14neg = ag.AGgate(indices_dict["CD14neg"],CD19neg,"BV786 CD14-A","APC CD16-A","non_CM")
        fileName="/home/wenxiaren/Part_3_ML/Prediction_Evaluation/Multiple_Classes_CD3neg_monocytes/plots/"+sampleName+"-CD14neg.png"
        ag.backGate(fcs,population=CD14neg,background_population=CD19neg,xCol="BV786 CD14-A", yCol="APC CD16-A",markersize=0.2, scale='bilog', T=500,filePlot=fileName)

    